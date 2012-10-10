
#from numpy import *
from numpy import *
from matplotlib.pyplot import *
import glob, os
from mayavi import mlab
import scitools.std as sci
import glob,os

for filename in glob.glob('wtmp*.png'):
    os.remove(filename)

I_f = lambda x,y: ones((len(x), len(y)), float)*1.2
V_f = lambda x,y: zeros((len(x)+2, len(y)+2), float)
q_f = lambda x,y: ones((len(x)+2, len(y)+2), float)*c
f_f = lambda x,y: zeros((len(x)+2, len(y)+2), float)

def gauss(x,y):
    sigma_x = 0.5
    sigma_y = 0.5
    x_0 = 0
    y_0 = 0
    return exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))

def standing_source(x,y,t):
    return cos(my*y*pi/Ly)*cos(mx*x*pi/Lx)*(exp(-b*t)*cos(w*t)*(q*pi**2*((mx/Lx)**2 + (my/Ly)**2) + exp(-b*t)*b*w*sin(w*t)))


def standing_V(x,y):
    return -b*standing_exact(x,y,0)
    #return -b*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
    
def standing_exact(x,y,t=0):
    return exp(-b*t)*cos(w*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)


c = 0.8
b = 0.2
my = 0
mx = 5.0
xstop = 10.0
xstart = 0
ystop = 10.0# + xstop/4
ystart = 0
tstop = 10.0
Ly = -ystart+ystop
Lx = xstop-xstart
w = 2.0
Nx = 30
Ny = Nx
#Nt = 30

#boundary conditions du/dn = 0, meaning u_Nx = u_Nx+1
x = linspace(xstart, xstop, Nx)
y = linspace(ystart,ystop,Ny)
#t = linspace(0,tstop, Nt)
#dt = t[1] - t[0]
dx = x[1] - x[0]
dy = y[1] - y[0]
if dx>dy:
    dt = dx/sqrt(2.0)
else:
    dt = dy/sqrt(2.0)
#dt = 0.5
Nt = int(tstop/dt)
t = linspace(0, tstop,Nt)
X1,Y1 = meshgrid(linspace(xstart,xstop,Nx),linspace(ystart,ystop,Ny))
X,Y = meshgrid(linspace(xstart,xstop,Nx+2),linspace(ystart,ystop,Ny+2))
q = q_f(x,y)
#V = V_f(x,y)
V = standing_V(X,Y)
I = standing_exact(x,y,0)
f0 = standing_source(X,Y,0)
u = zeros((Nx+2,Ny+2), float)
up = u.copy()
upp = up.copy()


#initial conditions
#u^-1ij = 2*dt*Vij + u^1ij, u^0 = I
C_x = dt/dx
C_y = dt/dy
up[1:-1,1:-1] = I.copy()

"""
a = up[1:-1]
for j in xrange(Ny):
    a[j][1:-1] = I[j].copy()

"""
up[0,:] = up[1,:].copy()
up[:,0] = up[:,1].copy()
up[-1,:] = up[-2,:].copy()
up[:,-1] = up[:,-2].copy()



#making u^1
for i in xrange(1,Nx+1):
    for j in xrange(1,Ny+1):
        x_para = ((q[i][j] + q[i+1][j])*(up[i+1][j] - up[i][j]) - (q[i-1][j] + q[i][j])*(up[i][j] - up[i-1][j]))
        y_para = ((q[i][j] + q[i][j+1])*(up[i][j+1] - up[i][j]) - (q[i][j-1] + q[i][j])*(up[i][j] - up[i][j-1]))
        rest = f0[i][j] + 4*up[i][j] + 2*dt*V[i][j]*(b*dt-2)
        u[i][j] = 1.0/(4.0)*(C_x**2*x_para + C_y**2*y_para + rest)

u[0,:] = u[1,:].copy()
u[:,0] = u[:,1].copy()
u[-1,:] = u[-2,:].copy()
u[:,-1] = u[:,-2].copy()

#making u^-1
upp = 2*dt*V + u



counter =0
count = 0
"""
test = standing_exact(X,Y,0)
#s = mlab.mesh(X,Y,test)
filename = "exact_%d.gif"%Nt
for k in xrange(Nt):
    #mlab.figure()
    mlab.mesh(X,Y,test)
    test = standing_exact(X,Y,t[k])
    #s.mlab_source.scalars = test
    mlab.savefig("wtmp%04d.png"%k)
    mlab.clf()
"""
#vectorized:
filename = "num_%d.gif"%Nt
for k in xrange(Nt):
    x_para = (q[1:-1,1:-1] + q[2:,1:-1])*(up[2:,1:-1] - up[1:-1,1:-1]) - (q[:-2, 1:-1] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[:-2,1:-1])
    y_para = (q[1:-1,1:-1] + q[1:-1,2:])*(up[1:-1,2:] - up[1:-1,1:-1]) - (q[1:-1,:-2] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[1:-1,:-2])
    f0 = standing_source(X,Y,t[k])
    #print C_x, C_y
    rest = f0[1:-1,1:-1] + 4*up[1:-1,1:-1] + upp[1:-1,1:-1]*(b*dt-2)
    #print shape(x_para), shape(y_para), shape(rest), shape(u[1,-1])
    u[1:-1,1:-1] = 1.0/(2+b*dt)*(C_x**3*x_para + C_y**2*y_para + rest)

    u[0,:] = u[1,:]
    u[:,0] = u[:,1]
    u[-1,:] = u[-2,:]
    u[:,-1] = u[:,-2]
    
    if k%5 == 0:
        5+5
        #f = mlab.figure()
        s = mlab.mesh(X,Y,u, color=(0.0,0.75,1.0))
        #s.mlab_source.scalars = u
        #h = mesh(X,Y,u)
        mlab.savefig("wtmp%04d.png" %k)
    mlab.clf()
    #mlab.draw()
    #print k
    upp = up.copy()
    up = u.copy()


#scalar-way

"""
#f = figure()
s = mlab.mesh(X,Y,u)
for k in xrange(Nt):
    for i in xrange(1,Nx):
        for j in xrange(1,Ny):
             x_para = ((q[i][j] + q[i+1][j])*(up[i+1][j] - up[i][j]) - (q[i-1][j] + q[i][j])*(up[i][j] - up[i-1][j]))
             y_para = ((q[i][j] + q[i][j+1])*(up[i][j+1] - up[i][j]) - (q[i][j-1] + q[i][j])*(up[i][j] - up[i][j-1]))
             rest = f0[i][j] + 4*up[i][j] + upp[i][j]*(b*dt-2)
             u[i][j] = 1.0/(2+b*dt)*(C_x**2*x_para + C_y**2*y_para + rest)

    u[0,:] = u[1,:]
    u[:,0] = u[:,1]
    u[-1,:] = u[-2,:]
    u[:,-1] = u[:,-2]
    
    if k%5 == 0:
        #f = mlab.figure()
        #s = mlab.mesh(X,Y,u)
        s.mlab_source.scalars = u
        mlab.savefig("wtmp%04d.png" %k)
    #mlab.clf()
    #mlab.draw()
    #print k
    upp = up.copy()
    up = u.copy()

"""

#mlab.show()
sci.movie("wtmp*.png",encoder='convert', fps=2, output_file=filename)

"""
for filename in glob.glob('wtmp*.png'):
    os.remove(filename)
"""

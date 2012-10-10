rom numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import glob, os
from mayavi import mlab
from scitools.std import movie

for filename in glob.glob('wtmp*.png'):
    os.remove(filename)

I_f = lambda x,y: ones((len(x), len(y)), float)*1.2
V_f = lambda x,y: zeros((len(x)+2, len(y)+2), float)
q_f = lambda x,y: ones((len(x)+2, len(y)+2), float)*c
f_f = lambda x,y: zeros((len(x), len(y)), float)
def gauss(x,y):
    sigma_x = 0.5
    sigma_y = 0.5
    x_0 = 0
    y_0 = 0
    return exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))

c = 0.8
b = 2.0
Nx = 50
Ny = Nx
#Nt = 30
xstop = 5.0
xstart = -5.0
ystop = 5.0
ystart = -5.0
tstop = 20.0
#boundary conditions du/dn = 0, meaning u_Nx = u_Nx+1
x = linspace(xstart, xstop, Nx)
y = linspace(ystart,ystop,Ny)
#t = linspace(0,tstop, Nt)
#dt = t[1] - t[0]
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = dx/sqrt(2.0)
#dt = 0.5
Nt = int(tstop/dt)
t = linspace(0, tstop,Nt)
X1,Y1 = meshgrid(linspace(xstart,xstop,Nx),linspace(ystart,ystop,Ny))
q = q_f(x,y)
V = V_f(x,y)
I = gauss(X1,Y1)
f0 = f_f(x,y)
u = zeros((Nx+2,Ny+2), float)
up = u.copy()
upp = up.copy()


#initial conditions
#u^-1ij = 2*dt*Vij + u^1ij, u^0 = I
C_x = dt/dx
C_y = dt/dy
a = up[1:-1]

for j in xrange(Ny):
    a[j][1:-1] = I[j].copy()

up[0,:] = up[1,:].copy()
up[:,0] = up[:,1].copy()
up[-1,:] = up[-2,:].copy()
up[:,-1] = up[:,-2].copy()


#making u^1
for i in xrange(1,Nx):
    for j in xrange(1,Ny):
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



X,Y = meshgrid(linspace(xstart,xstop,Nx+2),linspace(ystart,ystop,Ny+2))

counter =0
count = 0

#f = figure()
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
    
    if k%10 == 0:
        
        counter+=1
        #f = mlab.figure()
        s = mlab.mesh(X,Y,u)
        mlab.savefig("wtmp%04d.png" %k)
    mlab.clf()
    #mlab.draw()
    #print k
    upp = up.copy()
    up = u.copy()
print counter


#mlab.show()
movie("wtmp*.png")


for filename in glob.glob('wtmp*.png'):
    os.remove(filename)


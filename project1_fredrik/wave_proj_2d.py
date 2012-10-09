"""
Solves the 2.dim wave-equation with dampening. Optional commandline arguments are number of timesteps (T),length of area in x and y  direction (Lx,Ly), number of meshpoints in x and y direction (Nx,Ny), and the size of the timestep (dt). The user cal also specify if a vectorized (default) or scalar solver should be used by setting "-s" for "scalar", and if a simple wave solver ithout damping and underlying geography should be used by setting "-b" for "basic". A program call might look like this:

user@computer$: python wave_proj_2d.py -Lx 8 -T 350 -Nx 50 -b

"""
from numpy import *
from math import *
import sys, os, time, math,argparse
#import matplotlib.pyplot as plt
#from scitools.std import *
from mayavi import mlab
#--------Initialization---------------
parser = argparse.ArgumentParser()
parser.add_argument("-s",action="store_true") 
parser.add_argument("-b",action="store_true")
parser.add_argument("-Lx", type = int, dest="Lx", help="Size of area in x direction")
parser.add_argument("-Ly", type = int, dest="Ly", help="Size of area in y direction")
parser.add_argument("-T", type = int, dest="T", help="Number of timesteps")
parser.add_argument("-Nx", type = int, dest="Nx", help="Number of gridpoints in x direction")
parser.add_argument("-Ny", type = int, dest="Ny", help="Number of gridpoints i y direction")
parser.add_argument("-dt", type = float, dest="dt", help="timestep")
args = parser.parse_args()

Lx = args.Lx if args.Lx != None else 5
Ly = args.Ly if args.Ly != None else Lx
T = args.T if args.T != None else 100
Nx = args.Nx if args.Nx != None else 25
Ny = args.Ny if args.Ny != None else Nx

dx = Lx/float(Nx+1); dt = 1./float(T+1); dy = Ly/float(Ny+1);
if args.dt <=0 or args.dt==None:
	dt = dx/sqrt(2)
else:
	dt = args.dt
sigma_x = 0.5; sigma_y = 0.5;
X,Y = meshgrid(linspace(0,Ly,Nx/dx),linspace(0,Ly,Ny/dy))
n = len(X[0]);
h = zeros((n,n));
x_0 = Lx/2.0; y_0 = Ly/2.0;
u0 = zeros((n,n)); u1 = zeros((n,n));
uny = zeros((n,n));
f = zeros((n,n));
q = ones((n,n));
b = 0.3;	# dampening coefficient

#--------Initial conditions------------
def initial(x,y):
	"""
	Returns the initial shape of the wave
	"""
	return math.exp(-0.5*(((x-x_0)/(sigma_x))**2+((y-y_0)/(sigma_y))**2))

def geography(x,y):
	a = 0.3; c = -0.4;
	return a*x + c*y
def func(f,t):
	"""
	Returns the source term in a point at a given time
	"""
	g = -0.3; k = 13; p = 3.8; m = 7.4
	loc_sin = math.sin
	for i in xrange(n):
		for j in xrange(n):
			f[i,j] = g*loc_sin(k*X[i,j] + p*Y[i,j] - m*t)
	return f

for i in xrange(1,n-1):
	for j in xrange(1,n-1):
		h[i,j] = initial(X[i,j],Y[i,j]);
#q[:,:] = geography(X[:,:],Y[:,:]);
q *= 0.8; 
h[0,1:-1] = h[1,1:-1] 
h[1:-1,0] = h[1:-1,1] 
h[1:-1,n-1] = h[1:-1,n-2]
h[-1,1:-1] = h[-2,1:-1] 

u0 =copy(h);
A = 0;
B = 0;
C = 0;
#save some unnececary FLOPS
scale = ((dt*dt)/(1+0.5*b*dt))	
Dx = (1./(2*dx*dx))		
Dy = (1./(2*dy*dy))		
v = (2/(1+0.5*b*dt))		
r = ((1-0.5*b*dt)/(1+0.5*b*dt))	
Cx2 = (0.8*dt/dx)**2
Cy2 = (0.8*dt/dy)**2
dt2 = dt*dt
#--------Working loop---------------------
# u1 = up
# u0 = upp
if args.s and not args.b:
	#scalar version
	for j in xrange(1,n-1):
		for k in xrange(1,n-1):
			#Makes the first timestep 
			
			A= dt2*Dx*((u0[j+1,k] - u0[j,k])*(q[j+1,k] + q[j,k]) -(u0[j,k] - u0[j-1,k])*(q[j,k] + q[j-1,k]));
			B = dt2*Dy*((u0[j,k+1] - u0[j,k])*(q[j,k+1] + q[j,k]) - (u0[j,k] - u0[j,k-1])*(q[j,k] + q[j,k-1]));
			C = u0[j,k]
			u1[j,k] = (0.5*A + 0.5*B + C); 
	u1 = u0.copy()
	for i in xrange(T):
		for j in xrange(1,n-1):
			for k in xrange(1,n-1):
				A = Dx*( (u1[j+1,k]-u1[j,k])*(q[j+1,k]+q[j,k])-(u1[j,k]-u1[j-1,k])*(q[j,k]+q[j-1,k]) );
				B = Dy*( (u1[j,k+1]-u1[j,k])*(q[j,k+1]+q[j,k])-(u1[j,k]-u1[j,k-1])*(q[j,k]+q[j,k-1]) );
				C = v*u1[j,k] - r*u0[j,k];
				uny[j,k] = 0.5*scale*A + 0.5*scale*B + C;
		print i
		uny[0,1:-1] = uny[1,1:-1] 
		uny[1:-1,0] = uny[1:-1,1] 
		uny[1:-1,n-1] = uny[1:-1,n-2]
		uny[-1,1:-1] = uny[-2,1:-1] 

elif args.b and not args.s:
	u1 = u0.copy()
	for i in xrange(T):
		print i
		uny[1:-1,1:-1] = Cx2*(u1[:-2,1:-1] - 2*u1[1:-1,1:-1] + u1[2:,1:-1]) + \
		Cy2*(u1[1:-1,:-2] - 2*u1[1:-1,1:-1] + u1[1:-1,2:]) +\
		2*u1[1:-1,1:-1] - u0[1:-1,1:-1] 
		uny[0,1:-1] = uny[1,1:-1] 
		uny[1:-1,0] = uny[1:-1,1] 
		uny[1:-1,n-1] = uny[1:-1,n-2]
		uny[-1,1:-1] = uny[-2,1:-1] 
elif args.b and args.s:
	u1 = u0.copy()
	for i in xrange(T):	
		for j in xrange(1,n-1):
			for k in xrange(1,n-1):
				A = Cx2*(u1[j+1,k] - 2*u1[j,k] + u1[j-1,k]);
				B = Cy2*(u1[j,k+1] - 2*u1[j,k] + u1[j,k-1]);
				C = 2*u1[j,k] - u0[j,k];
				uny[j,k] = A + B + C;
			
		uny[0,1:-1] = uny[1,1:-1] 
		uny[1:-1,0] = uny[1:-1,1] 
		uny[1:-1,n-1] = uny[1:-1,n-2]
		uny[-1,1:-1] = uny[-2,1:-1] 
	
		u0 = u1.copy();
		u1 = uny.copy();
		print i
	#s = mlab.mesh(X, Y, u1)
	#mlab.show()
	#time.sleep(1)
else:
	#vectorized (default) version
	u1[1:-1,1:-1] = (0.5*dt2*Dx*((u0[2:,1:-1]-u0[1:-1,1:-1])*(q[2:,1:-1] + q[1:-1,1:-1])-(u0[1:-1,1:-1]-u0[:-2,1:-1])*(q[1:-1,1:-1])+q[:-2,1:-1]) \
	+0.5*dt2*Dy*((u0[1:-1,2:]-u0[1:-1,1:-1])*(q[1:-1,2:]+q[1:-1,1:-1])-(u0[1:-1,1:-1]-u0[1:-1,:-2])*(q[1:-1,1:-1]+q[1:-1,:-2])) \
	+u0[1:-1,1:-1])
	u1[0,1:-1] = u1[1,1:-1] 
	u1[1:-1,0] = u1[1:-1,1] 
	u1[1:-1,n-1] = u1[1:-1,n-2]
	u1[-1,1:-1] = u1[-2,1:-1] 
	u1 = u0.copy()
	for i in xrange(T):
		#f = f(f,i) # updates the source term
		uny[1:-1,1:-1] = scale*Dx*((u1[2:,1:-1]-u1[1:-1,1:-1])*(q[2:,1:-1]+q[1:-1,1:-1])-(u1[1:-1,1:-1]-u1[:-2,1:-1])*(q[1:-1,1:-1]+q[:-2,1:-1])) \
		+scale*Dy*((u1[1:-1,2:]-u1[1:-1,1:-1])*(q[1:-1,2:]+q[1:-1,1:-1])-(u1[1:-1,1:-1]-u1[1:-1,:-2])*(q[1:-1,1:-1]+q[1:-1,:-2])) \
		+ v*u1[1:-1,1:-1] - r*u0[1:-1,1:-1]
		print(i)
		uny[0,1:-1] = uny[1,1:-1] 
		uny[1:-1,0] = uny[1:-1,1] 
		uny[1:-1,n-1] = uny[1:-1,n-2]
		uny[-1,1:-1] = uny[-2,1:-1] 		
		u0 = u1.copy();
		u1 = uny.copy();
#print uny
s = mlab.mesh(X, Y, u1)
mlab.show()
#print u1

'''
This part is wrong for some reason...
u1[1:-1,1:-1] = (Dx*((u0[2:,1:-1]-u0[1:-1,1:-1])*(q[2:,1:-1] + q[1:-1,1:-1])-(u0[1:-1,1:-1]-u0[:-2,1:-1])*(q[1:-1,1:-1])+q[:-2,1:-1]) \
	+Dy*((u0[1:-1,2:]-u0[1:-1,1:-1])*(q[1:-1,2:]+q[1:-1,1:-1])-(u0[1:-1,1:-1]-u0[1:-1,:-2])*(q[1:-1,1:-1]+q[1:-1,:-2])) \
	+v*u0[1:-1,1:-1])
u1 /= (1+((2-b*dt)/(2+b*dt)))
#-----Langtangen code

u1[1:-1,1:-1] = u0[1:-1,1:-1] +  \
0.5*Cx2*(u0[:-2,1:-1] - 2*u0[1:-1,1:-1] + u0[2:,1:-1]) +0.5*Cy2*(u0[1:-1,:-2] - 2*u0[1:-1,1:-1] + u0[1:-1,2:])
u1[0,1:-1] = u1[1,1:-1] 
u1[1:-1,0] = u1[1:-1,1] 
u1[1:-1,n-1] = u1[1:-1,n-2]
u1[-1,1:-1] = u1[-2,1:-1] 
'''
'''
for j in xrange(1,n-1):
	for k in xrange(1,n-1):
		A = dt2*Dx*( (u1[j+1,k]-u1[j,k])*(q[j+1,k]+q[j,k])-(u1[j,k]-u1[j-1,k])*(q[j,k]+q[j-1,k]) );
		B = dt2*Dy*( (u1[j,k+1]-u1[j,k])*(q[j,k+1]+q[j,k])-(u1[j,k]-u1[j,k-1])*(q[j,k]+q[j,k-1]) );
		C = v*u1[j,k] - r*u0[j,k];
		uny[j,k] = A + B + C;
'''	

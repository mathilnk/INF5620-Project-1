"""
Solves the 2.dim wave-equation with dampening. Commandline arguments are 
"""
from numpy import *
from math import *
import sys, os, time
#import matplotlib.pyplot as plt
from scitools.std import *
from mayavi import mlab

#--------Initialization---------------
#C = float(sys.argv[1]); # The courrant number
T = int(sys.argv[1]);
Nx = int(sys.argv[2]);
Ny = Nx;
Lx = float(sys.argv[3])
Ly = Lx;
dx = Lx/float(Nx+1); dt = 1./float(T+1); dy = Ly/float(Ny+1);
sigma_x = 0.5; sigma_y = 0.5;
X,Y = meshgrid(linspace(0,Ly,Nx/dx),linspace(0,Ly,Ny/dy))
n = len(X[0]);
h = zeros((n,n));
#v = C; # v is the wave velocity C = v*dt/dx; dt = dx
x_0 = Lx/2.0; y_0 = Ly/2.0;
u0 = zeros((n,n)); u1 = zeros((n,n));
uny = zeros((n,n));
f = zeros((n,n));
q = ones((n,n));
b = 0.35;	# dampening coefficient
#--------Initial conditions------------
def initial(x,y):
	"""
	Returns the initial shape of the wave
	"""
	return exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))

def geography(x,y):
	a = 0.3; c = -0.4;
	return a*x + c*y

for i in range(n):
	for j in range(n):
		h[i][j] = initial(X[i][j],Y[i][j]);
		#q[i][j] = geography(X[i][j],Y[i][j]);
		q[i,j] *= 0.8; 




u0 =copy(h);
#--------Velocities at t=0---------------


for i in range(1,n-1):
	for j in range(1,n-1):
		u1[i][j] = h[i][j] #Makes the first timestep with forward euler
	u1[0][i] = u1[1][i] # 0;
	u1[n-1][i] = u1[n-2][i] # 0;
	u1[i][0] = u1[i][1] # 0;
	u1[i][n-1] = u1[i][n-2]

#--------Working loop---------------------
# u1 = up
# u0 = upp
A = 0;
B = 0;
C = 0;
scale = ((dt*dt)/(1+0.5*b*dt))	#save some unnececary FLOPS 
Dx = (1./(2*dx*dx))		#save some unnececary FLOPS
Dy = (1./(2*dy*dy))		#save some unnececary FLOPS
v = (2/(1+0.5*b*dt))		#save some unnececary FLOPS
r = ((2-b*dt)/(2+b*dt))		#save some unnececary FLOPS

for i in range(T):		
	for j in range(1,n-1):
		for k in range(1,n-1):
			A= Dx*((u1[j+1][k] - u1[j][k])*(q[j+1][k] + q[j][k]) -(u1[j][k] - u1[j-1][k])*(q[j][k] + q[j-1][k]));
			B = Dy*((u1[j][k+1] - u1[j][k])*(q[j][k+1] + q[j][k]) - (u1[j][k] - u1[j][k-1])*(q[j][k] + q[j][k-1]));
			C =  scale*f[j][k] + v*u1[j][k] - r*u0[j][k];
			uny[j][k] = scale*A + scale*B + C;		
		uny[0][j] = 0#uny[1][j] # 0;
		uny[n-1][j] = 0#uny[n-2][j] # 0;
		uny[j][0] = 0#uny[j][1] # 0;
		uny[j][n-1] = 0#uny[j][n-2] # 0;

	u0 = copy(u1);
	u1 = copy(uny);
	print i
	#s = mlab.mesh(X, Y, u1)
	#mlab.show()
	#time.sleep()
s = mlab.mesh(X, Y, u1)
mlab.show()
#print u1


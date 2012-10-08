from numpy import *
from matplotlib.pyplot import *
import time
import glob, os
for filename in glob.glob('tmp*.png'):
    os.remove(filename)

def gauss(x,mu=0,sigma=1.0):
    return 1/(sqrt(2*pi)*sigma)*exp(-0.5*((x-mu)/sigma)**2)

def V_0(x):
    return 0.01*x

def source_0(t):
    return x*0

def I_construct(x):
    return x*L-x**2

def V_construct(x):
    return 0.5*(L*x-x**2)

def f_construct(x,t):
    return 2*(1+0.5*t)*c**2

cour = 0.7
c = cour
#dt =dx
dt = 0.1
dx = 0.1


L = 10.0
T = 40.0
Nt = T/dt
Nx = int(L/dx)


x = linspace(-L/2,L/2,Nx)
t = linspace(0,T,Nt)

u = zeros(Nx)
cour = c*dt/dx


V = V_construct(x)
f = f_construct(x,t)


up = I_construct(x)
upp = up - V*2*dt


ion()
lines = plot(x,u)
axis([x[0], x[-1],-100, 1000])
xlabel('x')
ylabel('u')
#show()
counter = 0



for n in xrange(len(t)):
    u[1:-1] = -upp[1:-1] + 2*(1-cour**2)*up[1:-1] + cour**2*(up[2:]+ up[:-2]) + dt**2*f[n]

    
    lines[0].set_ydata(u)
    draw()
    #savefig('tmp_%40d.png'%counter)
    counter+=1
    upp = up.copy()
    up = u.copy()
    #time.sleep(2)
#movie('tmp_*.png')

    


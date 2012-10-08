from numpy import *
from matplotlib.pyplot import *
import time
import glob, os
for filename in glob.glob('tmp*.png'):
    os.remove(filename)

def gauss(x,mu=0,sigma=1.0):
    return 1/(sqrt(2*pi)*sigma)*exp(-0.5*((x-mu)/sigma)**2)


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
up = gauss(x)
upp = up.copy()

dt = abs(t[0] - t[1])
dx = abs(x[0] - x[1])

cour = c*dt/dx


ion()
lines = plot(x,u)
axis([x[0], x[-1],-max(up), max(up)])
xlabel('x')
ylabel('u')
#show()
counter = 0



for n in xrange(len(t)):
    """
    up = u.copy()
    for j in xrange(1,Nx-1):
        u[j] = -upp[j] + 2*(1-cour**2)*up[j] + cour**2*(up[j+1] +  up[j-1])
    upp = up.copy()
    up = u.copy()
    """
    u[1:-1] = -upp[1:-1] + 2*(1-cour**2)*up[1:-1] + cour**2*(up[2:]+ up[:-2])
    lines[0].set_ydata(u)
    draw()
    savefig('tmp_%40d.png'%counter)
    counter+=1
    upp = up.copy()
    up = u.copy()
    #time.sleep(2)
#movie('tmp_*.png')

    


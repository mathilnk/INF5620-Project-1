from numpy import *
from matplotlib.pyplot import *
import time
import glob, os

class wave_1d:


    def __init__(self, I, V, f, L, T, cour, exact):
        self.I = I
        self.V = V
        self.f = f
        self.exact = exact
        self.L = L
        self.T = T
        self.cour = cour
        self.dt = 0.1
        self.dx = 0.1
        self.c = cour
        
        self.t = linspace(0,T,T/self.dt)
        self.x = linspace(0,L,L/self.dx)

        self.u = zeros(len(self.x))



    def solve(self):
        V = self.V(self.x)
        up = self.I(self.x)
        upp = up - V*2*self.dt

        ion()
        lines = plot(self.x,self.u)
        axis([self.x[0], self.x[-1], -100, 400])
        counter = 0
        for n in xrange(len(self.t)):
            self.u[1:-1] = -upp[1:-1] + 2*(1-self.cour**2)*up[1:-1] + self.cour**2*(up[2:]+ up[:-2]) + self.dt**2*self.f(self.x, self.t)[n][1:-1]
            lines[0].set_ydata(self.u)
            draw()
            #savefig('tmp_%40d.png'%counter)
            counter+=1
            upp = up.copy()
            up = self.u.copy()

    def plot_exact(self):
        ion()
        lines = plot(self.x,self.exact(self.x, t=0))
        axis([self.x[0], self.x[-1], -100, 1000])
        counter = 0
        for n in xrange(len(self.t)):
            lines[0].set_ydata(self.exact(self.x,self.t[n]))
            draw()
            #savefig('tmp_%40d.png'%counter)
            counter+=1
        
        
            
        

    

    


    
    

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
    f = zeros((len(t), len(x)))
    for i in xrange(len(t)):
        f[i] = 2*(1+0.5*t[i])*c**2
    return f

def ex(x,t):
    return x*(L-x)*(1+ 0.5*t)

L = 10.0
T = 10.0
cour = 0.7
c = cour
S = wave_1d(I_construct, V_construct, f_construct, L, T, cour, ex)

#S.plot_exact()
S.solve()




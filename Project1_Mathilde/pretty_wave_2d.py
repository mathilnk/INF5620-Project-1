import scitools.std as sci
import glob,os,sys,argparse
from numpy import *
from mayavi import mlab

V_f_glob = lambda x,y: x*0
q_f_glob = lambda x,y: ones(shape(x))*sqrt(2)
f_f_glob = lambda x,y,t: x*0



class wave_2d:

    def __init__(self, Lx, Ly, T, Nx, Ny, dt, b = 0,I_f=None, V_f=V_f_glob, q_f=q_f_glob, f_f=f_f_glob,exact=None, gauss=False, standing=False):
        """
        I_f is the initial state-function, x,y
        V_f is the initial velocity-function, x,y
        q_f is the velocity of the wave at position x,y
        f_f is the source function , x,y,t
        gauss = True => I_f is set to a gauss function
        standing = True => I_f, V_f, q_f, f_f is set to fit a standing wave solution.
        """
        self.exact = exact
        if I_f==None and gauss==False and standing==False:
            print "You failed to give wave properties"
            sys.exit(1)
        elif(gauss):
            self.I_f = self.gauss_f
            self.f_f = f_f
            self.V_f = V_f
            self.q_f = q_f
        elif(standing):
            self.I_f = self.standing_exact
            self.f_f = self.standing_source
            self.V_f = self.standing_V
            self.q_f = q_f
            self.exact = self.standing_exact
            self.w = 2.0
            self.mx = 5
            self.my = 5
            self.standing = standing
        else:
            self.I_f = I_f

        self.Lx = Lx
        self.Ly = Ly
        self.T = T
        self.Nx = Nx
        self.Ny = Ny
        self.dt = dt
        self.b = b
        self.f_f = f_f


        self.x = linspace(0,Lx,Nx)
        self.y = linspace(0,Ly,Ny)
        self.dx = self.x[1]-self.x[0]
        self.dy = self.y[1]-self.y[0]
        self.Nt = int(T/dt)
        self.t = linspace(0,T,self.Nt)

        self.X1,self.Y1 = meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny))

        self.X,self.Y = meshgrid(linspace(0,Lx,Nx+2),linspace(0,Ly,Ny+2))
        X = self.X
        Y = self.Y
        self.V = V_f(X,Y)
        self.I = self.I_f(X,Y)
        self.q = q_f(X,Y)
        self.f = f_f(X,Y,0)
        #print "hallo"
        #self.q = q_f(X,Y)
        self.C_y = dt/self.dy
        self.C_x = dt/self.dx

    def gauss_f(self,x,y):
        sigma_x = 0.5
        sigma_y = 0.5
        x_0 = 0
        y_0 = 0
        return exp(-(((x-x_0)/(2*sigma_x))**2+((y-y_0)/(2*sigma_y))**2))

    def standing_source(self,x,y,t):#,my=0.5, mx=0.5):
        w = self.w
        mx = self.mx
        my = self.my
        return cos(my*y*pi/self.Ly)*cos(mx*x*pi/self.Lx)*(exp(-self.b*t)*cos(w*t)*(self.q*pi**2*((mx/self.Lx)**2 + (my/self.Ly)**2) + exp(-self.b*t)*self.b*w*sin(w*t)))
    
    
    def standing_V(self,x,y):
        return -self.b*self.standing_exact(x,y,0)
    #return -b*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
    
    def standing_exact(self,x,y,t=0):
        w = self.w
        mx = self.mx
        my = self.my
        return exp(-self.b*t)*cos(w*t)*cos(mx*x*pi/self.Lx)*cos(my*y*pi/self.Ly)

    def make_exact(self):
        if(self.exact!=None):
            sol = self.exact(self.X,self.Y,0)
            for k in xrange(self.Nt):
                mlab.mesh(self.X,self.Y,sol, color=(0.0, 0.3, 0.6))
                sol = self.exact(self.X, self.Y, self.t[k])
                mlab.savefig("wtmp%04d.png" %k)
                mlab.clf()
            filename = "exact_dt%2.1f.gif" %self.dt
            sci.movie("wtmp*.png",encoder='convert', fps=5, output_file=filename)
            
        
    def solve_num(self):
        Nx = self.Nx
        Ny = self.Ny
        Nt = self.Nt
        I = self.I
        q = self.q
        f = self.f
        dt = self.dt
        V = self.V
        b = self.b
        C_x = self.C_x
        C_y = self.C_y
        f_f = self.f_f
        X = self.X
        Y = self.Y
        t = self.t
        
        u = zeros((Nx+2,Ny+2), float)
        up = u.copy()
        upp = up.copy()
        
        up[1:-1,1:-1] = I[1:-1,1:-1].copy()
        
        
        up[0,:] = up[1,:].copy()
        up[:,0] = up[:,1].copy()
        up[-1,:] = up[-2,:].copy()
        up[:,-1] = up[:,-2].copy()

        
        


        #making u^1
        for i in xrange(1,Nx+1):
            for j in xrange(1,Ny+1):
                x_para = ((q[i][j] + q[i+1][j])*(up[i+1][j] - up[i][j]) - (q[i-1][j] + q[i][j])*(up[i][j] - up[i-1][j]))
                y_para = ((q[i][j] + q[i][j+1])*(up[i][j+1] - up[i][j]) - (q[i][j-1] + q[i][j])*(up[i][j] - up[i][j-1]))
                rest = f[i][j] + 4*up[i][j] + 2*dt*V[i][j]*(b*dt-2)
                u[i][j] = 1.0/(4.0)*(C_x**2*x_para + C_y**2*y_para + rest)

        u[0,:] = u[1,:].copy()
        u[:,0] = u[:,1].copy()
        u[-1,:] = u[-2,:].copy()
        u[:,-1] = u[:,-2].copy()
        
        #making u^-1
        upp = 2*dt*V + u
                


        #vectorized:
        filename_2 = "num_dt%2.1f.gif"%dt
        if self.standing:
            filename_2 = "standing_wave_dt%2.1f.gif"%dt
        for filename in glob.glob('wtmp*.png'):
            os.remove(filename)
        
        for k in xrange(Nt):
            x_para = (q[1:-1,1:-1] + q[2:,1:-1])*(up[2:,1:-1] - up[1:-1,1:-1]) - (q[:-2, 1:-1] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[:-2,1:-1])
            y_para = (q[1:-1,1:-1] + q[1:-1,2:])*(up[1:-1,2:] - up[1:-1,1:-1]) - (q[1:-1,:-2] + q[1:-1,1:-1])*(up[1:-1,1:-1] - up[1:-1,:-2])
            f = f_f(X,Y,t[k])
            #print C_x, C_y
            rest = f[1:-1,1:-1] + 4*up[1:-1,1:-1] + upp[1:-1,1:-1]*(b*dt-2)
            
            
            u[1:-1,1:-1] = 1.0/(2+b*dt)*(C_x**3*x_para + C_y**2*y_para + rest)
            
            u[0,:] = u[1,:]
            u[:,0] = u[:,1]
            u[-1,:] = u[-2,:]
            u[:,-1] = u[:,-2]
    
            if k%3 == 0:
                
                s = mlab.mesh(X,Y,u, color=(0.0,0.75,1.0))
                
                
                mlab.savefig("wtmp%04d.png" %k)
                mlab.clf()
            upp = up.copy()
            up = u.copy()

        

        sci.movie("wtmp*.png",encoder='convert', fps=2, output_file=filename_2)
 
    
        
        
    
class Standing_wave(wave_2d):
    def __init__(self,w=2.0,mx=5,my=5,dt=0.2,Lx=10.0,Ly=10.0,T=80.0,Nx=30,Ny=30,b=0.0):
        self.w=w
        self.mx=mx
        self.my=my
        print "hei"
        wave_2d.__init__(self,Lx,Ly,T,Nx,Ny,dt,b,I_f=self.standing_exact,V_f=self.standing_V,f_f=self.standing_source, exact=self.standing_exact)
       
   


    def standing_source(self,x,y,t):#,my=0.5, mx=0.5):
        w = self.w
        mx = self.mx
        my = self.my
        return cos(my*y*pi/self.Ly)*cos(mx*x*pi/self.Lx)*(exp(-self.b*t)*cos(w*t)*(self.q*pi**2*((mx/self.Lx)**2 + (my/self.Ly)**2) + exp(-self.b*t)*self.b*w*sin(w*t)))
    
    
    def standing_V(self,x,y):
        return -self.b*self.standing_exact(x,y,0)
    #return -b*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
    
    def standing_exact(self,x,y,t=0):
        w = self.w
        mx = self.mx
        my = self.my
        return exp(-self.b*t)*cos(w*t)*cos(mx*x*pi/self.Lx)*cos(my*y*pi/self.Ly)

    def make_exact(self):
        if(self.exact!=None):
            sol = self.exact(self.X,self.Y,0)
            for k in xrange(self.Nt):
                mlab.mesh(self.X,self.Y,sol, color=(0.0, 0.3, 0.6))
                sol = self.exact(self.X, self.Y, self.t[k])
                mlab.savefig("wtmp%04d.png" %k)
                mlab.clf()
            filename = "exact_dt%2.1f.gif" %self.dt
            sci.movie("wtmp*.png",encoder='convert', fps=5, output_file=filename)
            
            

        
class Gauss_wave(wave_2d):
    def __init__(self,dt=0.2,Lx=10.0,Ly=10.0,T=80.0,Nx=30,Ny=30,b=0.0):
        wave_2d.__init__(self,Lx,Ly,T,Nx,Ny,dt,b, I_f = self.gauss)

def plug_I(x,y):
    I = zeros(shape(x))
    a = shape(x)
    dx = a[0]/20.0
    I[a[0]/2-dx:a[0]/2+dx,:] = I[a[0]/2-dx:a[0]/2+dx,:]+2
    return I
    

parser = argparse.ArgumentParser()
parser.add_argument("-s",action="store_true", help="use scalar solver") 
parser.add_argument("-b",action="store_true", help ="use the basic solver without damping and subsea terrain")
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
	


#plug_I(zeros((10,5),float), 8)

#p = wave_2d(10,10,80,30,30,10.0/(29*sqrt(2)),I_f = plug_I)
#p.solve_num()
w = wave_2d(Lx,Ly,T,Nx,Ny,dt,standing=True)
#w.make_exact()

w.solve_num()
#S = Standing_wave(0.1)
"""
def plug_I(x,y):
    return zeros(shape(x))
"""

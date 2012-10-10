from mayavi import mlab
from numpy import *
from scitools.std import *
import os,glob,time,re

infile = open('initial.txt','r')
varlist = [line for line in infile]
infile.close()
Nx=0;Ny=0;Lx=0;Ly=0
for i in varlist:
	if i.startswith('Nx'):
		Nx = int(i.split('=')[-1])
	elif i.startswith('Ny'):
		Ny = int(i.split('=')[-1])
	elif i.startswith('Lx'):
		Lx = int(i.split('=')[-1])
	elif i.startswith('Ly'):
		Ly = int(i.split('=')[-1])
	else:
		print "hummm huuuumm... something is not right.."
		
dx = Lx/float(Ny+1);
dy = Ly/float(Ny+1)
f = open('u0.txt','r')
u0 = [ map(float,line.split()) for line in f ]
u0 = array(u0)
f.close()
i=0
X,Y = meshgrid(linspace(0,Lx,(Nx/dx)),linspace(0,Ly,(Ny/dy)))
s = mlab.mesh(X[1:-1,1:-1], Y[1:-1,1:-1], u0[1:-1,1:-1])
for matrix in sorted(glob.glob('texttmp*.txt')):
	f = open(matrix,'r');l=[]
	l = [ map(float,line.split()) for line in f ]
	l=array(l)
	f.close()
	os.remove(matrix)
	print"hei paa du!"
	mlab.figure()
	mlab.mesh(X[1:-1,1:-1], Y[1:-1,1:-1], l[1:-1,1:-1], color=(0.0,0.75,1.0))
	mlab.savefig("wtmp%.4d.png" %i)
	i+=1

movie("wtmp*.png")
time.sleep(1)
os.remove('initial.txt');os.remove('u0.txt')
for i in glob.glob("wtmp*.png"):
	os.remove(i)

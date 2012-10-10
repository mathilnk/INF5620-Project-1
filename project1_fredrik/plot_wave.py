from mayavi import mlab
from numpy import *
from scitools.std import *
import os,glob,time,re

infile = open('initial.txt','r')
varlist = [line for line in infile]
infile.close()
Nx=0;Ny=0;Lx=0;Ly=0;dx=0;dy=0
for i in varlist:
	if i.startswith('Nx'):
		Nx = int(i.split('=')[-1])
	elif i.startswith('Ny'):
		Ny = int(i.split('=')[-1])
	elif i.startswith('Lx'):
		Lx = int(i.split('=')[-1])
	elif i.startswith('Ly'):
		Ly = int(i.split('=')[-1])
	elif i.startswith('dy'):
		dy = float(i.split('=')[-1])
	elif i.startswith('dx'):
		dx = float(i.split('=')[-1])
	else:
		print "hmm... something is not right.."

f = open('u0.txt','r')
u0 = [ map(float,line.split()) for line in f ]
f.close()
print Ly/dy,Ly

X,Y = meshgrid(linspace(dx,Lx-dx,(Nx/dx)-2),linspace(dy,Ly-dy,(Ny/dy)-2))
print X
s = mlab.mesh(X, Y, u0)
for matrix in glob('texttmp*.txt'):
	f = open(matrix,'r');l=[]
	l = [ map(float,line.split()) for line in f ]
	f.close()
	sys.remove(matrix)
	s.mlab_source.scalars = l
	mlab.savefig("wtmp%04d.png" %i)

movie("wtmp*.png")
time.sleep(1)
os.remove('initial.txt');os.remove('u0.txt')
for filename in glob('wtmp*.png'):
    os.remove(filename)
		

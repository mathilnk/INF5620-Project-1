from mayavi import mlab
from scitools.std import *
import os,glob,time

infile = open('somefile.txt','r')
varlist = [line for line in infile]
infile.close()
eval(varlist)
X,Y = meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny))

s = mlab.mesh(X,Y,u0)
for matrix in glob('texttmp*.txt'):
	f = open(matrix,'r');l=[]
	l = [ map(float,line.split(',')) for line in f ]
	f.close()
	sys.remove(matrix)
	s.mlab_source.scalars = l
	mlab.savefig("wtmp%04d.png" %i)

movie("wtmp*.png")
time.sleep(1)
for filename in glob('wtmp*.png'):
    os.remove(filename)
		

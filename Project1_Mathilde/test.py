
        
import glob,os,sys
from numpy import *

V_f_glob = lambda x,y: x*0
q_f_glob = lambda x,y: ones(shape(x))*0.7
f_f_glob = lambda x,y,t: x*0

class test:
    def __init__(self, a):
        self.a = a
        self.c = 2.0

class Test(test):
    def __init__(self,a):
        test.__init__(self,a)
        self.v = 5.0

b = Test(1.0)
print b.v

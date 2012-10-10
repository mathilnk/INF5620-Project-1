from scitools.std import *
x = linspace(0,10,100)
plot(x,cos(5*x*pi/10))
hold ('on')
plot(x, cos(5*x*pi/10+3))

t = raw_input(" ")

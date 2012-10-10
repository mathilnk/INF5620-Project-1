from scitools.std import *

movie("wtmp*.png")

print glob('wtmp*.png')
for filename in glob('wtmp*.png'):
    os.remove(filename)
print "done"
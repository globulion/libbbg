#!/usr/bin/python
"""
Compute tcf
"""
from sys import argv, exit
print __doc__
if len(argv)==1: exit()
from utilities import get_tcf
from scitools import filetable as ft
from numpy import *
from pylab import *

f = open(argv[1])
r = ft.read(f)[:,:2]
f.close()

np = 25000
o = open('o','w')
out = 'tcf.out'

log = ''
for l in range(np):
    log += '%6i %13.5E\n' % tuple(r[l])
log = log.replace('E','D')
o.write(log)
o.close()

print help(get_tcf)

get_tcf('o',nmax=25000,norgns=25000,ndels=10000,
            nskip=0,lprint=False,
            save=True,outfile=out)

# plot results
f = open(out)
r = ft.read(f)
f.close()

x    = r[:,0]
y    = r[:,1]
f = figure()
plt.plot(x, y   , '-' , label='simul')
plt.legend()
plt.savefig(out[:-3]+'.eps')


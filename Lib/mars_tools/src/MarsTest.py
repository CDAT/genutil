## Automatically adapted for numpy Apr 13, 2006 by convertcode.py

import numpy
from numpy import *
#from Numeric import *

from mars import mars as _mars
from mars import plot as _plot
from mars import catprt as _catprt
from mars import fmod as _fmod
#from mars import *

from VVutil.ColumnIO import *
from MarsResponseSurface import *
import shelve
import Sampling.LHS
import vvcollections.timer as T
#import timeit

inputFile = "MarsTest.dat"

n=98
print "number of observations: ", n
p=2
print "number of predictor variables (dimensions): ", p
nk=15
print "max number of basis functions: ", nk
mi=2
print "max_interaction: ", mi
m=2
print "model flag; 1: plot linear mars model, 2: plot cubic mars model: ", m
#x = numpy.array(range(npts*nparams), numpy.Float32)
ngc=100
print "number of raster points for computing curve estimates: ", ngc
ngs=40
print "number of raster points on each axis for computing curve estimates: ", ngs
icx=1
print "convex hull flag; 0: plot over entire range of arg limits, >0: plot only inside convex hull: ", icx

w = numpy.array(n*[1], dtype=float32)
#print "weights: ", weights
lx = numpy.array(p*[1], dtype=int32)
#print "lx: ", lx

dataDict = readData(inputFile)
vars = dataDict['vars']
data = dataDict['data']

x = numpy.array([data[vars.index('x1')],data[vars.index('x2')]])
x = numpy.transpose(x)
print "x.shape: ", x.shape
print "x: ", x

y = numpy.array([data[vars.index('y')]])
y = numpy.transpose(y)
print "y.shape: ", y.shape
#y = data[vars.index('y')]
print "y: ", y

fm, im = _mars(n, p, x, y, w, nk, mi, lx)

nc, crv, ns, srf = _plot(m, n, p, x, fm, im, ngc, ngs, icx, nk)

_catprt(m, nk, fm, im)

f = _fmod(m, n, p, x, fm, im)

print "f: ", f

mrs = MarsResponseSurface((fm,im))


try:
    s = shelve.open('oldmars.she','c')
except:
    print 'error'
    sys.exit()

#s['mars']=mrs
oldmrs = s['mars']
s.close()

print (mrs==oldmrs) #no problem

#check the evaluation at random points
dim = 2
npts = 100

box = []
for i in range(dim):
    box = box + [(0., 1.)]
d = Sampling.LHS.genUniformLHSdata(npts, box)

d = array(d).T
d = d.tolist()

t=T.Timer()
t.start()
newOut = map(mrs, d)
t.stop()
etnew = t.getElapsedTime()
print etnew

t.start()
oldOut = map(oldmrs, d)
t.stop()
etold = t.getElapsedTime()
print etold

print etold/etnew
print newOut==oldOut
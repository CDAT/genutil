import os, sys, logging, numpy, time, pdb, scipy, genutil
from genutil.pymars import logger, setlogging
from genutil.mars import logit as F_logit
from genutil.mars import xvalid as F_xvalid
from genutil.mars import mars as F_mars
from genutil.mars import fmod as F_fmod
from genutil.mars import setlog as F_setlog
from genutil.mars import setdf as F_setdf
from genutil.marstools.MarsResponseSurface import MarsResponseSurface

data = numpy.loadtxt("MarsTest.dat", skiprows=1)
x = data[:, 0:-1]
y = data[:, -1]

nk = 5
mi = 2
m = 1
lx = [1, 1]
logisticRegression, crossValidation = 0,0

print 'logisticRegression=', logisticRegression
print 'crossValidation=', crossValidation

# set the mars parameters
n, p = x.shape
print "x.shape: ", str(x.shape)
print "y.shape: ", str(y.shape)

w = numpy.array(n * [1], dtype=numpy.float64)
if lx == None:
    lx = p * [1]
lx = numpy.array(lx, dtype=numpy.int32)
print 'lx= ', str(lx)

F_setlog('rs_mars.log')
# logistic regression
F_logit(logisticRegression)
# cross validation
F_xvalid(crossValidation)
# reset df
F_setdf(3.0)

#run Friedman mars
fm, im = F_mars(n, p, x, y, w, nk, mi, lx)

#create response surface
#mrs = MarsResponseSurface(logfile, fm=fm, im=im)
mrs = MarsResponseSurface((fm,im))

#list of points for evaluation
pnt=numpy.array([(1.,1.)])
print mrs
print mrs(pnt)
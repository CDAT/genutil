import Sampling
from MARStools.ResponseSurfaceTools import *
from math import *

npts = 6
dim = 2
box = genBox(dim, 0., 6.0)
#s = Sampling.LHS.genUniformLHSdata(npts, box) 
s = Sampling.LHS.genGeometricLHSdata(npts, box, 3)     

for i in range(npts):
    print i, s[0][i]
for i in range(npts):
    print i, s[1][i]
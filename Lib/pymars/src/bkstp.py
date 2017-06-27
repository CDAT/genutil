import numpy
from .sweep import sweep
from .border import border
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def bkstp (d, nBF, xb, yb, al, dp):
    """backward stepwise elimination - algorithm 3"""
    #print 'bkstp'
    a = numpy.zeros(shape=nBF+1, dtype=FLOAT_DTYPE)
    rss = 9.9e30
    k = 0
    #determine the effect of eliminating each basis function?
    for i in range(1, nBF):
        if d[i,i] < 0.0: 
            s = 0.0
            for j in range(1, nBF):
                if d[j,j] < 0.0:
                    if j != i:
                        if j < i:
                            a0 = d[j,i]
                        else:
                            a0 = d[i,j]
                        s = s + dp[j]*(d[j,nBF] - a0*d[i,nBF]/d[i,i])**2
            s = d[nBF, nBF]-d[i,nBF]**2/d[i,i]-al*s
            if s <= rss: 
                rss = s
                k = i
    if k > 0: 
        d = sweep(d, nBF, k, 1.0)
    #same as the end of lsf1
    a0 = yb
    rss = 0.0
    for i in range(1, nBF):
        a[i] = 0.0
        if d[i,i] < 0.0: 
            a[i] = d[i,nBF]
            a0 = a0 - a[i]*xb[i]
            rss = rss + dp[i]*a[i]**2
    rss = d[nBF, nBF] - al*rss
    return rss, a, a0, k

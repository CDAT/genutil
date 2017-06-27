import numpy
from .sweep import sweep
from .border import border
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def lsf1(d, nBF, xb, yb, al):
    #print 'lsf1'

    dp = numpy.zeros(shape=nBF+2, dtype=FLOAT_DTYPE)
    a = numpy.zeros(shape=nBF+2, dtype=FLOAT_DTYPE)
    
    eps = 1.e-4
    for i in range(1, nBF):
        dp[i] = d[i,i]
        d[i,i] = d[i,i]*(1.0 + al)
    for i in range(1, nBF):
        if dp[i] > 0.0:
            im1 = i-1
            s = dp[i]
            j = 1
            while j <= im1:
                if d[j,j] < 0.0:
                    s = s+dp[j]*d[j,i]**2
                j = j+1
            #(j,) = numpy.where(d.diagonal()[1:i] < 0.0)
            #s = dp[i] + numpy.inner(dp[j],d[j,i]**2)
            if (d[i,i] - al*s)/dp[i] >= eps: 
                d = sweep(d, nBF, i, -1.0)
    #same as end of bkstp
    rss = 0.0
    a0 = yb
    for i in range(1, nBF):
        a[i] = 0.0
        if d[i,i] < 0.0:
            a[i] = d[i,nBF]
            a0 = a0 - a[i]*xb[i]
            rss = rss + dp[i]*a[i]**2
    rss = d[nBF, nBF] - al*rss
    return d, rss, a, a0, dp

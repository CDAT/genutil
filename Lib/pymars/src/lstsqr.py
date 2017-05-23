#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from .border import border
from .phi import phi
from .lsf import lsf
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def lstsqr(n, x, y, w, nk, yb, sw, tb, cm):
    #print 'lstsqr'
    sc = numpy.zeros(shape = (n+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)     
    d = numpy.zeros(shape = (nk+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)

    k=0
    for i in range(1, n+1):
        k=0
        for m in range(1, nk+1):
            if tb[1,m] != 0.0:
                k = k+1
                sc[i,k] = phi(m, x[i,:], tb, cm)

    for m in range(1, k+1):
        b = numpy.inner(w[1:], sc[1:,m])/sw
        for l in range(1, m):
            d[l,m] = numpy.inner(w[1:], (sc[1:,m] - b)*sc[1:,l])

        pp = sc[:,m] - b
        s = numpy.inner(w[1:], pp[1:]**2)
        a = numpy.inner(w[1:], pp[1:]*y[1:])
 
        d[m,m] = s
        d[m,k+1] = a
        d[m,k+2] = b
    a0 = a
    a = numpy.array(d[1:k+1,k+3])
    a = a.flatten()
    a = border(a)
    a, a0, gf = lsf(nk, k, k+1, yb, d, a, a0, s, 1)
    gof = gf/sw
    return gof

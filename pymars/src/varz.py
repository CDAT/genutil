#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars.jf import jf
def varz(j, nk, tb, ub, cst, nd):

    ub = tb.copy()
    nd = 0
    if j > 0:
        for m in range(1, nk+1):
            if ub[1,m] != 0.0:
                if jf(m, j, ub) != 0:
                    ub[1,m] = 0.0
                    nd = nd+1
    (m,) = numpy.where(ub[1,1:] != 0.0)
    m = m+1
    cst = 1.0 + ub[5,m].sum()

    return ub, cst, nd

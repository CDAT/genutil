#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars.nord import nord
from pymars.jfv import jfv
from pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def coll(nk, tb, lp, lv):
    #print 'coll'
    jv = numpy.zeros(shape = nk+1, dtype = INT_DTYPE)
    mo = 0
    for m in range(1, nk+1):
        if tb[1,m] != 0.0:
            mo = max(mo, nord(m, tb))
    if mo == 0: 
        lp[1,1] = 0
        return lp, lv, jv
    l1 = 1
    l2 = l1
    for mt in range(1, mo+1):
        l10 = l1
        for m in range(1, nk+1):
            if tb[1,m] != 0.0 and nord(m,tb) == mt:
                jv = jfv(m, tb, jv)
                jg = 0
                i = l10
                #taken from collf
                while i <= l1-1:
                    k = lp[2,i]-1
                    ig = (jv[1:mt+1] == lv[k+1:k+mt+1])
                    if ig.all():
                        jg = 1
                        lp[3,i] = lp[3,i] + 1
                        break
                    else:
                        i = i+1
                if jg == 0: 
                    lp[1,l1] = mt
                    lp[2,l1] = l2
                    lp[3,l1] = 1
                    k = l2-1
                    for i in range(1, mt+1):
                        lv[i+k] = jv[i]
                    l1 = l1 + 1
                    l2 = l2 + mt
    lp[1,l1] = 0
    return lp, lv, jv

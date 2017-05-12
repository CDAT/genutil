#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.nord import nord
from pymars.jf import jf
def efp(l, jv, nk, tb):
    #print 'efp'
    efp = 0.0
    for m in range(1, nk+1):
        if tb[1,m] != 0.0:
            if nord(m, tb) == l:
                k = 0
                for j in range(1, l+1):
                    if jf(m, jv[j], tb) != 1:
                        k = 1
                        break 
                if k != 1:
                    efp = efp + tb[5,m]
    return efp
#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from .icf import icf
from .nordc import nordc
from .jf import jf
def knts(l, nt, jv, jl, kv, nk, tb, cm, x, js):

    l1 = 0
    for m in range(1, nk+1):
        ICF, js = icf(m, tb, cm, jl, kv, js)
        if ICF != 0: 
            if nordc(1, m, tb, cm) == l:  
                j = 1
                FOUND = True
                while FOUND and j < l+1:
                    FOUND = jf(m, jv[j], tb)
                    j = j+1
                if FOUND:
                    ip = m
                    l1 = l1+1
                    while ip > 0:
                        t = tb[2,ip]
                        j = int(abs(t)+.1)
                        if cm[2*j] != 0.0: 
                            ip = int(tb[4,ip]+.1)
                        else:   
                            (k,) = numpy.where(jv[1:] == j)
                            k=k[0]+1

                            x[l1,k] = tb[3,ip]
                            sign = 1.0
                            if t < 0:
                                sign = -1.0
                            x[l1,l+k] = sign#(1.00,t)
                            ip = int(tb[4,ip]+.1)
    return x, js

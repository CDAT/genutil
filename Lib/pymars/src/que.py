import numpy
from .cue import cue
def que(jp, l, nt, jv, n, x, tc, t):

    l2 = l+l
    l3 = l2+l
    l4 = l3+l
    
    T = numpy.array(t)
    for i in range(1, n+1):
        if t[i] != 0.0: 
            q = 1.0
            for k in range(1, l+1):
                j = int(jv[k])
                q = q*cue(x[i,j], tc[jp,k+l], tc[jp,k], tc[jp,k+l2], tc[jp,k+l3], tc[jp,k+l4])
                if q == 0.0:
                    break 
            T[i] = q
    return T

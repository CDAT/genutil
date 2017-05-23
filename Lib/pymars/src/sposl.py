#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
def sposl(a, m, n, b):

    for k in range(1, n+1):
        t = 0.0
        if k > 1:
            t = numpy.inner(a[1:k,k], b[1:k])
        b[k] = (b[k] - t)/a[k,k]
    
    for kb in range(1, n+1):
        k = n+1-kb
        b[k] = b[k]/a[k,k]
        t = -b[k]
        if k > 1: 
            if t != 0.0: 
                b[1:k] = b[1:k] + t*a[1:k,k]
    return b

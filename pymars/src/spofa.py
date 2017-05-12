#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
def spofa(a, m, n, info):
    #print 'spofa'
    j1 = info
    for j in range(j1, n+1):
        info = j
        s = 0.0
        if j >= 2: 
            for k in range(1, j):
                u = 0.0
                if k > 1:
                    u = numpy.inner(a[1:k,k], a[1:k,j])
                t = a[k,j] - u
                t = t/a[k,k]
                a[k,j] = t
                s = s + t*t
        s = a[j,j] - s
        if s <= 0.0:
            return a, info
        a[j,j] = numpy.sqrt(s)
    info = 0
    return a, info
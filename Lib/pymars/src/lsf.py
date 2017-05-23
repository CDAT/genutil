#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy

def lsf(nk, m, M, yb, d, a, a0, gf, k1):
    from genutil.pymars import parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
    from .spofa import spofa
    from .sposl import sposl

    #print 'lsf'
    dp = numpy.zeros(shape = (m+1, m+1), dtype = FLOAT_DTYPE)
    
    EPS = parameters['lsf']['eps']
    gf = 9.9e30
    if d[m,m] <= 0.0:
        return a, a0, gf
    
    for k in range(k1, m+1):
        dp[1:k+1,k] = d[1:k+1,k]
        dp[k,k] = dp[k,k]*(1.0 + EPS)
    a[1:m+1] = d[1:m+1, M]

    info = k1
    dp, info = spofa(dp, nk, m, info)

    if info != 0:
        return a, a0, gf
    a = sposl(dp, nk, m, a)

    diag = d.diagonal()
    a0 = yb - numpy.inner(a[1:m+1], d[1:m+1,M+1])
    gf = -numpy.inner(a[1:m+1], d[1:m+1,M] + EPS*diag[1:m+1]*a[1:m+1])
    return a, a0, gf

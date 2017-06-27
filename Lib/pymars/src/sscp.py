import numpy
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE, debug
def sscp (n, nBF, points, y, w, yb, yv, sw):
    #print 'sscp'
    da = numpy.zeros(shape=nBF+1, dtype=FLOAT_DTYPE)
    d = numpy.zeros(shape=(nBF+1, nBF+1), dtype=FLOAT_DTYPE)
#    debug.info('    points='+repr(nBF))
#    for k in range(1, nBF):
#        debug.info(' k='+repr(k))
#        debug.info('    '+repr(points[1:n+1,k]))
    for k in range(1, nBF):
        #debug.info('    '+repr(points[1:n+1,k]))
        #debug.info(' ')
        s = numpy.inner(w[1:n+1], points[1:n+1,k])/sw
#        debug.info('    s='+repr(s))
        da[k] = s
        points[1:n+1,k] = points[1:n+1,k] - s
        for j in range(1, k+1):
            d[j,k] = numpy.inner(w[1:n+1], points[1:n+1,j]*points[1:n+1,k])
        d[k,nBF] = numpy.inner(w[1:n+1], (y[1:n+1]-yb)*points[1:n+1,k])
    d[nBF,nBF] = sw*yv
    return d, da

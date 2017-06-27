import numpy
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def sweep (a, nBF, k, fl):
    u = numpy.zeros(shape=nBF+1, dtype=FLOAT_DTYPE)
    c = a[k,k]
    u[1:k+1] = a[1:k+1, k]
    a[1:k+1, k] = 0.0

    u[k:nBF+2] = a[k, k:nBF+1]
    a[k, k:nBF+1] = 0.0

    u[k] = fl
    for i in range(1, nBF+1):
        for j in range(1, nBF+1):
            a[i,j] = a[i,j] - u[i]*u[j]/c
    return a

#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars import parameters
def sclato(n, p, x, xm, xs, cm):

    z = numpy.array(x)
    for j in range(1, p+1):
        j1 = int(cm[2*j]+.1)
        if j1 == 0:
            if xs[j] > 0.0:
                for i in range(1, n+1):
                    z[i,j] = xs[j]*x[i,j] + xm[j]
        else:        
            j1 = j1-1
            for i in range(1, n+1):
                l = int(x[i,j]+.1)
                z[i,j] = cm[l+j1]
    parameters['fmrs']['ifg'] = 0
    parameters['cmrs']['ifg'] = 0
    return z
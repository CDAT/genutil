#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars.icat import icat
from pymars.basisFunction import defaultBasisFunction, categoricalBasisFunction
from pymars import parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE

def fmrs(n, x, nk, az, tb, cm):
    """Compute the response surface at requested values."""
    y = numpy.zeros(n+1, dtype=FLOAT_DTYPE)
    for i in range(1, n+1):
        s = az
        for m in range(1, nk+1):
            if tb[1,m] != 0.0:
                phi = 1.0
                ip = m
                while ip > 0:
                    variable = tb[2,ip]
                    knot = tb[3,ip]
                    j = int(abs(variable)+.1) 
                    #very similar to phi
                    if cm[2*j] > 0.0:
                        if parameters['fmrs']['ifg'] == 0:
                            k = icat(x[i,j], j, cm)
                        else:   
                            k = int(x[i,j]+.1)
                        u = categoricalBasisFunction(variable, knot, cm, k, x[i,j])
                    else:
                        u = defaultBasisFunction(variable, knot, 1.0, x[i,j])
                    if u == 0.0:
                        phi = 0.0
                        break
                    else:
                        phi = phi*u
                        ip = int(tb[4,ip]+.1)
                s = s + tb[1,m]*phi
        y[i] = s
    return y
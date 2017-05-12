#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.basisFunction import defaultBasisFunction,categoricalBasisFunction
def phi(m, point, tb, cm):
    """ This is equation 38 from Friedman paper. """
    phi = 1.0
    ip = m
    while ip > 0:
        orient = tb[2,ip]
        knot = tb[3,ip]
        j = int(abs(orient)+.1)
        if cm[2*j] > 0.0:
            u = categoricalBasisFunction(orient, knot, cm, int(point[j]+.1), point[j])
        else:
            u = defaultBasisFunction(orient, knot, 1.0, point[j])
        if u <= 0.0:
            phi = 0.0
            return phi
        phi = phi*u
        ip = int(tb[4,ip]+.1)
    return phi
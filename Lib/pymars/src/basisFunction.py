import numpy
from .border import border
class FriedmanBasisFunction():
    """currently unused"""
    def __init__(self):
        self.hockeySticks = []
        self.DF = None
class HockeyStick():
    """currently unused"""
    def __init__(self, coef, variable, knot, orient):
        self.coefficient = coef
        self.variable = variable
        self.knot = knot
        self.orient = orient 
    def __call__(self, coef, x):
        """ hockey stick basis function: needs work"""
        sign = 1.0
        if self.orient < 0:
            sign = -1.0
        y = sign*(x - self.knot)
        y = y.clip(min=0.0)
        y = y*coef
        return y
def defaultBasisFunction(orient, knot, coef, x):
    """ hockey stick basis function"""
    sign = 1.0
    if orient < 0:
        sign = -1.0
    y = sign*(x - knot)
    y = y.clip(min=0.0)
    y = y*coef
    return y
def categoricalBasisFunction(orient, knot, cm, k, point):
    #k = int(point+.1)
    if k == 0:
        u = 0.0
    else:
        u = cm[k + int(knot+.1)]
    if orient < 0.0:
        if u == 0.0:
            u = 1.0
        else:
            u = 0.0
    return u
def flatBasisFunction(orient, knot, coef, cm, w, x):
    y = numpy.zeros(shape = x.shape, dtype = type(x))
    
    i = numpy.array(x[1:]+.1, dtype=numpy.int32)
    u = cm[i + int(knot+.1)]
    u = border(u)
    
    (c,) = numpy.where(coef[1:] > 0.0)
    if orient < 0.0:
        (j,) = numpy.where(u[1:][c] == 0.0)
    else:
        (j,) = numpy.where(u[1:][c] > 0.0)
    j = j+1 #increment since j starts at 0
    y[j] = coef[j]
    
    (j,) = numpy.where(w[1:][c] > 0.0)
    nw = len(j)
    (l,) = numpy.where(u[1:][j] == 0.0)
    n0 = len(l)
    
    return y, nw, n0
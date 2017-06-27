import numpy, time, genutil.pymars
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from .border import border
from .basisFunction import defaultBasisFunction, flatBasisFunction
def update(il, n, kr, x, y, w, sw, yb, orient, knot, cm, cmVAL, partialEval, d, dy):
    START = time.time()
    db = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    
    eps = 1.e-4  
    kp = kr+1
    
    coef = partialEval.clip(min=0.0)
    
    if il == 1:
        currentEval = coef*(x - knot)
    else:
        if cmVAL > 0.0:  
            currentEval, nw, n0 = flatBasisFunction(orient, knot, coef, cm, w, x)
            if n0 == 0 or n0 == nw:
                return kr, currentEval, d, dy
        else:
            currentEval = defaultBasisFunction(orient, knot, coef, x)
            
    #weighted mean of the new point       
    b = numpy.inner(w[1:], currentEval[1:])
    b = float(b/sw)

    d[1:,kp] = currentEval[1:] - b
    (i,) = numpy.where(currentEval[1:] > 0.0)
    q = w[1:]*currentEval[1:]
    q = border(q)
    s = numpy.inner(q[1:][i], currentEval[1:][i] - b)
    v = numpy.inner(q[1:][i], y[1:][i] - yb)
    i = i+1 #increment for relative address
    for j in range(1,kr+1):
        db[j] = numpy.inner(q[i], d[i,j])
    END = time.time()
    if 'update' in genutil.pymars.TIME.keys(): genutil.pymars.TIME['update'] += [(START,END)]
    if s <= 0.0:
        return kr, currentEval, d, dy
    
    dv = s
    s = s - numpy.inner(db[1:kr+1], db[1:kr+1])
    if s < eps*dv:
        return kr, currentEval, d, dy

    for j in range(1, kr+1):
        d[1:n+1,kp] = d[1:n+1,kp] - db[j]*d[1:n+1,j]

    s = 1.0/numpy.sqrt(s)
    dy[kp] = s*(v + numpy.inner(-db[1:kr+1], dy[1:kr+1]))
    d[1:n+1, kp] = s*d[1:n+1,kp]
 
    kr = kp
    return kr, currentEval, d, dy

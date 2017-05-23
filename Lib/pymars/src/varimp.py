#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import LOG, logger, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from .varz import varz
from .vp import vp
from .numprt import numprt
from .generalizedCrossValidation import generalizedCrossValidation as GCV
from .initialize import init
def varimp(n, p, x, y, w, nk, il, az, tb, cm):
    vip = numpy.zeros(shape=p+1, dtype = FLOAT_DTYPE)
    sw, wn, yb, yv = init(w, y)
    yv = yv/sw
    ip = nk + 4
    
    #dummy initialization
    cst = 0.0
    nd = 0

    #I'm creating a new array named ub which satisfies
    #the needs in a call to varz & vp. ub will not be
    #copied back to sc primarily because there is a 
    #mismatch in dimension & I don't see any subsequent
    #need for the array. This array seems to be only
    #useful in varimp.
    ub = numpy.zeros(shape = (5+1, nk+1), dtype = FLOAT_DTYPE)
    ub, cst, nd = varz(0, nk, tb, ub, cst, nd)
    if cst == 1.0: 
        g0 = 0.0
        if il > 0:
            g0 = yv
    else:  
        g0, ub = vp(n, x, y, w, nk, il, yb, sw, az, ub, cm)
    cst = GCV(1.0, cst, wn) 
    if il == 0: 
        g0 = (g0+yv)*cst
    else:  
        g0 = g0*cst
    
    for j in range(1, p+1):
        ub, cst, nd = varz(j, nk, tb, ub, cst, nd)
        if nd == 0: 
            vip[j] = g0
        else:    
            if cst == 1.0:
                  g = 0.0
                  if il > 0: 
                      g = yv
            else:   
                g, ub = vp(n, x, y, w, nk, il, yb, sw, az, ub, cm)

            cst = GCV(1.0, cst, wn)
            if il == 0: 
                g = (g + yv)*cst
            else:  
                g = g*cst
            vip[j] = g  
    if LOG:
        logger.info('\n gcv removing each variable:')
        numprt(p, vip)

    az = 0.0
    for j in range(1, p+1):
        vip[j] = numpy.sqrt(max(0.0, vip[j]-g0))
        az = max(az, vip[j])
    if az <= 0.0:
        return vip

    vip = 100.0*vip/az
    if LOG: 
        logger.info(' relative variable importance:')
        numprt(p, vip)
    return vip
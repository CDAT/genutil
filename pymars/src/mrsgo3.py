#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.decideNextVariable import decideNextVariable
from pymars.newb import isNewHS
from pymars import debug
def keepHS(tb,  ict, fln, bfIndex, bfIndexLast, 
           kcp, kcp0, nc, k1, kr, ssq, rsq, sw, defaultHS, proposedHS):
    """
    This function defines a new hockey stick function after some checks.
    """
    proposedHS[0] = (rsq + proposedHS[0])/sw #coefficient

    if ict == 0 and defaultHS[0] <= fln*proposedHS[0]:
        index = bfIndexLast
        HS = defaultHS
    else:
        index = bfIndex
        HS = proposedHS
    DECIDE = isNewHS(HS[1:], tb[2:5, 1:index])
    COMPARE = '<'
    CAT_CHECK= (True, ict, kcp0, kcp, nc)
    
    if ict == 0:
        bfIndex = bfIndexLast
        bfIndexLast = bfIndex-1
        kr = k1
        rsq = ssq

    return HS, DECIDE, COMPARE, CAT_CHECK, bfIndex, bfIndexLast,  kr, rsq

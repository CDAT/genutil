#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from nordc import *
from jfvc import *
from pymars.basisFunction import categoricalBasisFunction
from pymars.border import border
from pymars import debug
def icf(m, tb, cm, jl, kv, jv):

    icf = 0
    if tb[1,m] == 0.0 or nordc(2, m, tb, cm) != jl:
        return icf, jv
    if jl == 0: 
        icf = 1
        return icf, jv

    #The following requires more debugging
    #print 'icf'
    nv = 0
    JP = jv[jl+1:]
    JP = border(JP)
    #jv, jp = jfvc(2, m, tb, cm, nv, jv, jv[jl+1])#needs work
    nv, jv, JP = jfvc(2, m, tb, cm, nv, jv, JP)
    jv[jl+1:] = JP[1:]

    for j in range(1, jl+1):
        if abs(jv[j]) != abs(kv[1,j]):
            return icf, jv

    for j in range(1, jl+1):
        l1 = kv[2,j]
        l2 = jv[jl+j]
        k = 2*abs(jv[j])
        kk = jv[j]*kv[1,j]
        nc = int(cm[k+1]+.1) - int(cm[k]+.1) + 1
        for i in range(1, nc+1):
            z = cm[i+l2]
            if kk < 0:
                if z == 0.0: 
                    z = 1.0
                else: 
                    z = 0.0
            if cm[i+l1] != z: 
                return icf, jv
    icf = 1
    return icf, jv
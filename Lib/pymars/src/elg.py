from .nnord import nnord
from .nordc import nordc
from genutil.pymars import parameters
def elg(NEST, SETINT, jv, l, lx, tb, cm):
    #print 'elg'
    ic = parameters['elg']['ic']
    elg = False
    kx = abs(lx[jv])
    if kx == 0:
        return elg
    if l == 0: 
        elg = True
        return elg
    if kx == 2 or kx == 3: 
        if nnord(l,tb) > 0:
            return elg
    ip = l
    while ip > 0:
        jl = int(abs(tb[2,ip])+.1)
        ip = int(tb[4,ip]+.1)
    k = abs(lx[jl])
    jb = NEST.isnstr(jl)
    if (k == 2 or k == 3) and jb == 0:
        return elg
    if ic == 1: 
        if lx[jv] < 0 and nordc(1,l,tb,cm) > 0: 
            return elg
        if lx[jv] > 0 and nordc(2,l,tb,cm) > 0: 
            return elg
    else:
        if ic == 2:
            if lx[jv] > 0 and nordc(1,l,tb,cm) >= 2:
                 return elg
    ip = l
    while ip > 0:
        jl = int(abs(tb[2,ip])+.1)
        k = SETINT.intalw(jv, jl)
        if k == 0:
            return elg
        ip = int(tb[4,ip]+.1)
    elg = True
    return elg
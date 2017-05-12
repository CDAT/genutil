#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.nordc import nordc
from pymars.icf import icf
from pymars.jfvc import jfvc
def collf(nk, tb, cm, jl, kv, l1, l2, lp, lv, jv):

    mo = 0
    for m in range(1, nk+1):
        ICF, jv = icf(m, tb, cm, jl, kv, jv)
        if ICF != 0:
            mo = max(mo, nordc(1, m, tb, cm))
    if mo == 0:
       return l1, l2, lp, lv, jv
    for mt in range(1, mo+1):
        l10 = l1
        for m in range(1, nk+1):
            ICF, jv = icf(m, tb, cm, jl, kv, jv)
            if ICF != 0:
                if nordc(1, m, tb, cm) == mt:
                    nv = 0
                    nv, jv, jp = jfvc(1, m, tb, cm, nv, jv, jv)
                    jg = 0
                    i = l10
                    while i <= l1-1:
                        k = lp[2,i]-1
                        ig = (jv[1:mt+1] == lv[k+1:k+mt+1])
                        if ig.all():
                            jg = 1
                            lp[3,i] = lp[3,i] + 1
                            break
                        else:
                            i = i+1
                    if jg == 0:
                        lp[1,l1] = mt
                        lp[2,l1] = l2
                        lp[3,l1] = 1
                        k = l2-1
                        for i in range(1, mt+1):
                            lv[i+k] = jv[i]
                        l1 = l1+1
                        l2 = l2+mt
    return l1, l2, lp, lv, jv
#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars import LOG, logger, debug
from pymars.icf import icf
from pymars.nord import nord
from pymars.border import border
def purcat(nk, tb, cm, kp, kv, li, jv):
    #print 'purcat'
    lm = 1
    while kp[1, lm] >= 0:
        lm = lm + 1
    ll = 1
    li = 0

    while kp[1, ll] >= 0:
        jl = kp[1,ll]
        if jl <= 0:
            ll = ll + 1
        else: 
            ifg = 0
            jfg = ifg
            #debug.info('in purcat   kv='+repr((kv[1,1])))
            for m in range(1, nk+1):
                KV = numpy.array(kv[1:3, kp[2,ll]:kp[2,ll]+jl])
                KV = border(KV)
                #debug.info('before icf    kv='+repr((kv[1,1])))
                ICF, jv = icf(m, tb, cm, jl, KV, jv)
                #debug.info('after icf    ICF, kv='+repr((ICF, kv[1,1])))
                if ICF !=0:
                #if icf(m, tb, cm, jl, kv[1,kp[2,ll]], jv) != 0: 
                    if nord(m,tb) == jl: 
                        ifg = 1
                    else:
                        jfg = 1
            if ifg == 0: 
                if jfg == 0: 
                    logger.info(' bug in purcat - term not found.')
                    raise "stop"
                ll = ll + 1
            else:
                li = li + 1
                j = lm
                #print 'purcat 2'
                #debug.info('purcat 2, lm='+repr((lm)))
                while j >= li:
                    for i in [1,2,3,4,5]:
                        kp[i,j+1] = kp[i,j]
                    j = j-1
                lm = lm + 1
                ll = ll + 1
                for i in [1,2,3,4,5]:
                    kp[i,li] = kp[i,ll]
                kp[3,li] = 0
                kp[4,li] = 1
                kp[5,li] = 0
                if jfg == 1: 
                    ll = ll + 1
                else:
                    j = ll + 1
                    while j <= lm:
                        for i in [1,2,3,4,5]:
                            kp[i,j-1] = kp[i,j]
                        j = j + 1
                    lm = lm - 1

    return kp, li

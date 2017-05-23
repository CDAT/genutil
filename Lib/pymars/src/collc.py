#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from .nordc import *
from .jfvc import *
from .border import border
from genutil.pymars import debug
def collc (nk, tb, cm, kp, kv, jv):
    nv = 0
    kp[1,1] = 0
    kp[2,1] = 1
    l1 = 2
    l2 = 1
    mc = 0
    for m in range(1, nk+1):
        if tb[1,m] != 0.0:
            mc = max(mc, nordc(2, m, tb, cm))

    mt = 1
    while mt <=  mc:
        l10 = l1
        for m in range(1, nk+1):
            if tb[1,m] != 0.0 and nordc(2, m, tb, cm) == mt:
                #print 'collc'
                #needs work
                JP = numpy.array(jv[mt+1:])
                JP = border(JP)
                #nv, jv, jv[mt+1] = jfvc(2, m, tb, cm, nv, jv, JP)
                nv, jv, JP = jfvc(2, m, tb, cm, nv, jv, JP)
                jv[mt+1:] = JP[1:]
                jg = 0
                l1m1 = l1-1
                i = l10
                while i <=  l1m1:
                    k = kp[2,i]-1
                    #possible replacement
                    #ig = (abs(jv[1:mt+1]) == abs(lv[k+1:k+mt+1]))
                    ig = 0
                    for j in range(1, mt+1):
                        if abs(jv[j]) != abs(kv[1,k+j]):
                            ig = 1
                            break 
                    if ig == 0:
                        for j in range(1, mt+1):
                            m1 = kv[2,k+j]
                            m2 = jv[mt+j]
                            jj = abs(jv[j])
                            nc = int(cm[2*jj+1]+.1) - int(cm[2*jj]+.1) + 1
                            kk = jv[j]*kv[1,k+j]
                            for jk in range(1, nc+1):
                                #categorical basis function
                                z = cm[jk+m2]
                                if kk < 0:
                                    if z == 0.0:
                                        z = 1.0
                                    else:
                                        z = 0.0
                                if cm[jk+m1] != z:
                                    ig = 1
                                    break 
                            if ig == 1:
                                break 
                    if ig != 0:
                        i = i+1 
                    else:
                        jg = 1
                        break

                if jg == 0:
                    kp[1,l1] = mt
                    kp[2,l1] = l2
                    k = l2-1
                    for i in range(1, mt+1):
                        kv[1,i+k] = jv[i]
                        kv[2,i+k] = jv[i+mt]
                    l1 = l1+1
                    l2 = l2+mt
                #debug.info('in collc, m, kp='+repr((m, kp[1, m])))
        mt = mt+1
    kp[1,l1] = -1
    return kp, kv, jv
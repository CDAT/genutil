#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import LOG, logger, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE, debug
from .collc import *
from .purcat import *
from .collf import *
def ccoll(n, p, mi, nk, tb, cm, kp, kv, lp, lv):
    #print 'ccoll'
    jv = numpy.zeros(shape = ARRAY_SIZE, dtype = FLOAT_DTYPE)
#    debug.info('before collc, kp='+repr((kp[1, 1:nk+1])))
    kp, kv, jv = collc(nk, tb, cm, kp, kv, jv)

    #dummy assignment
    li = 0
#    debug.info('before purcat, kp='+repr((kp[1, 1:nk+1])))
    kp, li = purcat(nk, tb, cm, kp, kv, li, jv)
#    debug.info('after purcat, kp='+repr((kp[1, 1:nk+1])))
    ll = li+1
    l1 = 1
    l2 = l1
    while kp[1, ll] >= 0:
      	kp[4,ll] = l1
        jl = kp[1,ll]
        #taken from cubic
        KV = []
        if jl > 0:
            KV = numpy.array(kv[1:,kp[2,ll]:].flatten())
            KV = KV[0:2*jl]
            KV.shape = 2,jl
            KV = border(KV)
      	l1, l2, lp, lv, jv = collf(nk, tb, cm, jl, KV, l1, l2, lp, lv, jv)
      	kp[3,ll] = l1-kp[4,ll]
      	ll = ll+1
    lp[1,l1] = 0

    return kp, kv, lp, lv

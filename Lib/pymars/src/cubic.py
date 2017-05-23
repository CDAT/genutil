#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import LOG, logger, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE, debug
from .border import *
from .knts import *
from .side import *
from .que import *
from .lsf import lsf
from .initialize import init
from .generalizedCrossValidation import generalizedCrossValidation as GCV
def cubic (n, p, x, y, w, nk, tb, cm, kp, kv, lp, lv, bz, tc):

    t = numpy.zeros(shape = (n+1,nk+1), dtype = FLOAT_DTYPE)
    z = numpy.zeros(shape = (2+1,p+1), dtype = FLOAT_DTYPE)
    sc = numpy.zeros(shape = n+1, dtype = FLOAT_DTYPE)
    d = numpy.zeros(shape = (nk+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)
    js = numpy.zeros(shape = ARRAY_SIZE, dtype = FLOAT_DTYPE)
    
    big = 9.9e30
    sw, wn, yb, yv = init(w, y)
    yv = yv/sw
    
    (m,) = numpy.where(tb[1,1:] != 0.0)
    ni = len(m)
    if ni == 0:
      	bz = yb
      	u = GCV(yv, 1.0, wn)
     	if LOG:
            logger.info('\n piecewise cubic fit on %3i  basis functions, gcv = %12.4g' %(ni,u))
            return tb, cm, kp, kv, lp, lv, bz, tc
    nkp1 = nk+1
    nkp2 = nk+2
    nkp3 = nk+3
    lm = nkp3+nk
    for i in range(1, p+1):
      	xl = big
      	xr = -xl
      	for j in range(1, n+1):
      		xl = min(xl, x[j,i])
      		xr = max(xr, x[j,i])
      	z[1,i] = xl
      	z[2,i] = xr
    #similar to logitc
    ll = 1
    la = ll
    l1 = la
    lt = 0
    while kp[1, ll] >= 0:
        sc[1:] = 1.0

        FLAG1 = (kp[1,ll] > 0)
        FLAG2 = (kp[3,ll] <= 0)
#        debug.info('in cubic  ll, FLAG1, FLAG2='+repr((ll, FLAG1, FLAG2)))
        if FLAG1:
            jl = kp[1,ll]
            for il in range(1, jl+1):
      			k = kp[2,ll]+il-1
      			jj = kv[1,k]
      			j = abs(jj)
      			kk = kv[2,k]
      			for i in range(1, n+1):
	      			if sc[i] != 0.0:
	      				ic = int(x[i,j]+.1)
	      				sc[i] = cm[ic+kk]
                        if jj < 0:
    	      				if sc[i] == 0.0:
    	      					sc[i] = 1.0
    	      				else: 
    	     				    sc[i] = 0.0
            if FLAG2:
                lt = lt+1
                kp[5,ll] = 0
                t[1:n+1, lt] = sc[1:n+1]

        if not FLAG2: 
            kp3 = kp[3,ll]
            kp[5,ll] = la
            for m in range(1, kp3+1):
                l = lp[1,l1]
                nt = lp[3,l1]

                jv = numpy.array(lv[lp[2,l1]:lp[2,l1]+l])
                jv = border(jv)
                jl = kp[1,ll]
                KV = []
                if jl > 0:
                    KV = kv[1:,kp[2,ll]:].flatten()
                    KV = KV[0:2*jl]
                    KV.shape = 2,jl
                    KV = border(KV)
                TC_dim = len(tc[la:])/nt
                TC = tc[la:la+nt*TC_dim].copy()
                TC.shape = nt,TC_dim
                TC = border(TC)
                TC, js = knts(l, nt, jv, jl, KV, nk, tb, cm, TC, js)
                TC = side(l, nt, jv, z, TC)

                for jp in range(1, nt+1):
                    lt = lt+1
                    t[1:n+1,lt] = sc[1:n+1]
                    T = numpy.array(t[1:n+1,lt])
                    T = border(T)
                    T = que(jp, l, nt, jv, n, x, TC, T) 
                    t[1:n+1,lt] = T[1:]
                TC = TC[1:,1:].T.flatten()
                tc[la:la+nt*TC_dim] = TC
                l1 = l1 + 1
                la = la + nt*(5*l+1)
        ll = ll + 1
    
#    debug.info('in cubic   lt='+repr((lt)))
    for j in range(1, lt+1):
        u = 0.0
        s = numpy.inner(w[1:], t[1:n+1,j])/sw
        d[j,nkp2] = s
        t[1:n+1,j] = t[1:n+1,j] - s

        s = numpy.inner(w[1:], (y[1:] - yb)*t[1:n+1,j])
        d[j,nkp1] = s
        for k in range(1, j+1):
            s = 0.0
            s = numpy.inner(w[1:], t[1:n+1,k]*t[1:n+1,j])
            d[k,j] = s
#        debug.info('j, d='+repr((j, d[1:j+1,j])))

    A = numpy.array(d[1:lt+1,lm])
    A = border(A)
    A, s, u = lsf(nk, lt, nk+1, yb, d, A, s, u, 1)
#    debug.info('from lsf u, lt, d='+repr((u, lt, d[lt,lt])))
    (i,) = numpy.where(tb[1,1:] != 0.0)
    eft = 1.0 + tb[5,i+1].sum()
    u = GCV(u/sw + yv, eft, wn)
    
    bz = s
    ll = 1
    l1 = ll
    le = la-1
    la = 0
    lt = la
    while kp[1, ll] >= 0:
        if kp[1,ll] != 0  or kp[3,ll] > 0:
             if kp[3,ll] <= 0:
                  le = le + 1
                  kp[3,ll] = -le
                  lt = lt + 1
                  tc[le] = A[lt]
             else:
                 kp3 = kp[3,ll]
                 for m in range(1, kp3+1):
                     nt = lp[3,l1]
                     la = la + 5*lp[1,l1]*nt
                     for i in range(1, nt+1):
                         lt = lt + 1
                         tc[i+la] = A[lt]
                     la = la + nt
                     l1 = l1 + 1
        ll = ll + 1
    if LOG: 
        logger.info('\n piecewise cubic fit on %3i  basis functions, gcv = %12.4g' %(lt,u))
    return tb, cm, kp, kv, lp, lv, bz, tc

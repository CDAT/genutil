import numpy
from genutil.pymars import LOG, logger, debug, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from .border import border, copyback
from .phi import phi
from .coll import coll
from .nord import nord
from .jf import jf
from .exch import exch
from .efp import efp
from .varf import varf
from .lsf import lsf
from .initialize import init
from .generalizedCrossValidation import generalizedCrossValidation as GCV
def anova (n, x, y, w, nk, tb, cm):

    if not LOG:
        return 
    t = numpy.zeros(shape = (n+1,nk+1), dtype = FLOAT_DTYPE)        
    d = numpy.zeros(shape = (nk+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)
    lv = numpy.zeros(shape = ARRAY_SIZE, dtype = INT_DTYPE) 
    lp = numpy.zeros(shape = (3+1, nk+3), dtype = INT_DTYPE)
    a = numpy.zeros(shape = nk+1, dtype = FLOAT_DTYPE)
    TB = numpy.zeros(shape = nk+1, dtype = FLOAT_DTYPE)
    
    sw, wn, yb, u = init(w, y)
    yv = u/sw

    (m,) = numpy.where(tb[1,1:] != 0.0)
    m = m+1
    eft = 1.0 + tb[5,m].sum()

    ni = 0
    for m in range(1, nk+1):
        if tb[1,m] != 0.0:
            ni = ni+1
            for j in range(1, n+1):
                t[j,ni] = phi(m, x[j,:], tb, cm)
            s = numpy.inner(w[1:n+1], t[1:n+1,ni])/sw
            t[1:n+1, ni] = t[1:n+1, ni] - s

            for i in range(1, ni+1):
                d[i,ni] = numpy.inner(w[1:n+1], t[1:n+1,i]*t[1:n+1,ni])
            d[ni,nk+1] = numpy.inner(w[1:n+1], t[1:n+1,ni]*(y[1:n+1]-yb))
            #d[ni,nk+2] = tb[1,m]
            TB[ni] = tb[1,m]
            lv[ni] = m

    if ni == 0:
      	logger.info(' estimated optimal model = response mean.')
      	return
    for m in range(1, ni+1):
      	t[m,1] = lv[m]
    logger.info(' anova decomposition on %3i basis functions:\n'  \
    '  fun. std. dev.     -gcv    #bsfns  #efprms  variable(s)' %(ni))

    lp, lv, JV = coll(nk, tb, lp, lv)
    (m,) = numpy.where(lp[1,1:] != 0)
    na = len(m)
    
    m = 1
    if na == 1:
        k2 = lp[2,m]
      	i2 = lp[1,m] + k2 - 1
      	efm = eft-1.0
      	u = GCV(yv, 1.0, wn)
      	s = numpy.sqrt(varf(nk, d, TB, sw, 1, ni))
        lv_str = ' %3i  %10.4g  %10.4g  %2i    %4.2g    '%(m, s, u, lp[3,m], efm)
        for K in lv[k2:i2+1]:
            lv_str = lv_str + '%4i'%(K)
      	logger.info(lv_str)
      	return 

    for m in range(1, na+1):
        k2 = lp[2,m]
      	l = lp[1,m]
      	i2 = l + k2 - 1
      	ll = k2 - 1
      	np = ni
      	for im in range(1, ni+1):
            i = int(t[im,1] + .1)
            if nord(i,tb) != l:
      			t[im,2] = 0.0
            else:
                k = 0
                for j in range(1, l+1):
                    if jf(i, lv[ll+j], tb) != 1:
                        k = 1
                        break
                if k == 1:
	      			t[im,2] = 0.0
                else:
                    t[im,2] = 1.0
                    np = np-1
        EXECUTE = True
        while EXECUTE:
            k = 0
            for i in range(1, ni):
                if t[i,2] > t[i+1,2]:
                    k = 1
                    #the following code was changed by some
                    #renigade from the original fortran
                    if ni < nk+2:
                        T = numpy.array(t[1:nk+2+1, 1])
                        T = border(T)
                        b = numpy.array(t[1:nk+2+1, 2])
                        b = border(b)
                        d, T, b = exch(nk, nk+2, i, d, T, b) #t[1,2])
                        t[1:nk+2+1, 1] = T[1:]
                        t[1:nk+2+1, 2] = b[1:]
                    else:
                        T = numpy.array(t[1:ni+1, 1])
                        T = border(T)
                        b = numpy.array(t[1:ni+1, 2])
                        b = border(b)
                        d, T, b = exch(nk, ni, i, d, T, b)#t[1,2])
                        t[1:ni+1, 1] = T[1:]
                        t[1:ni+1, 2] = b[1:]
      		EXECUTE = (k == 1)

    	A, s, u = lsf(nk, np, nk+1, 0.0, d, a[0:np+1], s, u, 1)

        L = lp[1,m]
        JV = numpy.array(lv[lp[2,m]:lp[2,m]+L])
        JV = border(JV)
      	efm = efp(L, JV, nk, tb)
        
      	u = GCV(u/sw + yv, eft-efm, wn)
      	s = numpy.sqrt(varf(nk, d, TB, sw, np+1, ni))
        lv_str = ' %3i %10.4g %10.4g     %2i    %4.2g    ' %(m, s, u, lp[3,m], efm)
        for K in lv[k2:i2+1]:
            lv_str = lv_str + '%4i'%(K)
      	logger.info(lv_str)
    return 

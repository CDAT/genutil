#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import debug, LOG, logger, FLOAT_DTYPE
from .border import border
from .coll import coll
from .cptb import CPTB
from .nord import nord
from .jf import jf
from .vp import vp
from .efp import efp
from .initialize import init
from .generalizedCrossValidation import generalizedCrossValidation as GCV
def anoval(n, x, y, w, nk, il, az, tb, cm, lp, lv):
    
    if not LOG:
	    return 

    TB = numpy.zeros(shape = (5+1,nk+1), dtype = FLOAT_DTYPE)

    sw, wn, yb, yv = init(w, y)
    yv = yv/sw

    (m,) = numpy.where(tb[1,1:] != 0.0)
    ni = len(m)
    eft = 1.0 + tb[5,m+1].sum()

    if ni == 0:
     	logger.info('\n estimated optimal model = response mean.')
      	return 

    logger.info('\n logit anova decomposition on %3i  basis functions:'\
				'\n  fun.    -gcv    #bsfns  #efprms  variable(s)' %(ni))

    lp, lv, JV = coll(nk, tb, lp, lv)
    m = 1
    while lp[1,m] != 0:
      	m = m+1
    na = m-1
    m = 1
    if na == 1:
        k2 = lp[2,m]
        i2 = lp[1,m] + k2 - 1 
        efm = eft-1.0
        u = GCV(yv, 1.0, wn)
        LV = ''
        for i in range(k2, i2+1):
            LV = LV + '%4i'%(lv[i])
	    logger.info('%3i   %12.4g     %2i     %4.1d       '
	    			%(m,u,lp[3,m],efm)+LV)
	    return 
    ip = nk+4
    CPTB0 = CPTB(nk)
    for m in range(1, na+1):
      	k2 = lp[2,m]
      	l = lp[1,m]
      	i2 = l+k2-1
      	ll = k2-1
        CPTB0.ub = TB
      	CPTB0.cptb(nk, tb)
      	for i in range(1, nk+1):
      		if tb[1,i] != 0.0:
	      		if nord(i,tb) == l:
		      		k = 0
		      		for j in range(1, l+1):
		      			if jf(i, lv[ll+j], tb) != 1:
		      				k = 1
		      				break
		     		if k != 1:
		      			CPTB0.setz(i)    
      	u, TB = vp(n, x, y, w, nk, il, yb, sw, az, TB, cm)

        L = lp[1,m]
        LV = lv[lp[2,m]:lp[2,m]+L+1]
        LV = border(LV)
      	efm = efp(L, LV, nk, tb)
      	u = GCV(u, (eft-efm), wn)
        LVstr = ''
        for i in range(k2, i2+1):
            LVstr = LVstr + '%4i'%(lv[i]) #should this be LV not lv
     	logger.info('%3i   %12.4g     %2i     %4.1f       ' 
					%( m,u,lp[3,m],efm)+ LVstr)
    return

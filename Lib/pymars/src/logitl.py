import numpy
from .phi import phi
from .que import que
from .lsf import lsf
from .border import border

class LOGITL:
    def __init__(self):
        #in the original code niter was set to 30; if you set it to 30
        #it appears that the calculation becomes unstable, that is 
        #the answers between fortran & python are different.
        
        self.niter = 25
        self.wm  = .0001
        self.thr = .0001
        return
    def logitl(self, n, x, y, w, nk, il, az, tb, cm):
        from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
        
        sc = numpy.zeros(shape = (n+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)     
        d = numpy.zeros(shape = (nk+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)
        TB = numpy.zeros(shape = nk+1, dtype = numpy.float64)
        
        k = 0
        for i in range(1, n+1):
            k = 0
            for m in range(1, nk+1):
    	      	if tb[1,m] != 0.0:
                    k = k+1
                    sc[i,k] = phi(m, x[i,:], tb, cm)
        if k == 0:
	      	az = numpy.log(az/(1.0-az))
	      	return az, tb, sc, d
        mk = k
        a = az
        jnt = 1
 
        k = 0
        for m in range(1, nk+1):
            if tb[1,m] != 0.0:
                k = k+1
                TB[k] = tb[1,m]
                
        az, tb, TB = self.calc(il, n, y, w, nk, mk, a, tb, sc, d, TB)

        k = 0
        for m in range(1, nk+1):
            if tb[1,m] != 0.0:
                k = k+1
                tb[1,m] = TB[k]

        return az, tb, sc
	      
    def logitc(self, n, x, y, w, nk, il, cm, tb, kp, kv, lp, lv, bz, tc):
        from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
        
        sc = numpy.zeros(shape = (n+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)
        d = numpy.zeros(shape = (nk+1, ARRAY_SIZE), dtype = FLOAT_DTYPE)
        ss = numpy.zeros(shape = n+1, dtype = FLOAT_DTYPE)
        TB = numpy.zeros(shape = nk+1, dtype = FLOAT_DTYPE)
        #print 'logitc'
        ll = 1
        la = ll
        l1 = la
        lt = 0
        while kp[1, ll] >= 0:
            for i in range(1, n+1):
                ss[i] = 1.0
	        FLAG1 = (kp[1,ll] > 0)
	      	FLAG2 = (kp[3,ll] <= 0)

	      	if FLAG1:
		      	jl = kp[1,ll]
		      	for il in range(1, jl+1):
		      		k = kp[2,ll]+il-1
		      		jj = kv[1,k]
		      		j = abs(jj)
		      		kk = kv[2,k]
		      		for i in range(1, n+1):
			      		if ss[i] != 0.0:
			      			ic = int(x[i,j] + .1)
			      			ss[i] = cm[ic+kk]
			      			if jj < 0:
			      				if ss[i] == 0.0:
			      					ss[i] = 1.0
			      				else:
			     					ss[i]=0.0
            if FLAG1 and FLAG2:
                lt = lt+1
                for i in range(1, n+1):
				    sc[i,lt] = ss[i]
            if not FLAG2: 
                kp3 = kp[3,ll]
                for m in range(1, kp3+1):
                    l = lp[1,l1]
                    nt = lp[3,l1]
                    for jp in range(1, nt+1):
                        lt = lt+1
                        for i in range(1, n+1):
				            sc[i,lt] = ss[i]

                        JV = lv[lp[2,l1]:lp[2,l1]+l]
                        JV = border(JV)
                        
                        TC_dim = len(tc[la:])/nt
                        TC = tc[la:la+nt*TC_dim]
                        TC.shape = TC_dim, nt
                        TC = TC.T
                        TC = border(TC)
                        
                        T = sc[1:n+1, lt]
                        T = border(T)
                        T = que(jp, l, nt, JV, n, x, TC, T)
                        sc[1:n+1, lt] = T[1:]
                        
                    l1 = l1 + 1
                    la = la + nt*(5*l+1)
            ll = ll+1
	    if lt == 0:
	      	bz = numpy.log(bz/(1.0-bz))
	     	return bz, tb, tc
	    mk = lt
	    a = bz
	    jnt = 2
        TB = self.cin(mk, kp, lp, tc, TB)
        bz, tb, TB = self.calc(il, n, y, w, nk, mk, a, tb, sc, d, TB)
        tc = self.cout(mk, kp, lp, tc, TB)

        return bz, tb, tc
    def cin(self, mk, kp, lp, tc, TB):

        ll = 1
        l1 = ll
        la = 0
        lt = la
        while kp[1,ll] >= 0:
            if kp[1,ll] != 0 or kp[3,ll] > 0:
                if kp[3,ll] <= 0:
                    lt = lt+1
                    TB[lt] = tc[-kp[3,ll]]
                else:
                    kp3 = kp[3,ll]
                    for m in range(1, kp3+1):
                        nt = lp[3,l1]
                        la = la + 5*lp[1,l1]*nt
                        for i in range(1, nt+1):
                            lt = lt+1
                            TB[lt] = tc[i+la]
                        la = la+nt
                        l1 = l1+1
            ll = ll+1
        return TB
    def cout(self, mk, kp, lp, tc, TB):
     
        ll = 1
        l1 = ll
        la = 0
        lt = la
        while kp[1,ll] >= 0:
            if kp[1,ll] != 0 or kp[3,ll] > 0:
                if kp[3,ll] <= 0:
                    lt = lt + 1
                    tc[-kp[3,ll]] = TB[lt]
                else:
                    kp3 = kp[3,ll]
                    for m in range(1, kp3+1):
                        nt = lp[3,l1]
                        la = la + 5*lp[1,l1]*nt
                        for i in range(1, nt+1):
                            lt = lt + 1
                            tc[i+la] = TB[lt]
                        la = la + nt
                        l1 = l1 + 1
            ll = ll + 1
        return tc
    def calc(self, il, n, y, w, nk, mk, a, tb, sc, d, TB):        
        mkp1 = mk+1
        mkp2 = mk+2
        mkp3 = mk+3
        mkp4 = mk+4
        iter = 0

        EXECUTE = True
        while EXECUTE:
            iter = iter+1
            b = 0.0
            sw = b
            yb = sw
            for i in range(1, n+1):
                sc[i,mkp3] = a + numpy.inner(TB[1:mk+1],sc[i,1:mk+1])
                pp = 1.0/(1.0 + numpy.exp(-sc[i,mkp3]))
                ww = max(pp*(1.0-pp), self.wm)
                sc[i,mkp3] = sc[i,mkp3] + (y[i]-pp)/ww
                if il == 2:
                    ww = ww**2
                ww = ww*w[i]
                sc[i,mkp2] = ww
                sw = sw + ww
                yb = yb + ww*sc[i,mkp3]
                if iter > 1:
                    b = b + abs(pp-sc[i,mkp1])
                sc[i,mkp1] = pp
            EXECUTE = not(iter > self.niter or (iter > 1 and b/n < self.thr))
            if EXECUTE:
                yb = yb/sw      
                for m in range(1, mk+1):
                    b = numpy.inner(sc[1:n+1, mkp2], sc[1:n+1,m])/sw
                    for l in range(1, m):
				      	d[l,m] = numpy.inner(sc[1:n+1, mkp2], (sc[1:n+1,m]-b)*sc[1:n+1,l])
                    pp = sc[1:n+1, m] - b
                    s = numpy.inner(sc[1:n+1, mkp2], pp**2)
                    a = numpy.inner(sc[1:n+1, mkp2], pp*sc[1:n+1, mkp3])
                    d[m,m] = s
                    d[m,mkp1] = a
                    d[m,mkp2] = b
                TB, a, s = lsf(nk, mk, mk+1, yb, d, TB, a, s, 1)
        return a, tb, TB

#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def atoscl(n, p, w, x, lx, mm,  cm):
    xm = numpy.zeros(shape = p+1, dtype = FLOAT_DTYPE)
    xs = numpy.zeros(shape = p+1, dtype = FLOAT_DTYPE)
    z = numpy.zeros(shape = x.shape, dtype = FLOAT_DTYPE)
    
    sw = 0.0
    sw = w[1:].sum()
    ip = 0
    nct = ip
    for i in range(1, p+1):
        if lx[i] == 0:
            xm[i] = 0.0
            xs[i] = xm[i]
        elif lx[i] < 0: 
      		nc = 0
      		xm[i] = ip
      		j = 1
      		nct = nct+1
      		while j <= n:
	    		j0 = j
	      		if j < n:
	     			while x[mm[j+1,i],i] <= x[mm[j,i],i]:
	      				j = j+1
	      				if j >= n:
	      					break 
	     		ip = ip+1
	      		cm[ip] = x[mm[j,i],i]
	      		nc = nc+1
	      		for k in range(j0, j+1):
	      			z[mm[k,i],i] = nc
	      		j = j+1
     		xs[i] = nc
        else:
            t = 0.0
            s = numpy.inner(w[1:], x[1:,i])/sw
            xm[i] = s
            z[1:,i] = x[1:,i] - s
            t = numpy.inner(w[1:], z[1:,i]**2)
            xs[i] = 1.0
            
            if t > 0.0:
                t = numpy.sqrt(t/sw)
                xs[i] = t
                t = 1.0/t
                z[1:,i] = t*z[1:,i]

    n2 = 2*p+1
    if nct == 0:
        cm[1:n2+1] = 0.0
      	return z, xm, xs, cm
   
    n2p1 = n2+1
    i = ip
    while 1 <= i:
        cm[i+n2] = cm[i]
        i = i-1 
    j = 0
    i = 2
    while i <= n2:
        j = j+1
        if lx[j] < 0: 
            cm[i] = xm[j]+n2p1
            cm[i+1] = cm[i]+xs[j]-1.0
        else:
            cm[i] = 0.0
            cm[i+1] = cm[i]
        i = i+2
    cm[1] = nct
    parameters['fmrs']['ifg'] = 1
    parameters['cmrs']['ifg'] = 1
    return z, xm, xs, cm

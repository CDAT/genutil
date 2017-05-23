#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import LOG, logger, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def rspnpr(il, n, y, w):
    wm = numpy.zeros(shape = 3, dtype = FLOAT_DTYPE)
    if not LOG:
		return
    if il == 1:
    	wm[1] = 0.0
      	wm[2] = wm[1]

      	for i in range(1, n+1):
      		k = int(y[i] + 1.1)
      		wm[k] = wm[k] + w[i]
      	wt = wm[1] + wm[2]
      	wm[1] = wm[1]/wt
      	wm[2] = wm[2]/wt
      	logger.info('(/ binary (0/1) response:  mass(0) = %12.4g,   mass(1) = %12.4g)'%(wm[1],wm[2]))
      	return
      	
    logger.info(' ordinal response:\n')
    logger.info('          min          n/4          n/2          3n/4         max')
    m = y.argsort()

    logger.info( '%12.4g %12.4g %12.4g %12.4g %12.4g ' %(y[m[1]], y[m[n/4]], y[m[n/2]], y[m[n-n/4]], y[m[n]]))
    return

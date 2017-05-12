#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars import LOG, logger
def ordpr(n,p,x,lx,m):
	if not LOG:
		return

	(j,) = numpy.where(lx[1:] > 0)
	no = len(j)

	if no == 0:
		return

	logger.info('\n there are %3i ordinal predictor variables.' %(no))
	logger.info('  var          min         n/4         n/2        3n/4            max ')
	n1=n/4
	n2=n/2
	n3=n-n1
	for j in range(1, p+1):
		if lx[j] > 0:
			logger.info(' %3i %12.4g %12.4g %12.4g %12.4g %12.4g' 
					    %(j, x[m[1,j],j], x[m[n1,j],j], x[m[n2,j],j], x[m[n3,j],j], x[m[n,j],j]))
	return

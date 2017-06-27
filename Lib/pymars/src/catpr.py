import numpy
from genutil.pymars import LOG, logger, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def catpr(n, p, x, cm):

    if not LOG:
        return
    nct = int(cm[1]+.1)
    if nct == 0: 
        return
    
    mm = numpy.zeros(shape = ARRAY_SIZE, dtype = INT_DTYPE)
    n2 = 2*p+1
    np = 0
    logger.info('\n there are %3i categorical predictor variables.'%(nct))

    i = 2
    while i <= n2:
        np = np+1
        j1 = int(cm[i]+.1)
        if j1 != 0:
            j2 = int(cm[i+1]+.1)
            nv = j2-j1+1
            for j in range(1, nv+1):
                mm[j] = 0
            for j in range(1, n+1):
                ic = int(x[j,np]+.1)
                mm[ic] = mm[ic]+1
            logger.info(' categorical variable %3i has %3i values.'%(np,nv))
            logger.info(' value     internal code     counts')
            k = 0
            for j in range(j1, j2+1):
                k = k+1
                logger.info('%6.0f %13i %15i' %(cm[j], k, mm[k]))
        i=i+2
    return
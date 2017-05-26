import numpy
from genutil.pymars import LOG, logger
def computeErrors(id, x, y, n):
    print id + ' L1 Absolute error: %10g' % L1AbsError(x, y, n)
    
    print id + ' L1 Relative error: %10g' % L1RelError(x, y, n)
            
    print id + ' L2 Absolute error: %10g' % L2AbsError(x, y, n)
            
    print id + ' L2 Relative error: %10g' % L2RelError(x, y, n)
            
    print id + ' Linf Absolute error: %10g' % LinfAbsError(x, y, n)
            
    logger.info('L1 Absolute error: %10g' % L1AbsError(x, y, n))
    logger.info('L1 Relative error: %10g' % L1RelError(x, y, n))
    logger.info('L2 Absolute error: %10g' % L2AbsError(x, y, n))
    logger.info('L2 Relative error: %10g' % L2RelError(x, y, n))
    logger.info('Linf Absolute error: %10g' % LinfAbsError(x, y, n))   
def L1RelError(x, y, n):
    n = min(len(x),len(y))
    sum = 0.0
    for i in range(0,n-1):
        z = max(x[i],y[i])
        if z != 0.0:
            u = (y[i] - x[i]) / z
            sum = abs(u) + sum
            u = 0.0
    return numpy.sqrt(sum/n)

#L1 Absolute Error
def L1AbsError(x, y, n):
    n = min(len(x),len(y))
    sum = 0.0
    for i in range(0,n-1):
        z = max(x[i],y[i])
        if z != 0.0:
            u = y[i] - x[i]
            sum = abs(u) + sum
            u = 0.0
    return numpy.sqrt(sum/n)

#L2 Absolute Error
def L2AbsError(x, y, n):
    n = min(len(x),len(y))
    sum = 0.0
    for i in range(0,n-1):
        z = max(x[i],y[i])
        if z != 0.0:
            u = y[i] - x[i]
            sum = u**2 + sum
            u = 0.0
    return numpy.sqrt(sum/n)
        
#L2 Relative Error
def L2RelError(x, y, n):
    n = min(len(x),len(y))
    sum = 0.0
    for i in range(0,n-1):
        z = max(x[i],y[i])
        if z != 0.0:
            u = (y[i] - x[i]) / z
            sum = u**2 + sum
            u = 0.0
    return numpy.sqrt(sum/n)

#Linf Absolute Error
def LinfAbsError(x, y, n):
    n = min(len(x),len(y))
    mx = 0.0
    for i in range(0,n-1):
        z = max(x[i],y[i])
        if z != 0.0:
            u = y[i] - x[i]
            if abs(u) > mx:
                mx = abs(u)
            u = 0.0
    return mx
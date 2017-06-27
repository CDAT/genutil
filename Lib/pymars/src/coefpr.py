import numpy
from .org import org
from .border import border
from genutil.pymars import LOG, logger
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def coefpr(nk, az, tb, cm, xs):
    #print 'coefpr'
    a = numpy.zeros(shape=6 + 1, dtype=FLOAT_DTYPE)
    i2 = 0
    limit = 7
    cnt = -1
    while i2 < nk:
        if i2 == 0: 
            i1 = 0
            i2 = min(5,nk)
            l2 = i2+1
            a[1] = az
            A = numpy.array(a[2:])
            A = border(A)
            A = org(1, i2, tb, cm, xs, A)
            a[2:] = A[1:]
        else:
            i1 = i2+1
            i2 = i2+6
            if i2 > nk:
               i2 = nk
            l2 = i2-i1+1
            a = org(i1, i2, tb, cm, xs, a)
        st = ''
        for i in range(i1,i2+1):
            st = st + '    %4i    '%(i)
        logger.info('\n bsfn:' + st)

        st = ''
        limit = min(7, nk-cnt+1)
        for i in range(1, limit):
            st = st + '%12.4f'%(a[i])
        logger.info(' coef:' + st)

        cnt = cnt + 6
    return 

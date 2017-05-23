#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from genutil.pymars import LOG, logger
def numprt(n, a):

    i2 = 0
    while i2 < n:
        i1 = i2+1
        i2 = i2+6
        if i2 > n:
            i2 = n
        s1 = ''
        s2 = ''
        for j in range(i1, i2+1):
            s1 = s1 + '    %4i'%(j)    
            s2 = s2 + '%10.4g'%(a[j])
        logger.info(s1+'\n')
        logger.info(s2+'\n')
    return

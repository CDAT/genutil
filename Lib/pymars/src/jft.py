import numpy
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def jft(m, j, tb):
    #no longer needed in checkVariable
    #print 'jft'
    x = numpy.array(abs(tb[2,1:m+1])+.1, dtype = INT_DTYPE)
    (k,) = numpy.where(x == j)
    return (len(k) > 0)
#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def stseed(seed):
    parameters['stseed']['i'] = seed
    return
def rnms(n):
    x = numpy.zeros(shape = n+1, dtype = FLOAT_DTYPE)
    for j in range(1, n+1):
        #self.i = dmod(self.i*16807.d0,2147483647.d0)
        parameters['stseed']['i'] = parameters['stseed']['i']*16807 % 2147483647
        u = parameters['stseed']['i']
        u = u*.465661287e-9
        x[j] = u
    return x

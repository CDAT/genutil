#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
def array(p, n):
    i = numpy.fmod(p,n)
    if i == 0:
    	i=n
    j = (p-i)/n + 1
    return i, j

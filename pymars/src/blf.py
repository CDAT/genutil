#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def blf(l, n, sc):

    if l <= 0:
        bl = numpy.array([0.0]+n*[1.0], FLOAT_DTYPE)
    else:
        #this was changed so that I wasn't moving large
        #amounts of data; see where it's used
        bl = sc#[:,l]
    return bl

#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
def border(x, val=0):
    import numpy
    
    dim = x.shape
    if len(dim) == 2:
        m,n = dim
        size = (m+1)*(n+1)
        y = numpy.zeros(shape = (m+1, n+1), dtype=type(x[0,0]))
        y[1:m+1, 1:n+1] = x
        y[0,:] = val
        y[:,0] = val
    elif len(dim) == 1:
        m = len(x)
        y = numpy.zeros(shape=m+1, dtype=type(x[0]))
        y[1:m+1] = x
        y[0] = val
    return y
def copyback(x):
    m,n = x.shape
    y = []
    for i in range(1, n):
        y = y + x[1:m, i].tolist()
    return numpy.array(y)
if __name__ == '__main__':   
    import numpy, scipy
    x = numpy.array([[1,2,3],[4,5,6]], dtype=float)
    print border(x, val = -scipy.inf)
    
    y = numpy.array([1,2,3.])
    print border(y, val=-1)
    
    
    z = copyback(border(x.T, val=-1))
    print z
    
    u=numpy.zeros(shape=5, dtype=float)
    print border(u)
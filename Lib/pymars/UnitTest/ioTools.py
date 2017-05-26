import numpy
def getData(testName):
    if testName.startswith('Cat'):
        try:
            print 'reading:' + "logit.csv"
            data = numpy.loadtxt("logit.csv", dtype = 'float', skiprows=1, delimiter=',')
            input = data[:,1:].copy()
            output = data[:,0]     
        except:
            raise 'problem reading input file'   
    elif testName.startswith('Cont'):
        try:
            data = numpy.loadtxt("MarsTest.dat", skiprows=1)
            input = data[:,0:-1]
            output = data[:,-1]
        except:
            raise 'problem reading input file'
    elif testName.startswith('Mixed'):
        n = 10
        x = numpy.array(range(n))*.1
        y = numpy.array(5*[0.,1])
        z = numpy.array([2., 1., 0, 2., 1., 0, 2., 1., 0, 2])

        input = cartesian((x,y,z))
        output = input[:,0]**2 + input[:,1]**2 + input[:,2]**2  

    elif testName.startswith('Multi'):
        n = 10
        x = numpy.array([2., 1., 0, 3, 4., 0, 2., 1., 3, 4])
        y = numpy.array(5*[0.,1])
        z = numpy.array([2., 1., 0, 2., 1., 0, 2., 1., 0, 2])

        input = cartesian((x,y,z))
        output = input[:,0]**2 + input[:,1]**2 + input[:,2]**2  
    else:
        input = []
        output = []
    return input, output      
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [numpy.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = numpy.prod([x.size for x in arrays])
    if out is None:
        out = numpy.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = numpy.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

if __name__ == '__main__':
    I,O = getData('Multi')
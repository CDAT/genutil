import string, numpy
class ParseMarsMemory:
    def __init__(self, (fm, im)):
        #Remember when comparing the points into fm & im
        #Python/C counts from 0, while Fortran counts from 1
        #Also Fortran arrays are row major while Python/C
        #arrays are column major
        
        inputs = {}
        inputs['nTrainingPnts'] = im[2]
        inputs['dim'] = im[3]
        inputs['nKnots'] = im[4]
        inputs['maxInteraction'] = im[5]
        self.inputs = inputs
        
        constPntr = im[10]-1
        bfPntr = im[11]-1
        nBFs = self.inputs['nKnots']
        basisFunctions=numpy.array(fm[bfPntr:5*nBFs + 1])
        basisFunctions.shape=(nBFs, 5)
        BFs = [fm[constPntr]]
        count = 0
        for bf in basisFunctions:
            coef = bf[0]
            if coef != 0.0:
                count = count+1
            var = int(abs(bf[1]))
            knot = float(bf[2])
            parent = int(bf[3])
            sign = bf[1]/abs(bf[1])
            sign = int(sign)
            bf = {'var':var, 'knot':knot, 'coef':coef, 'parent': parent, 'sign':sign}
            BFs = BFs + [bf] 
        self.basisFunctions = BFs
        
        outputs = {}
        outputs['nBFused'] = count
        outputs['gcv'] = 'not yet known, its buried in mars'
        self.outputs = outputs
if __name__ == '__main__': 
    from VVutil.ColumnIO import * 
    import mars
    
    inputFile = "MarsTest.dat"
    print "Input File: %s" % inputFile
    n=98
    print "number of observations: ", n
    p=2
    print "number of predictor variables (dimensions): ", p
    nk=15
    print "max number of basis functions: ", nk
    mi=2
    print "max_interaction: ", mi
    m=2
    print "model flag; 1: plot linear mars model, 2: plot cubic mars model: ", m
    #x = numpy.array(range(npts*nparams), numpy.Float32)
    ngc=100
    print "number of raster points for computing curve estimates: ", ngc
    ngs=40
    print "number of raster points on each axis for computing curve estimates: ", ngs
    icx=1
    print "convex hull flag; 0: plot over entire range of arg limits, >0: plot only inside convex hull: ", icx
    

    w = numpy.array(n*[1], numpy.float32)
    #print "weights: ", weights
    lx = numpy.array(p*[1], numpy.int32)
    #print "lx: ", lx

    dataDict = readData(inputFile)
    vars = dataDict['vars']
    data = dataDict['data']

    x = numpy.array([data[vars.index('x1')],data[vars.index('x2')]])
    x = numpy.transpose(x)
    print "x.shape: ", x.shape
    #print "x: ", x

    y = numpy.array([data[vars.index('y')]])
    y = numpy.transpose(y)
    print "y.shape: ", y.shape
    #y = data[vars.index('y')]
    #print "y: ", y
    fm, im = mars.mars(n, p, x, y, w, nk, mi, lx)
    test = ParseMarsMemory((fm, im))
    nk=im[4]
    knots=numpy.array(fm[1:5*nk + 1])
    knots.shape=(nk, 5)
    for i in range(nk):
        print knots[i, 0:5]
    
    for bf in test.basisFunctions:
        print bf
    print test.inputs
    print test.outputs
    b=test.basisFunctions[1]
    print type(b['var'])
    print type(b['sign'])
    print type(b['parent'])
    print type(b['knot'])
    print type(b['coef'])
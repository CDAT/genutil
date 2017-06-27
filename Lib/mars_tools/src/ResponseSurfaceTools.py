import Sampling.LHS, Sampling.MonteCarlo, numpy, math, mars, scipy
from VVutil.ColumnIO import *
from MarsResponseSurface import *
def buildHullMask(samples, box, nbasisFncs = 20, maxInteraction = 4):
    """this function builds a Mars surface on the function
    which is 1 on the data provided and 0 on the boundary
    of the box.  It's purpose is to compute a stochastic
    convex hull for a set of points."""
    
    npts = len(samples[0])
    dim = len(box)

    #create samples in the faces and append to the existing data
    #total number of points per face is approximately the number
    #of interior points.
    nptsPerFace = npts/(2*dim)
    faces = getSamplesInFaces(dim, nptsPerFace, box)
    d=list(samples)
 
    for face in faces:
        for i in range(dim):
            d[i] = d[i] + face[i]
    nfaces = 2*dim

    allData = numpy.array(d, order='C')   
    allData = allData.transpose()
 
    #allData = Numeric.transpose(Numeric.array(d))       
 
    #make z as an array of constant value 1.0
    #and equal zero on the faces
    z=npts*[1.0]
    z = z + nptsPerFace*nfaces*[0.0]
    z = numpy.array(z, order='C')
    nPts = len(z)

    #compute MARS model
    weights = numpy.array(nPts*[1], numpy.float32)
    lx = numpy.array(dim*[1], numpy.int32)
    mars.setdf(1.)
    fm, im = mars.mars(nPts, dim, allData, z, weights, nbasisFncs, maxInteraction, lx)

    #mrs = MarsResponseSurface('mars.log') 
    mrs = MarsResponseSurface((fm, im))
    
    return mrs
def getDigits(nDigits, base):
    """this method computes all of the digits between
    0 and base**nDigits-1 for any specified base"""
    n = base**nDigits
    digits = []
    for i in range(nDigits):
        digits.append([])
    for i in range(n):
        quotient = i
        for j in range(nDigits):
            remainder = math.fmod(quotient, base)
            digits[j].append(remainder)
            quotient = quotient/base
    return digits    

def computeMonteCarloVolume(box, mrs, nsamples, lb, ub):
    """compute an estimate of the volume as a ratio of the
    number of random points inside the region with the volume 
    of the region"""
    
    samples = Sampling.MonteCarlo.genMC(nsamples, box)
    dim = len(box)

    count = 0
    for i in range(nsamples):
        x=[]
        for j in range(dim):
            x = x + [samples[j][i]]
        y = mrs(tuple(x))
        if lb <= y and y <= ub:
            count = count+1
    vol = 1.0
    for (blow, bhgh) in box:
        vol = vol*(bhgh - blow)
    volume = vol*count/nsamples
    return volume
def genBox(dim, low, high):
    box = []
    for i in range(dim):
        box = box +[(low, high)]
    return box
def getMinMax(mrs, data):
    """compute the min & max of the response surface at the original data"""
    upper = -10.
    lower = 10.
    dim = data['nvars']
    npts = data['npts']
    d = data['data']
    for i in range(npts):
        x=[]
        for j in range(dim):
            x = x + [d[j][i]]
        x=tuple(x)
        val = mrs(x)
        upper = max(val, upper)
        lower = min(val, lower)

    return (lower, upper)
def getStatisticalMinMax(data, width = 1.0):
    """compute the statistical min & max of the response surface 
    from the original data"""

    mu = scipy.mean(data)
    sd = scipy.std(data)
    lower = mu-width*sd
    upper = mu+width*sd
    return (lower, upper)
def getConstrainedPnts(samples, mrs, lb, ub):
    """get those points on the respose surface between lb and ub"""
    
    dim = len(samples)
    nsamples = len(samples[0])
    
    pnts = initDict()
    pnts['data']=[]
    for i in range(dim):
        pnts['data'] = pnts['data'] + [[]]
   
    for i in range(nsamples):
        x=[]
        for j in range(dim):
            x = x + [samples[j][i]]
        y = mrs(tuple(x))

        if lb <= y and y <= ub:
            for j in range(dim):
                pnts['data'][j].append(x[j])
        
    pnts['nvars'] = dim
    pnts['npts'] = len(pnts['data'][0])
    return pnts
def getPntsInBox(samples, box):
    """get those points that belong to the box"""
    
    dim = len(samples) 
    nsamples = len(samples[0])
       
    pnts=initDict()
    pnts['data']=[]
    for i in range(dim):
        pnts['data'] = pnts['data'] + [[]]
        
    for i in range(nsamples):
        keep = True
        x=[]
        for j in range(dim):
            (low, hgh) = box[j]
            val = samples[j][i]
            x = x + [val]
            keep = keep and (low <= val and val <= hgh)

        if keep:
            for j in range(dim):
                pnts['data'][j].append(x[j])
        
    pnts['nvars'] = dim
    pnts['npts'] = len(pnts['data'][0])

    return pnts
def getSamplesInFaces(dim, npts, box):
    """get LHS samples in each of the 2*dim faces"""

    faceSamples = []
    for i in range(dim):
        face = []
        for j in range(dim):
            if j != i:
                face = face + [box[j]]
        (low, hgh)=box[i]

        lhsdata = Sampling.LHS.genUniformLHSdata(npts, face)

        f = list(lhsdata)
        f.insert(i, npts*[hgh])
        faceSamples = faceSamples + [f]
        
        f = list(lhsdata)
        f.insert(i, npts*[low])
        faceSamples = faceSamples + [f]

    return faceSamples
def boxVolume(box):
    """compute the volume of the box"""
    vol = 1.0
    for (L,H) in box:
        vol = vol*(H-L)
    return vol
def crossSection(mrs, box, smallBox, index, npts, dump = False):
    """This function computes a cross section of the response surface
    in the index dimension and the midpoint of the other dimensions.
    This is used for diagnostics of the response surface."""
    
    f = []
    dim = len(box)
    slice = []
    for i in range(dim):
        (L, H) = smallBox[i]
        if i != index:
            slice = slice + [(H-L)/2]
    (L, H) = box[index]

    delta = (H - L)/(npts-1)
    for i in range(npts):
        val = L + i*delta
        v=list(slice)
        v.insert(index, val)
        x = tuple(v)
        y = mrs(x)
        if dump:
            print x, y
        f = f + [y]
        
    return f
def buildResponseSurfaceOnRegularGrid(mrs, fm, im, nreq, box):    
    """compute response surface on a regular grid; used for graphics purposes"""
    
    dim = len(box)
    digits = getDigits(dim, nreq)
    for i in range(dim):
        (low, hgh) = box[i]
        delta = (hgh - low)/(nreq - 1)
        digits[i] = low + delta*numpy.array(digits[i])

    RS = []
    for i in range(dim+1):
        RS = RS + [[]]
    for i in range(nreq**dim):
        x = []
        for j in range(dim):
            x = x + [digits[j][i]]
        X = numpy.array(x)
        X.shape= 1,dim
        val = mrs(X, fm, im)
        x = x +[val]
        for j in range(dim+1):
            RS[j] = RS[j] + [x[j]]
    return RS        
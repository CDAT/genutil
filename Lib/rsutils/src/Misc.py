import mars, string
from numpy import *
import string

def reshuffle(ind, x):

    y=numpy.array(x)
    i=0
    for j in ind:
        y[j]=x[i]
        i=i+1
        
    return y
 
def computeMarsModel(npts, nparams, x, y, nbasis_fncs, max_interaction, x_request, nrequests):
    "This is a general purpose method to first compute the mars model\r\n \
    for the input data x vs y and return the model evaluated at the entries\r\n\
    in x_request.  x_request is an nparams*npts**nparams array."
    
    #use piecewise cubic model
    m=2

    weights = array(npts*[1], Float32)
    lx = array(nparams*[1], Int32)
    fm = array(50000, Float32)
    im = array(50000, Int32)
    sp = array(200000, Float32)
    dp = array(200000, Float64)
    mm = array(50000, Int32)
    f = array(nrequests*[0], Float32)

    #compute MARS model
    fm, im, sp, dp, mm = mars.mars(npts, nparams, x, y, weights, nbasis_fncs, max_interaction, lx)
    
    #compute the values of the mars model requested in x_request 
    #note: it may seem that only part of x_request is passed to fmod.  The entire array is 
    #referenced deeper in the mars code.  The only reason part of the array is passed is so 
    #that python stops bitching about an array length.  At this stage there does not seem to 
    #be an elegant solution to this problem.  Once we abandon the garbage fortran code,  this 
    #issue disappears.
    f = mars.fmod(m, nrequests, nparams, x_request, fm, im, sp)
    
    return f
def getFunctions(fn):
    f=open(fn, 'r')
    pairs = []
    for line in f.readlines():
        if line[0] != ";":
            line=string.strip(line)
            line=string.split(line, ",")
            pairs.append((line[0], string.strip(line[1])))
    f.close()    
    return pairs       
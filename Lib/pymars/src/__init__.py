#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
"""
pymars is a pure python version of mars, originally written in fortran.
python INPUTS
mars(n, p, x, y, w, nk, mi, lx) has a number of inputs
input data: 
    x is an nxp dimensional array for the independent variables
    y is an n-dimensional array for the dependent variable
weights:
    w is an n-dimensional array that weights each value of y
knots:
    nk is the number of requested knots 
max interaction:
    mi is an integer that indicates how many of the basis functions 
    can appear in a product
input variable flags:
    lx is a p-dimensional integer array used to direct the algorithm for each variable
    value  indicates
    3      variable has only linear effect
    2      variable is treated as only additively
    1      unrestricted variable
    0      exclude
    -1     unrestricted categorical variable
    -2     categorical variable: same as 2
    NOTE: the only capability tested thus far is value of 1 or -1
pymars directives:
from genutil.pymars import parameters, marsParameters
marsParameters() sets the defaults.
regression type:
    parameters['mars']['il'] = value
    value  indicates
    0      ordinary regression
    1      logistic regression; this assumes the data takes values of 0 or 1.
cross validation:
    parameters['cvmars']['ix'] = value
    value  indicates
    0      turned off
    >1     turned on, the value is used as the number of cross-fold validation
    NOTE: value < 0 is acceptable but untested; it does a validation where the
    the value represents array entry used.
degrees of freedom:
    parameters['mars']['df'] = value   
    df is a parameter that represents how hard the algorithm should work to
    get a fit, specifically in the knot optimization.  The default value is 3.0.
    Smaller values can be used for a better and more time consuming fit. Testing
    only used the default.
other directives are available but remain untested.   
"""
import ctypes
import logging
import numpy

from genutil.pymars.addpar import ADDPAR
from genutil.pymars.logitl import LOGITL
from genutil.pymars.nest import NEST
from genutil.pymars.setint import SETINT

__all__=['addpar', 'anova', 'anoval', 'array', 'atoscl', 'basisFunction', 'bkstp', 'blf', 
        'blf0', 'border', 'catpr', 'ccoll', 'cmrs', 'coefpr', 'coll', 'collc', 'collf',
        'cptb', 'csp', 'cubic', 'cue', 'cvmars', 'cvmod', 'efp', 'elg', 'exch', 'fmod',
        'fmrs', 'generalizedCrossValidation', 'holl', 'ibfext', 'icat', 'icf', 'ieq', 
        'ieqbf', 'initialize', 'jf', 'jft', 'jfv', 'jfvc', 'knts', 'lcm', 'logitl', 
        'lsf', 'lsf1', 'lstsqr', 'mars', 'marsgo', 'mnspan', 'mrsgo1', 'mrsgo2', 'mrsgo3',
        'nest', 'newb', 'nnord', 'nord', 'nordc', 'numprt', 'ordpr', 'org', 'orgpc',
        'orgpl', 'phi', 'pr', 'purcat', 'que', 'rspnpr', 'sclato', 'scpc', 'setint',
        'side', 'spofa', 'sposl', 'sscp', 'stseed', 'sweep', 'update', 'varf', 'varimp',
        'varz', 'vp', 'LOG', 'TIME', 'ARRAY_SIZE', 'FLOAT_DTYPE', 'INT_DTYPE', 'logger',
        'debug', 'parameters', 'SETINT0', 'NEST0', 'ADDPAR0', 'LOGITL0', 'cmars']
def marsParameters():
    #parameters contains all of the mars parameters that have a data statement
    #it is a dictionary of dictionaries; one for routine, the parameters are the keys
  
    parameters['mars'] = {'ms':0, 'df':3.0, 'il':0, 'fv':0.0, 'ic':0, 'z00001':0}
    parameters['marsgo']= {'alr':1.e-7, 'eps':1.e-4, 'big':9.9e30, 'fln':-1.0, 'nmin':5, 
                           'alf':0.0500000007, 'vcst':[0, 1.0,.666667,.333333]}
    parameters['cvmars']= {'ix':0, 'eps':1.E-6, 'big':9E30, 'dfs':0.0, 'cvm':0.0, 'im':0}
    parameters['lsf'] = {'eps':1.e-05}
    parameters['stseed'] = {'i':987654321}
    parameters['fmrs'] = {'ifg':0}
    parameters['cmrs'] = {'ifg':0}
    parameters['elg'] = {'ic':0}
    #parameters['cptb'] = {'ub':numpy.zeros(shape=(5+1, nk+1), dtype=numpy.float64)}
    return 

class simplelogging():
    """This logger is being used because I could not get python loggong to work right.
    It keeps messing up the first pymars logger when executed."""
    def __init__(self,fn):
        self.buf = []
        self.fn = fn
        return
    def setlog(self,fn):
        self.buf = []
        self.fn = fn
        return
    def info(self, s):
        self.buf = self.buf + [s + '\n']
        return
    def close(self):
        file = open(self.fn, 'w')
        file.writelines(self.buf)
        file.close()
        return

def setlogging(fn):
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(fn, 'w')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)  
    return fh
def closefh():
    fh.close()
    return
def xvalid(ix):
    parameters['cvmars']['ix'] = ix
    return
def setdf(x):
    parameters['mars']['df'] = x
    return
def speed(i):
    lque = [9999, 20, 20, 10, 5]
    freq = [9.e30, 9.e30, 0.2, 0.2, 0.2]
    j = i
    if i < 0: 
        j=0
    if i > 4:
        j=4
    ADDPAR0.setmpr(lque[j])
    ADDPAR0.setfrq(freq[j])
    return
def timeTotals(PRINT=True):
    totals={}
    for key in TIME.keys():
        total = 0.0
        for s,e in TIME[key]:
            total += e-s
        totals[key] = total
        if PRINT:
            print 'total time in '+ key + ' is ' + str(total)
    return totals
#define some package parameters
TIME = {}
ARRAY_SIZE = 10000           #a large number that may not be large enough
FLOAT_DTYPE = numpy.float64  #change to higher precision if necessary
INT_DTYPE = numpy.int64      #change to higher precision if necessary
LOG = True
#logger = logging.getLogger('pymarsLog' )
logger = simplelogging('pymarsLog')
#logging.basicConfig(filename= 'pymars.log', filemode='w')
#fh = logging.FileHandler('pymars.log', 'w')
#logger.addHandler(fh)

debug = logging.getLogger('debugLog' )
debug.setLevel(logging.DEBUG)
dfh = logging.FileHandler('debug.log', 'w')
dfh.setLevel(logging.DEBUG)
debug.addHandler(dfh)

parameters = {}
marsParameters()

SETINT0 = SETINT()
NEST0 = NEST()
ADDPAR0 = ADDPAR()
LOGITL0 = LOGITL()


try:
    import mpi
    #from NodePool import *
    #np = NodePool()
    MPI = True
except:
    MPI = False
MPI = False

#SITE_PACKAGES = "/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/"
#SITE_PACKAGES = "/Users/mcenerney1/anaconda/envs/2.8/lib/python2.7/site-packages/"
#cmars = Ctypes.cdll.LoadLibrary(SITE_PACKAGES + "cmars.so")
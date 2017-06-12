import logging
import numpy


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

class mylogging():
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
    return

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
TIME = None
ARRAY_SIZE = 10000           #a large number that may not be large enough
FLOAT_DTYPE = numpy.float64  #change to higher precision if necessary
INT_DTYPE = numpy.int64      #change to higher precision if necessary
LOG = True
logger = mylogging('marsLog')
#logging.basicConfig(filename= 'pymars.log', filemode='w')
#fh = logging.FileHandler('pymars.log', 'w')
#logger.addHandler(fh)

debug = logging.getLogger('debugLog')
debug.setLevel(logging.INFO)
dfh = logging.FileHandler('debug.log', 'w')
dfh.setLevel(logging.INFO)
debug.addHandler(dfh)

parameters = {}
marsParameters()
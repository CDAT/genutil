import os, sys, numpy, scipy, time, pdb
import genutil.mars, genutil.pymars
from genutil.pymars import LOG, logger, parameters, marsParameters
from genutil.mars import logit as F_logit
from genutil.mars import setdf as F_setdf
from genutil.mars import xvalid as F_xvalid
from genutil.mars import mars as F_mars
from genutil.mars import fmod as F_fmod
from genutil.mars import setlog as F_setlog

from genutil.pymars.mars import mars as py_mars
from genutil.pymars.fmod import fmod as py_fmod
from genutil.pymars.border import border
from errorTools import *
from ioTools import *
if genutil.pymars.MPI:
    import mpi

class RunMarsAndPymars():
    def __init__(self):
        return
    def execute(self, testName, x, y, nk, mi, m, lx, logisticRegression, crossValidation,
                logfile, Flogfile):
        genutil.pymars.TIME = {}
        
        logger.setlog(logfile)
        F_setlog(Flogfile)
            
        logger.info(testName)
        logger.info('logisticRegression='+repr(logisticRegression))
        logger.info('crossValidation='+repr(crossValidation))
        logger.info('mars output')
        
        #set the mars parameters
        n,p = x.shape
        print "x.shape: ", x.shape
        print "y.shape: ", y.shape
        
        w = numpy.array(n*[1], dtype=numpy.float64)
        if lx == None:
            lx = p*[1]
        lx = numpy.array(lx, dtype=numpy.int32)   
        print 'lx=',lx   
        #logistic regression
        F_logit(logisticRegression)
        #cross validation
        F_xvalid(crossValidation)
        
        #run Friedman Fortran mars
        start = time.time()
        fm, im = F_mars(n, p, x, y, w, nk, mi, lx)
        end = time.time()
        F_time = end-start
        #evaluate Friedman response surface at input variable values
        f = F_fmod(m, n, p, x, fm, im)
        #print _fmod.__doc__
        print 'error in evaluation of fortran response surface on input data'    
        computeErrors('fortran', y, f, n)
    
        #run python mars        
        #put -infinity in dummy location for future sorting in these parameters
        logger.info('pymars output')
        x = border(x, val = -scipy.inf)
        y = border(y, val = -scipy.inf)
        w = border(w)
        lx = border(lx)
            
        #logistic regression
        parameters['mars']['il'] = logisticRegression
        #cross validation
        parameters['cvmars']['ix'] = crossValidation
        
        #there is some sort of overflow with exp with logistic regression
        numpy.seterr(over='ignore')
        start = time.time()
        #try:
        x, az, tb, cm, kp, kv, lp, lv, bz, tc = py_mars(n, p, x, y, w, nk, mi, lx)
        #except:
        #    pass
        end = time.time()
        PY_time = end-start
        
        print "Returned from python mars"
        logger.info("Returned from python mars")      
                  
        #evaluate python response surface at input variable values
        pf = py_fmod(m, n, nk, x, az, tb, cm, kp, kv, lp, lv, bz, tc)

        print 'error in evaluation of python response surface on input data'
        computeErrors('Pymars', y[1:], pf[1:], n)
        logger.close()
        
        print 'fortran time=', F_time
        print 'python time=', PY_time
        
        return

                
if __name__ == '__main__':
    master = 0
    if genutil.pymars.MPI:
        master = mpi.rank

    if master == 0:
        marsParameters()
        runMarsAndPymars = RunMarsAndPymars()
        args = sys.argv[1:]
        print args
        testName = args[0]
        input, output = getData(testName)
        nk = int(args[1])
        mi = int(args[2])
        m  = int(args[3])
        logisticRegression = int(args[4])
        crossValidation = int(args[5])
        logfile = args[6]
        Flogfile = args[7]
        lx = eval(args[8])  #do this to convert string to a list

        runMarsAndPymars.execute(testName, input, output, nk, mi, m, lx,
                                 logisticRegression, crossValidation,
                                 logfile, Flogfile)
        #shutdown processors
        if genutil.pymars.MPI:
            mpi.scatter(None)

        #if pymars.MPI:
            #print 'node count is ', pymars.np.nodeCount
            #pymars.np.shutdown()
            #mpi.scatter(None)
    else:
        from genutil.pymars.mrsgo1 import computeProposedHockeyStick, makePreliminaryDecision
        #title = 'slave #'+ str(mpi.rank) + ' time ='
        #print 'from processor ', mpi.rank
        genutil.pymars.TIME = {}

        while True:
            input, status = mpi.recv()
            if input == None:
                break
            (variable, X, MM,
            parent, n, p, y, yb, w, sw, lx, tb, cm, db,
            bfIndex, bfIndexLast, parentEval, partialEval,
            ms, fv, nBasisFunctions, nep,  rsq, kcp0, catKnots,
            mn, me, mel, ja, lbf,  nnt, kr,   DY,
            fln, eps, big, nmin, alf) = input
    
            (CM, nc, ict,   RSQ,  ja, fvr,  prelimDecision, defaultHS, proposedHS) = \
                    computeProposedHockeyStick(input)
    
            hockeyStickState = makePreliminaryDecision(prelimDecision, tb,  nc, ict, bfIndex, fln, fvr, sw, 
                                                       RSQ,  defaultHS, proposedHS, CM, variable, ja)
            mpi.send(hockeyStickState, 0)

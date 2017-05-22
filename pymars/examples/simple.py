#simple example

import numpy, scipy, time, random, pymars
from pymars import LOG, logger, parameters, marsParameters, FLOAT_DTYPE, timeTotals
from pymars.border import border
from pymars.mars import mars as py_mars
from pymars.fmod import fmod as py_fmod

#marsParameters()
testName = 'simple'
mi = 1
m = 1
logisticRegression = 0
crossValidation = 0
logfile = 'simple_pymars.log'
Flogfile = 'simple_mars.log'
logger.setlog(logfile)
            
logger.info(testName)
logger.info('logisticRegression='+repr(logisticRegression))
logger.info('crossValidation='+repr(crossValidation))
logger.info('mars output')

#setup timing study
pymars.TIME = {}
pymars.TIME['marsgo2'] = []
pymars.TIME['update'] = []

p = 2
n = 100
nk = 10

x=[]
y=[]
for i in range(n):
    x = x + [random.random()]
    y = y + [random.random()]
x = numpy.array(x)*.01
y = numpy.array(y)*.01 
z = x**2 + y**2

x = numpy.array([x,y]).T
w = numpy.array(n*[1], dtype=numpy.float64)
lx = numpy.array(p*[1], dtype=numpy.int32)    
try:
    import mars
    from mars import logit as F_logit
    from mars import setdf as F_setdf
    from mars import xvalid as F_xvalid
    from mars import mars as F_mars
    from mars import fmod as F_fmod
    from mars import setlog as F_setlog

    F_setlog(Flogfile)
    #logistic regression
    F_logit(logisticRegression)
    #cross validation
    F_xvalid(crossValidation)
        
    #run Friedman Fortran mars
    start = time.time()
    fm, im = F_mars(n, p, x, z, w, nk, mi, lx)
    end = time.time()

    #evaluate Friedman response surface at input variable values
    f = F_fmod(m, n, p, x, fm, im)
    print F_fmod.__doc__
    print 'fortran time=', end-start
except:
    print 'mars is not installed'
   
#run python mars        
#put -infinity in dummy location for future sorting inside the algorithm
logger.info('pymars output')
x = border(x, val = -scipy.inf)
z = border(z, val = -scipy.inf)
w = border(w)
lx = border(lx)
            
#logistic regression
parameters['mars']['il'] = logisticRegression
#cross validation
parameters['cvmars']['ix'] = crossValidation
   
start = time.time()
x, az, tb, cm, kp, kv, lp, lv, bz, tc = py_mars(n, p, x, z, w, nk, mi, lx)
end = time.time()
PY_time = end-start
        
print "Returned from python mars"
logger.info("Returned from python mars")      
                  
#evaluate python response surface at input variable values
pf = py_fmod(m, n, nk, x, az, tb, cm, kp, kv, lp, lv, bz, tc)
logger.close()
print 'python time=', PY_time
timeTotals()
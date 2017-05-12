#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy, pymars
from pymars import LOG, logger, debug, parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from pymars import ADDPAR0, NEST0, SETINT0
from pymars.stseed import rnms
from pymars.marsgo import buildResponseSurface
from pymars.cvmod import cvmod
from pymars.border import border
from pymars.generalizedCrossValidation import generalizedCrossValidation as GCV

def crossValidation(n, p, x, y, w, nk, ms, df, fv, mi, lx, xm, xs, tb, cm, mm):
    """ Cross validation of mars"""
    
    wt = numpy.zeros(shape = (n,2), dtype = FLOAT_DTYPE, order='F')
    wt = border(wt)
        
    cv = numpy.zeros(shape = (nk,4), dtype = FLOAT_DTYPE, order='F')
    cv = border(cv)
        
    ix = parameters['cvmars']['ix']
    if LOG:
        logger.info('\n cvmars output')
        logger.info('\n sample reuse to estimate df:')
    if ix > 0:
        nr = ix
        nd = nr
        logger.info(' %3i - fold cross-validation.'%(ix))
    else:
        nr = 1
        nd = -ix
        logger.info('independent test set - every %4i observations.\n'%(nd))
        
    wt[1:,1] = w[1:]
    wt[1:,2] = range(1,n+1)
        
    #setup the split of data for cross validation
    for i in range(1, n+1):
        r = rnms(1)[1]
        k = (n-i+1)*r+i
        t = wt[i,2]
        wt[i,2] = wt[k,2]
        wt[k,2] = t

    cv0 = 0.0
    sw = cv0
    wn = sw
    yv = wn
    fc = yv
    for ir in range(1, nr+1):
        i = ir
        while i <= n:
            wt[int(wt[i,2]+.1),1] = 0.0
            i = i + nd
        #dummy assignment
        az = 0.0
        tb = numpy.zeros(shape = (5+1, nk+1), dtype = FLOAT_DTYPE)
        az, tb, cm, cvStuff = buildResponseSurface(n, p, x, y, wt[:,1], nk, ms, df, fv, 
                                                   mi, lx, xm, xs, az, tb, cm, mm)
        cvStuff = cvStuff[1:,1:].T.flatten()
        cvStuff = border(cvStuff)
        yv1 = cvStuff[3]
        yv = yv + yv1
        wn1 = cvStuff[2]
        wn = wn + wn1
        fc = fc + cvStuff[1]
        mk = int(cvStuff[(nk+1)**2+4] + .1)
        i = ir
        while i <= n:
            k = int(wt[i,2]+.1)
            wt[k,1] = w[k]
            sw = sw + w[k]
            CV = numpy.array(cv[1:, 3:])
            CV = border(CV)
            cv0, CV = cvmod(k, n, x, y, w, nk, mk, tb, cm, cvStuff, cv0, CV) #cv[1,3])
            cv[1:, 3:] = CV[1:, 1:]
            i = i + nd
        for m in range(1, nk+1):
            am = cvStuff[m+4]
            cv[m,2] = cv[m,2] + am
            am1 = yv1
            if m > 1: 
                am1 = cvStuff[m+3]
            if am1/yv1 > parameters['cvmars']['eps']:
                r = numpy.sqrt(am/am1)
            else:
       			r = 1.0
            cv[m,1] = cv[m,1] + ((wn1-1.0)*(1.0-r)/(m-r*(m-1))-1.0)/cvStuff[1]
    cv[1:nk+1,1] = cv[1:nk+1,1]/nr
    cv[1:nk+1,2] = cv[1:nk+1,2]/nr
    cv[1:nk+1,3] = cv[1:nk+1,3]/sw

    fc = fc/nr
    yv = yv/nr
    wn = wn/nr
    cv0 = cv0/sw

    if LOG: 
        logger.info('\n  #bsfns      df         asr           gcv           cv')
    parameters['cvmars']['im'] = 0
    parameters['cvmars']['cvm'] = cv0
    dmx = -parameters['cvmars']['big']
    cvl = cv[nk,1]
    m = nk

    while m >= 1:
        if cv[m,1] > dmx:
            dmx = cv[m,1]
            dfu = 0.5*(cvl + cv[m,1])
            cvl = cv[m,1]
            if cv[m,3] <= parameters['cvmars']['cvm']:
          		parameters['cvmars']['cvm'] = cv[m,3]
          		df = dfu
          		parameters['cvmars']['im'] = m
            gcv = GCV(cv[m,2], (dfu*fc+1.0)*m+1.0, wn)
            if LOG: 
    			logger.info(' %5i %10.2f %12.4g %12.4g %12.4g' %(m,dfu,cv[m,2],gcv,cv[m,3]))
        m = m-1
    if cv0 <= parameters['cvmars']['cvm']:
        parameters['cvmars']['cvm'] = cv0
        df = dmx
        parameters['cvmars']['im'] = 0
    parameters['cvmars']['dfs'] = df
    gcv = GCV(yv, 1.0, wn)
    if LOG:
        logger.info(' %5i %10.2f %12.4g %12.4g %12.4g' %(0,dmx,yv,gcv,cv0))
        logger.info('estimated optimal df( %3i ) = %7.2f with (estimated) pse = %12.4g'
                        %(parameters['cvmars']['im'], df, parameters['cvmars']['cvm']))
    return tb, cm, df

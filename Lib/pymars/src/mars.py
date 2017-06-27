import numpy, time, genutil.pymars, pdb
from genutil.pymars import LOG, logger, debug, parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from genutil.pymars import ADDPAR0, NEST0, SETINT0, LOGITL0
from .rspnpr import rspnpr
from .ordpr import ordpr
from .atoscl import atoscl
from .catpr import catpr
from .cvmars import crossValidation
from .sclato import sclato
from .orgpc import orgpc
from .orgpl import orgpl
from .varimp import varimp
from .cubic import cubic
from .ccoll import ccoll
from .anova import anova
from .marsgo import buildResponseSurface
from .fmrs import fmrs
from .anoval import anoval
from .cmrs import cmrs
from .generalizedCrossValidation import generalizedCrossValidation as GCV
from .generalizedCrossValidation import logisticGCV
from .initialize import init

def mars(n, p, x, y, w, nk, mi, lx):
        
    #initialization
    az = 0.0
    tb = numpy.zeros(shape = (5+1, nk+1), dtype = FLOAT_DTYPE)
    cm = numpy.zeros(shape = ARRAY_SIZE, dtype = FLOAT_DTYPE)
    kp = numpy.zeros(shape = (5+1, nk+2), dtype = INT_DTYPE)
    kv = numpy.zeros(shape = (2+1, nk*mi+1), dtype = INT_DTYPE)
    lp = numpy.zeros(shape = (3+1, nk+3), dtype = INT_DTYPE)
    d2 = 16 + 5*nk + 2*nk*mi + 3*(nk+2) + nk*mi-1
    lv = numpy.zeros(shape = ARRAY_SIZE, dtype = INT_DTYPE)    
    bz = 0.0 
    tc = numpy.zeros(shape = nk*(5*mi+1)+1, dtype = FLOAT_DTYPE)


    ms = parameters['mars']['ms']
    df = parameters['mars']['df']
    il = parameters['mars']['il']
    fv = parameters['mars']['fv']
    ic = parameters['mars']['ic']
       
    if LOG:
        logger.info(' MARS modeling, version 4.0 (7/13/09)')
        logger.info( ' input parameters (see doc.):\n ' + \
                    '    n    p    nk   ms    mi   df     il   fv    ic \n' + \
                    '  %5i %3i %4i %4i  %4i   %4.1f    %i    %3.1f   %2i' \
                        %(n,p,nk,ms,mi,df,il,fv,ic)
                     )

        logger.info('\n predictor variable flags:\n')

        predVarStr  = " var:  "
        predFlagStr = " flag: "
        for i in range(1, p+1):
            predVarStr  = predVarStr  + "  " + str(i)
            predFlagStr = predFlagStr + "  " + str(lx[i])
        logger.info(predVarStr + "\n")
        logger.info(predFlagStr + "\n")

    SETINT0.intlst()
    NEST0.nstlst()

    rspnpr(il, n, y, w)  
    mm = numpy.zeros(shape = (n+1, p+1), dtype = INT_DTYPE)
    for j in range(1, p+1):
        mm[:,j] = x[:,j].argsort()
    ordpr(n, p, x, lx, mm)

    x, xm, xs, cm = atoscl(n, p, w, x, lx, mm, cm)

    catpr(n, p, x, cm)
    NEST0.oknest(p, lx, cm)
    
    if parameters['cvmars']['ix'] not in [0,1]:
        tb, cm, df = crossValidation(n, p, x, y, w, nk, ms, df, fv, mi, lx, xm, xs, tb, cm, mm)

    az, tb, cm, cvStuff = buildResponseSurface(n, p, x, y, w, nk, ms, df, fv, 
                                               mi, lx, xm, xs, az, tb, cm, mm)

    if il >= 1:
        az, tb, junk = LOGITL0.logitl(n, x, y, w, nk, il, az, tb, cm)
        
        if LOG:
            sw, wn, junk1, junk2 = init(w, y)
            (k,) = numpy.where(tb[1,1:] != 0.0)
            k = k+1
            ef = 1.0 + tb[5,k].sum()
            ef = GCV(1.0, ef, wn)

            Y = fmrs(n, x, nk, az, tb, cm)
            s, t = logisticGCV(n, ef, w, sw, y, Y)
            logger.info('\n piecewise-linear logistic gcv = %12.4g   ave var=%12.4g'%(s, t))

    if LOG: 
        if il == 0:
            anova(n, x, y, w, nk, tb, cm)

        if il > 0:  
            anoval(n, x, y, w, nk, il, az, tb, cm, lp, lv)

    kp, kv, lp, lv = ccoll(n, p, mi, nk, tb, cm, kp, kv, lp, lv)

    tb, cm, kp, kv, lp, lv, bz, tc = \
        cubic(n, p, x, y, w, nk, tb, cm, kp, kv, lp, lv, bz, tc)

    if il >= 1: 
        bz, tb, tc = \
            LOGITL0.logitc(n, x, y, w, nk, il, cm, tb, kp, kv, lp, lv, bz, tc)
        if LOG:
            Y = cmrs(n, x, cm, kp, kv, lp, lv, bz, tc)
            s, t = logisticGCV(n, ef, w, sw, y, Y)
            logger.info(' piecewise-cubic logistic gcv = %12.4g   ave var =%12.4g'%(s,t))

    variableImportance = varimp(n, p, x, y, w, nk, il, az, tb, cm)

    tb = orgpl(xm, xs, nk, tb, cm)
    
    tc = orgpc(xm, xs, lp, lv, tc)

    x = sclato(n, p, x, xm, xs, cm)
 
    return x, az, tb, cm, kp, kv, lp, lv, bz, tc

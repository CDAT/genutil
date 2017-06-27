import numpy
from genutil.pymars import LOG, logger, parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from .setint import SETINT
from .nest import NEST
from .rspnpr import rspnpr
from .ordpr import ordpr
from .atoscl import atoscl
from .catpr import catpr
from .cvmars import cvmars
from .sclato import sclato
from .orgpc import orgpc
from .orgpl import orgpl
from .varimp import varimp
from .cubic import cubic
from .ccoll import ccoll
from .anova import anova
from .addpar import ADDPAR
from .logitl import LOGITL
from .fmrs import fmrs
from .anoval import anoval
from .cmrs import cmrs
from .generalizedCrossValidation import generalizedCrossValidation as GCV
from .generalizedCrossValidation import logisticGCV
from .initialize import init
def mars1 (n, p, x, y, w, nk, mi, lx,
           az, tb, cm, kp, kv, lp, lv,
           bz, tc):
    """this code was moved into mars"""
    
    ms = parameters['mars1']['ms']
    df = parameters['mars1']['df']
    il = parameters['mars1']['il']
    fv = parameters['mars1']['fv']
    it = parameters['mars1']['it']
    ic = parameters['mars1']['ic']
       
    if LOG:
        logger.info(' MARS modeling, version 4.0 (7/13/09)')
        logger.info( ' input parameters (see doc.):\n ' + \
                    '    n    p    nk   ms    mi   df     il   fv    ic \n' + \
                    '  %5i %3i %4i %4i  %4i   %4.1f    %i    %3.1f   %2i' \
                        %(n,p,nk,ms,mi,df,il,fv,ic)
                     )

        logger.info(' predictor variable flags:\n')

        predVarStr  = " var:  "
        predFlagStr = " flag: "
        for i in range(1, p+1):
            predVarStr  = predVarStr  + "  " + str(i)
            predFlagStr = predFlagStr + "  " + str(lx[i])
        logger.info(predVarStr + "\n")
        logger.info(predFlagStr + "\n")

    SETINT0 = SETINT()
    SETINT0.intlst()
    NEST0 = NEST()
    NEST0.nstlst()
    ADDPAR0 = ADDPAR()
    
    rspnpr(il, n, y, w)  
    mm = numpy.zeros(shape = (n+1, p+1), dtype = INT_DTYPE)
    for j in range(1, p+1):
        mm[:,j] = x[:,j].argsort()
    ordpr(n, p, x, lx, mm)

    x, xm, xs, cm = atoscl(n, p, w, x, lx, mm, cm)

    catpr(n, p, x, cm)
    NEST0.oknest(p, lx, cm)
    
    if parameters['cvmars']['ix'] not in [0,1]:
        tb, cm, df = cvmars(ADDPAR0, NEST0, SETINT0,
                            n, p, x, y, w, nk, ms, df, fv, mi, lx, 
                            xm, xs, tb, cm, mm)

    az, tb, cm, cvStuff = marsgo(ADDPAR0, NEST0, SETINT0,
                                 n, p, x, y, w, nk, ms, df, fv, mi, lx, 
                                 xm, xs, az, tb, cm, mm)
    if il >= 1:
        LOGITL0 = LOGITL()
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
        bz, tb, tc = LOGITL0.logitc(n, x, y, w, nk, il, cm, tb, 
                                    kp, kv, lp, lv, bz, tc)
        if LOG:
            Y = cmrs(n, x, cm, kp, kv, lp, lv, bz, tc)
            s, t = logisticGCV(n, ef, w, sw, y, Y)
            logger.info(' piecewise-cubic logistic gcv = %12.4g   ave var =%12.4g'%(s,t))

    variableImportance = varimp(n, p, x, y, w, nk, il, az, tb, cm)

    tb = orgpl(xm, xs, nk, tb, cm)
    
    tc = orgpc(xm, xs, lp, lv, tc)

    x = sclato(n, p, x, xm, xs, cm)

    return x, az, tb, cm, kp, kv, lp, lv, bz, tc
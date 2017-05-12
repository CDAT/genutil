#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy, logging, time, pymars
from pymars import LOG, logger, debug, parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
from pymars.border import *
from pymars.nnord import nnord
from pymars.blf import blf
from pymars.blf0 import blf0
from pymars.jf import jf
from pymars.elg import elg
from pymars.mnspan import mnspan
from pymars.sscp import sscp
from pymars.lsf1 import lsf1
from pymars.coefpr import coefpr
from pymars.update import update
from pymars.bkstp import bkstp
from pymars.mrsgo1 import findBestChildVariable
from pymars.array import array
from pymars.phi import phi
from pymars.holl import holl
from pymars.generalizedCrossValidation import generalizedCrossValidation as GCV
from pymars.initialize import init

def buildResponseSurface(n, p, x, y, w, nk, ms, df, fv, mi, lx, xm, xs, az, tb, cm, mm):
    """ This function is where the real computations take place """

    if LOG:
        logger.info('\n forward stepwise knot placement:\n')
        logger.info('  basfn(s)    gcv      #indbsfns  #efprms   variable      knot            parent')

    #nextVar is the next variable to work on; formerly jq
    nextVar = 0
    cvStuff= numpy.zeros(shape=(n+1, ARRAY_SIZE), dtype=FLOAT_DTYPE)
    #newBF is the new basis function; formerly tx
    newBF = numpy.zeros(shape=6, dtype=FLOAT_DTYPE)
    #evaluations keeps track of the evaluation of each coordinate
    #with a new basis function
    evaluations = numpy.zeros(shape=(n+1, nk+1), dtype=FLOAT_DTYPE)
    DY = numpy.zeros(shape=nk+1, dtype=FLOAT_DTYPE)
    db = numpy.zeros(shape = (n+1, nk+2), dtype = FLOAT_DTYPE)

    #initialize hockey sticks to infinity
    tb = numpy.array(6*(nk+1)*[numpy.inf], dtype = FLOAT_DTYPE)
    tb.shape = (5+1, nk+1)
        
    nop, jas, nc, k1, k, nep, kl, kcp0 = 8*[0]
    df1, ssq,  txt, xt,  tx1 = 5*[0.0]

    #nk is the number of requested knots
    nBasisFunctions = nk

    for j in range(1, p+1):
        if lx[j] != 0:
            if x[mm[1,j],j] < x[mm[n,j],j]:
                nep = nep + 1
                cst = parameters['marsgo']['vcst'][abs(lx[j])]
                if mi == 1: 
                    cst = min(cst, parameters['marsgo']['vcst'][2])
                df1 = df1 + cst
    if nep == 0:
        if LOG: 
            logger.info('('' no predictor variables.'')')
        raise 'stop'
    if nep == 1: 
        df1 = parameters['marsgo']['vcst'][3]
    cfac = df1/nep
    df1 = df*cfac

    sw, wn, yb, s = init(w, y)
    i=n+1

    yv = s/sw
    tcst = 1.0
    tcmx = wn - df1*parameters['marsgo']['vcst'][1] - 2.0
    
    if cm[1] > 0.0:
        i = 2
        while i <= 2*p:
            if cm[i] > 0.0: 
                kcp0 = int(cm[i+1]+.1)
            i = i+2
    
    #bfIndex is the basis function index; formerly m
    bfIndex = 0
    mtot = bfIndex
    txm = GCV(yv, 1.0, wn)
    #rsq is initialized as the error of using the mean as the response function
    #this error decreases as the algorithm progresses. This means subtract.
    #rsq/sw appears in the numerator of the GCV, I think it means residual square
    rsq = yv*sw
    kr = 0
    nopt = 0
    
    if LOG:
        logger.info('  %3i   %12.4g    %5.1f    %5.1f' %(bfIndex, txm, 0.0, 1.0))
    if parameters['marsgo']['fln'] < 0.0: 
        parameters['marsgo']['fln'] = 1.0 + 4.0/wn
    pymars.ADDPAR0.addpar(0)
    parent = 0
    #main loop
    START = time.time()

    while bfIndex < nBasisFunctions and tcst < tcmx:
        #debug.info('bfIndex = '+str(bfIndex))
        nopt = nopt+1
        pymars.ADDPAR0.itrpar(nopt)
        bfIndexLast = bfIndex
        bfIndex = bfIndex+1
        txi = parameters['marsgo']['big']
        kcp = kcp0
        asq0 = rsq/sw
        parent, nextVar = pymars.ADDPAR0.nxtpar(parent, nextVar)
        #debug.info('before buildNextBasisFunction '+repr(tb[2:5, bfIndex]))

        (tb, cm, mm, evaluations, db, bfIndex, bfIndexLast, 
        ms, fv, newBF, txl, tx1, txi, ssq, rsq, kcp0, 
        nop, nextVar, jas, kr, nc, kcp, k1, DY )= \
        buildNextBasisFunction(parent, nextVar,
                               n, p, x, y, yb, w, sw, nk, lx, tb, cm, mm, mi, 
                               db, bfIndex, bfIndexLast, evaluations, 
                               ms, fv, newBF, tx1, txi, 
                               nBasisFunctions, nep, ssq, rsq, asq0,
                               kcp0, nop, jas,
                               kr, nc, kcp, k1, DY, 
                               parameters['marsgo']['fln'], 
                               parameters['marsgo']['eps'], 
                               parameters['marsgo']['big'],
                               parameters['marsgo']['nmin'], 
                               parameters['marsgo']['alf'])

        variable = int(newBF[2]+.1)
        pymars.ADDPAR0.selpar(int(newBF[4]+.1))
        if cm[2*variable] > 0.:
            nc = int(cm[2*variable+1]+.1) - int(cm[2*variable]+.1) + 1
            kcp0 = kcp0 + nc
        #jas > 0 not yet debugged
        if jas > 0: 
            print 'in marsgo, jas=',jas
            tb, cm, kcp0, k1, kr, newBF, bfIndex, evaluations, 
            jn, rsq, db, DY, tcst, fjn, fkr, gcv = \
                processNestedData(n, x, y, w, sw, yb, tb, cm, jas, jn, 
                                  kcp0, kcp, kr, rsq, nBasisFunctions, 
                                  newBF, bfIndex, evaluations, partialEval, db, DY)
        if bfIndex <= nBasisFunctions:
            INDEX = int(newBF[4]+.1)
            partialEval = blf(INDEX, n, evaluations[:,INDEX])
            tb, evaluations, k1, kr, rsq, db, DY = \
                addBasisFunction(n, x, y, w, sw, yb, tb, cm, 
                                 evaluations, partialEval, k1, kr, rsq, 
                                 bfIndex, newBF, db, DY)
            if (bfIndex < nBasisFunctions and 
                (cm[2*variable] > 0.0  or newBF[3] > x[mm[1,variable],variable])):
                bfIndex = bfIndex + 1
                tb, evaluations, k1, kr, rsq, db, DY = \
                    addReflectedBasisFunction(n, x, y, w, sw, yb, tb, cm, variable, 
                                              evaluations, partialEval, k1, kr, rsq, 
                                              bfIndex, newBF, db, DY)

            tcst = nopt*df1 + kr + 1.0
            #at this stage the basis function has been decided; print it
            if LOG: 
                printBasisFunction(xm, xs, bfIndex, mtot, kr, tcst, cm, newBF, rsq, sw, wn)
            mtot = bfIndex
    #end main loop
    END = time.time()
    if 'marsgo' in pymars.TIME.keys(): pymars.TIME['marsgo'] += [(START, END)]         
    nBasisFunctions = min(bfIndex,nBasisFunctions) #actual number of basis functions found
    tb[1, nBasisFunctions+1:nk+1] = 0.0
    D, XB = sscp(n, nBasisFunctions, evaluations, y, w, yb, yv, sw)
    START = time.time()
    D, b, A, a, DP = lsf1(D, nBasisFunctions, XB, yb, parameters['marsgo']['alr'])
    END = time.time()
    if 'marsgo_lsf1' in pymars.TIME.keys(): pymars.TIME['marsgo_lsf1'] += [(START, END)]        
    (k,) = numpy.where(A[1:] != 0.0)
    nli = len(k)  #number of non-trivial basis functions
    df1 = df1*nopt + nli
    tcst = df1 + 1.0
    df1 = df1/nli
    tb[5,1:] = df1

    asm = GCV(b/sw, tcst, wn)
    tcsts = tcst
    az = a
    tb = copyAtoTB(A, tb)

    if parameters['cvmars']['ix'] != 0:
        cvStuff[1,1] = (cfac*nopt)/nli
        cvStuff[2,1] = wn
        cvStuff[3,1] = yv
        cvStuff[4,1] = yb
        for k in range(nli, nk+1): 
            i, j = array(k+4,n)
            cvStuff[i,j] = b/sw
            k1 = k*(nk+1)+3
            k1, cvStuff = copyAtoSC(n, nk, k1, nBasisFunctions, a, A, cvStuff)
        i, j = array((nk+1)**2+4, n)
        cvStuff[i,j] = nBasisFunctions
        kl = nli
        
    #backward stepwise elimination - algorithm 3
    for ll in range(2, nli+1):
        b, A, a, k = \
            bkstp(D, nBasisFunctions, XB, yb, parameters['marsgo']['alr'], DP)
        if k != 0:
            if parameters['cvmars']['ix'] != 0:
                i, j = array(kl+3, n)
                cvStuff[i,j] = b/sw
                kl = kl-1
                k1 = kl*(nk+1) + 3
                k1, cvStuff = \
                    copyAtoSC(n, nk, k1, nBasisFunctions, a, A, cvStuff)
            tcst = tcst-df1
            b = GCV(b/sw, tcst, wn)
            if b < asm:
                asm = b
                tcsts = tcst
                az = a
                tb = copyAtoTB(A, tb)

    if txm <= asm:
        asm = txm
        tcsts = 1.0
        az = yb
        tb[1,1:] = 0.0

    if LOG:
        logger.info('\n  final model after backward stepwise elimination:')
        coefpr(nBasisFunctions, az, tb, cm, xs)
        logger.info('   (piecewise linear) gcv = %12.4g   #efprms = %5.1f'
                    %(asm,tcsts))
    return az, tb, cm, cvStuff
def buildNextBasisFunction(parent0, nextVar,
                           n, p, x, y, yb, w, sw, nk, lx, tb, cm, mm, mi, 
                           db, bfIndex, bfIndexLast, evaluations, 
                           ms, fv, newBF, tx1, txi, 
                           nBasisFunctions, nep, ssq, rsq, asq0,
                           kcp0, nop, jas,
                           kr, nc, kcp, k1, DY, fln, eps, big, nmin, alf):
    """ 
    Basis functions are products of hockey stick functions where each variable
    is checked and the knot is selected from the data in that coordinate. There
    is no repetition of a variable and selection is subject to the constraint  
    of max interaction(mi), an input parameter.  
    """

    parent = parent0
    while parent >= 0:
        txl = big
        NNORD = nnord(pymars.NEST0, parent, tb)
        if NNORD >= mi:
            pymars.ADDPAR0.updpar(0,-1.0)
        else:
            nnt, partialEval = blf0(pymars.NEST0, parent, 0, n, x, w, cm, evaluations[:,parent])
            lbf = 0
            if nnt <= nmin:
                pymars.ADDPAR0.updpar(0,-1.0)
            else:
                nep = 0
                for variable in range(1, p+1):
                    #check if first point is smaller than last point
                    if x[mm[1,variable],variable] < x[mm[n,variable],variable]:
                        #check if variable is already included
                        if jf(parent, variable, tb) == 0:
                            ja = pymars.NEST0.isfac(parent, variable, bfIndexLast, tb, cm)
                            if ja >= 0:
                                if elg(pymars.NEST0, pymars.SETINT0, 
                                       variable, parent, lx, tb, cm):
                                    nep = nep+1               
                if nep == 0:
                    pymars.ADDPAR0.updpar(0,-1.0)
                else:
                    mn, me, mel = mnspan(ms, alf, nep, nnt)
                    if nnt <= max(me,mel):
                        pymars.ADDPAR0.updpar(0,-1.0)
                    else:
                        (nextVar, newHS, cm, db, fv, tx1, txi, ssq, rsq, 
                         kcp0, kcp, nop, ja, jas, lbf, nnt, nc, kr, k1, DY) = \
                        findBestChildVariable(parent, nextVar, 
                                          n, p, x, y, yb, w, sw, lx, tb, cm, mm, db,
                                          bfIndex, bfIndexLast, evaluations, partialEval,
                                          ms, fv, newBF, txl, tx1, txi, 
                                          nBasisFunctions, nep, ssq, rsq, kcp0, kcp,
                                          mn, me, mel, nop, 
                                          ja, jas, lbf, nnt, kr, k1, nc, DY,
                                          fln, eps, big, nmin, alf)
                        pymars.ADDPAR0.updpar(nextVar, asq0-tx1)
        parent, nextVar = pymars.ADDPAR0.nxtpar(parent, nextVar)
    return tb, cm, mm, evaluations, db, bfIndex, bfIndexLast, \
           ms, fv, newHS, txl, tx1, txi, ssq, rsq, kcp0, \
           nop, nextVar, jas, kr, \
           nc, kcp, k1, DY
def copyAtoSC(n, nk, k1, nBasisFunctions, a, A, cvStuff):
    l = 0
    while l <= nk:
        k1 = k1+1
        i, j = array(k1, n)
        if l == 0:
            cvStuff[i,j] = a
        elif (l > nBasisFunctions):
            cvStuff[i,j] = 0.0
        else:
            cvStuff[i,j] = A[l]#d[l,2]
        l = l+1
    return k1, cvStuff
def copyAtoTB(A, tb):
    tb[1,1:] = 0.0
    (i,) = numpy.where(A != 0.0)
    tb[1,i] = A[i]
    return tb        
def addBasisFunction(n, x, y, w, sw, yb, tb, cm,
                     evaluations, partialEval, k1, kr, rsq, bfIndex, newBF, db, DY):
    tb[1:, bfIndex] = newBF[1:]
    k1 = kr
    J = int(abs(tb[2,bfIndex])+.1)
    kr, currentEval, db, DY = update(2, n, kr, x[:,J], y, w, sw, yb, 
                                   tb[2,bfIndex], tb[3,bfIndex], 
                                   cm, cm[2*J], partialEval, db, DY)
    evaluations[:,bfIndex] = currentEval
            
    if kr > k1:
        rsq = rsq - DY[kr]**2
    pymars.ADDPAR0.addpar(bfIndex)      
    return tb, evaluations, k1, kr, rsq, db, DY
def addReflectedBasisFunction(n, x, y, w, sw, yb, tb, cm, variable, 
                              evaluations, partialEval, k1, kr, rsq, bfIndex, newBF, db, DY):
    """add the reflected basis function; only the orientation changes"""
    tb[1:5, bfIndex] = newBF[1:5]
    tb[2, bfIndex] = -tb[2, bfIndex]
    if cm[2*variable] > 0.0: 
        for i in range(1,n+1):
            evaluations[i,bfIndex] = phi(bfIndex, x[i,:], tb, cm)
    else:    
        k1 = kr
        J = int(abs(tb[2,bfIndex])+.1)
        kr, currentEval, db, DY = update(2, n, kr, x[:,J], y, w, sw, yb, 
                                       tb[2,bfIndex], tb[3,bfIndex], 
                                       cm, cm[2*J], partialEval, db, DY)

        evaluations[:,bfIndex] = currentEval
        if kr > k1:
            rsq = rsq - DY[kr]**2
    pymars.ADDPAR0.addpar(bfIndex)
    return tb, evaluations, k1, kr, rsq, db, DY
def printBasisFunction(xm, xs, bfIndex, mtot, kr, tcst, cm, newBF, rsq, sw, wn):
    mp = bfIndex-1
    variable = int(abs(newBF[2])+.1)
    fkr = kr
    gcv = GCV(rsq/sw, tcst, wn)
    if cm[2*variable] > 0.0:
        hol = holl(variable, cm, newBF[3])
        if bfIndex == mtot+1:
            logger.info('  %3i    %12.4g   %5.1f   %5.1f       %3.0d %28s %3.0d'
                        %(bfIndex, gcv, fkr, tcst, newBF[2], hol, newBF[4]))                         
        if bfIndex == mtot+2:
            logger.info('%3i  %3i %12.4g   %5.1f   %5.1f       %3.0d %28s %3.0d'
                        %(bfIndex, mp, gcv, fkr, tcst, newBF[2], hol, newBF[4]))
    else: 
        xk = xm[variable] + xs[variable]*newBF[3]
        if bfIndex == mtot+1: 
            logger.info('  %3i   %12.4g   %5.1f    %5.1f     %3.0f   %12.4g      %8.0d'
                        %(bfIndex, gcv, fkr, tcst, newBF[2], xk, newBF[4]))
        if bfIndex == mtot+2: 
            logger.info('%3i  %3i %12.4g   %5.1f    %5.1f    %3.0d   %12.4g      %8.0d'
                        %(bfIndex, mp, gcv, fkr, tcst, newBF[2], xk, newBF[4]))
    return
def processNestedData(n, x, y, w, sw, yb, tb, cm, jas, jn, kcp0, kcp, kr, rsq, 
                      nBasisFunctions, newBF, bfIndex, evaluations, partialEval, db, DY):
    VALS = cm[kcp0+1:kcp+1]
    VALS = border(VALS)
    jn, VALS = NEST0.getnst(jas, cm, jn, kcp, VALS)#needs work
    cm[kcp0+1:kcp+1] = VALS[1:]
    tb[2, bfIndex] = jn
    tb[3, bfIndex] = kcp0
    kcp0 = kcp0+kcp
    tb[4, bfIndex] = newBF[4]
    k1 = kr
    INDEX = int(newBF[4]+.1)
    partialEval = blf(INDEX, n, evaluations[:,INDEX])#sc)
    #sc[1:n+1,mkp1] = partialEval[1:]
    newBF[4] = bfIndex
            
    #DY = d[0:, 1]
    #DB = d[0:, 3]
    J = int(abs(tb[2, bfIndex])+.1)
    kr, currentEval, db, DY = update(2, n, kr, x[:,J], y, w, sw, yb, 
                                    tb[2,bfIndex], tb[3,bfIndex], cm, cm[2*J],
                                    partialEval, db, DY)
    #d[0:,1] = DY[0:]
    #d[0:,3] = DB[0:]
    evaluations[:,bfIndex] = currentEval

    if kr > k1:
        rsq = rsq - DY[kr]**2
    pymars.ADDPAR0.addpar(bfIndex)
    if bfIndex < nBasisFunctions: 
        bfIndex = bfIndex + 1
        tb[2, bfIndex] = -tb[2, bfIndex-1]
        for i in [3,4]:
            tb[i, bfIndex] = tb[i, bfIndex-1]
        if ibfext(bfIndex, tb, cm) > 0:
            bfIndex = bfIndex - 1
        else:    
            for i in range(1, n+1):
                evaluations[i, bfIndex] = phi(bfIndex, x[i,:], tb, cm)
            pymars.ADDPAR0.addpar(bfIndex)
    if LOG: 
        mp = bfIndex - 1
        tcst = (nopt-1)*df1 + kr + 1.0
        fjn = jn
        fkr = kr
        gcv = GCV(rsq/sw, tcst, wn) 
        hol = holl(jn, cm, tb[3,bfIndex])
        if bfIndex == mtot+1:
            logger.info('   %3i     %12.4g   %5.1f   %5.1f       %3.0d %28a %3.0d'
                        %(bfIndex, gcv, fkr, tcst, fjn, hol, tb[4, bfIndex]))
        if bfIndex == mtot+2:
            logger.info(' %3i   %3i    %12.4g   %5.1f   %5.1f       %3.0d %28a %3.0d'
                        %(bfIndex, mp, gcv, fkr, tcst, fjn, hol, tb[4, bfIndex]))
    mtot = bfIndex
    bfIndex = bfIndex + 1
    return tb, cm, kcp0, k1, kr, newBF, bfIndex, evaluations, \
           jn, rsq, db, DY, tcst, fjn, fkr, gcv
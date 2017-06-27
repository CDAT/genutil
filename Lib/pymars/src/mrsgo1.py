import genutil.pymars, time, sys, scipy, numpy, operator, pdb
from .mrsgo2_no_ctypes import findBestKnot
#from mrsgo2Cython import findBestKnot
from .elg import elg
from .mnspan import mnspan
from .jft import jft
from .jf import jf
from .csp import csp
from .finalDecision import *
from .border import *
from genutil.pymars import LOG, logger, debug, parameters, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE, MPI
from .newb import isNewHS
def findBestChildVariable(parent, nextVar, 
           n, p, x, y, yb, w, sw, lx, tb, cm, mm, db,
           bfIndex, bfIndexLast, parentEval, partialEval,
           ms, fv, newHS, txl, tx1, txi, 
           nBasisFunctions, nep,  rsq, kcp0, kcp,
           mn, me, mel,  
           ja, lbf, nnt, kr,  DY,
           fln, eps, big, nmin, alf):
    """
    Find the best child variable for the partially constructed response surface.
    """
    try:
        import mpi
        PARALLEL = True
        rank = mpi.rank
        nprocs = mpi.size
    except :
        PARALLEL = False

    START = time.time()
    if nextVar == 0:
        VARIABLES = range(1,p+1)
    else:
        VARIABLES = [nextVar]

    catKnots = categoricalKnots(p, VARIABLES, lx, cm, kcp, kcp0)

    commonInput = (parent, n, p, y, yb, w, sw, lx, tb, cm, db,
                 bfIndex, bfIndexLast, parentEval, partialEval,
                 ms, fv, nBasisFunctions, nep,  rsq, kcp0, catKnots,
                 mn, me, mel, ja, lbf,  nnt, kr,   DY,
                 fln, eps, big, nmin, alf)
    if PARALLEL:
        sent = []
        if len(VARIABLES) >= 2:
            for variable in VARIABLES[1:]: 
                proc = variable-1
                sent += [proc]
                X  = x[:,variable]  #data in one coordinate
                MM = mm[:,variable] #memory map in the same coordinate       
                input = (variable, X, MM) + commonInput           
                mpi.send(input, proc)    
            VARIABLES = [VARIABLES[0]]

    states = []
    for variable in VARIABLES: 
        X  = x[:,variable]  #data in one coordinate
        MM = mm[:,variable] #memory map in the same coordinate       
        input = (variable, X, MM) + commonInput
        
        (CM, nc, ict,   RSQ,  ja, fvr,  prelimDecision, defaultHS, proposedHS) = \
                computeProposedHockeyStick(input)

        hockeyStickState = makePreliminaryDecision(prelimDecision, tb,  nc, ict, bfIndex, fln, fvr, sw, 
                                                   RSQ,  defaultHS, proposedHS, CM, variable, ja)
        
        states += [hockeyStickState]

    if PARALLEL:
        received = []
        while received != sent:
            hockeyStickState, status = mpi.recv()
            states += [hockeyStickState]
            received += [status.source]
            received.sort()
                
    nextVar, newHS, tx1, txi, jas, cm = decideWinner(states, txl, tx1, txi, cm, kcp0, nextVar, newHS)
    END = time.time()
    if 'marsgo1' in genutil.pymars.TIME.keys():
        t0 = END-START
        genutil.pymars.TIME['marsgo1'][0] += t0
        genutil.pymars.TIME['marsgo1'][2] += 1
        if len(VARIABLES) > 1:
            genutil.pymars.TIME['marsgo1'][1] += t0
            genutil.pymars.TIME['marsgo1'][3] += 1
    #debug.info('END findBestVariable:newHS='+repr((newHS[1:5])))
    #debug.info(' ')
    return nextVar, newHS, cm, tx1, txi, jas
def computeProposedHockeyStick((variable, X, MM, parent, 
                 n, p, y, yb, w, sw, lx, tb, cm, db,
                 bfIndex, bfIndexLast, parentEval, partialEval,
                 ms, fv, nBasisFunctions, nep,  rsq, kcp0, catKnots,
                 mn, me, mel, ja, lbf,  nnt, kr, DY,
                 fln, eps, big, nmin, alf)):
    """
    Compute a hockey stick in this variable with the best knot location.
    """
    CM = None
    ict = 0
    fvr = 1.0
    prelimDecision = False
    nc = -1
    proposedHS = 4*[0]
    defaultHS = 4*[0]
    RSQ = rsq
    
    if X[MM[1]] < X[MM[n]]:
        #check if variable is already included
        if jf(parent, variable, tb) == 0: 
            ja = genutil.pymars.NEST0.isfac(parent, variable, bfIndexLast, tb, cm)
            if ja >= 0: 
                if elg(genutil.pymars.NEST0, genutil.pymars.SETINT0, variable, parent, lx, tb, cm):
                    if ja == 0:
                        if lbf != 0:
                            nnt, partialEval = blf0(genutil.pymars.NEST0, parent, 0, n, X, w, cm, parentEval)
                            lbf = 0
                            mn, me, mel = mnspan(ms, alf, nep, nnt)
                    else: 
                        nnt, partialEval = blf0(genutil.pymars.NEST0, parent, ja, n, X, w, cm, parentEval)
                        lbf = 1
                        if nnt > nmin: 
                            mn, me, mel = mnspan(ms, alf, nep, nnt)
                    if ja == 0 or (ja != 0 and nnt > nmin and nnt > max(me, mel)):
                        fvr = 1.0
                        #make it easer to use this variable if it is not already used.
                        if variable not in abs(tb[2,1:bfIndexLast+1]):
                            fvr = 1.0 + fv
                            
                        ict = 0                        
                        if lx[variable] < 0:
                            #categorical variable
                            ict = 1
                            nc = int(cm[2*variable+1]+.1) - int(cm[2*variable]+.1) + 1
                            coef, knot, CM, nop = \
                                csp(nc, n, X, y, w, catKnots, variable,
                                    yb, kr, nnt, sw, me, partialEval, db, DY)

                            prelimDecision = (nop != 0)
                            defaultHS = [coef, variable, knot, parent]
                            proposedHS = defaultHS
                        else:
                            #continuous variable

                            (db,  rsq,   prelimDecision, DY, coefAndKnot) = \
                                findBestKnot(n, X, y, w, lx,  cm, MM,  partialEval, db, 
                                             bfIndex, nBasisFunctions, yb, sw, rsq, 
                                             eps, big, parent, mn, me, mel,  
                                             variable, nnt, kr, DY)
                            coef, knot = coefAndKnot['default']
                            defaultHS = [coef, variable, knot, parent]
                            proposedHS = defaultHS
                            if 'proposed' in coefAndKnot.keys():
                                coef, knot = coefAndKnot['proposed']
                                proposedHS = [coef, variable, knot, parent]
                                #debug.info('after findBestKnot proposedHS=' + repr(proposedHS))

    return CM, nc, ict,  rsq, ja, fvr, prelimDecision, defaultHS, proposedHS
def makePreliminaryDecision(prelimDecision, tb,  nc, ict, bfIndex, fln, fvr, sw, 
                            RSQ, defaultHS, proposedHS, CM, variable, ja):

    """
    This function makes some changes to the proposed hocky stick before 
    making a preliminary acceptance decision.
    COMPARE is either <= for continuous or < for categorical variable.
    """

    makeDecision = False
    HS = None
    COMPARE = None

    if prelimDecision:
        proposedHS[0] = (RSQ + proposedHS[0])/sw #coefficient
    
        if ict == 0: 
            if defaultHS[0] <= fln*proposedHS[0]:
                index = bfIndex
                HS = defaultHS
            else:
                index = bfIndex + 1
                HS = proposedHS                
        else:
            index = bfIndex
            HS = proposedHS
        #In the case of categorical variables makeDecision will always be true.
        #The reason is that the knot is kcp0+nc and kcp0 is incremented once
        #a hockey stick for a categorical variable is found.
        makeDecision = isNewHS(HS[1:], tb[2:5, 1:index])
        COMPARE = operator.lt
    else:
        if ict == 0:
            HS =defaultHS
            HS[0] = RSQ/sw
            makeDecision = isNewHS(HS[1:], tb[2:5, 1:bfIndex])
            COMPARE = operator.le

    hockeyStickState = hsState(makeDecision, variable, HS, COMPARE, fvr, ja, ict, nc, CM)
    return hockeyStickState
def categoricalKnots(p, VARIABLES, lx, cm, kcp, kcp0):
    """
    This function computes the potential knots locations for categorical variables.
    """
    catKnots = numpy.zeros(p+1)
    firstCat = True
    for variable in VARIABLES:
        if lx[variable] < 0:
            if firstCat:
                nc = int(cm[2*variable+1]+.1) - int(cm[2*variable]+.1) + 1 
                knot = kcp
                ncLast = nc
                firstCat = False
            else:
                nc = int(cm[2*variable+1]+.1) - int(cm[2*variable]+.1) + 1  
                knot = kcp0+ncLast
                ncLast = nc
        else:
            knot = None

        catKnots[variable] = knot

    return catKnots
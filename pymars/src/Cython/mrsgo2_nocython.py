#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy, time, pymars
from pymars.update import update
from pymars.newb import isNewHS
from pymars.decideNextVariable import decideNextVariable
from pymars import debug, ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
if pymars.MPI:
    import mpi

def findBestKnot(n, X, y, w, lx, cm, MM,  partialEval, db, 
           bfIndex, nBasisFunctions, yb, sw, rsq, 
           eps, big, parent, mn, me, mel, variable, nnt, kr,  DY):
    """ 
    This function computes the best knot & coefficient in a specified variable.
    """
    #START = time.time()
    
    knot = X[MM[1]]
    
    k1 = kr
    kr, currentEval, db, DY = update(1, n, kr, X, y, w, sw, yb, variable, knot,
                                     cm, cm[2*variable], partialEval, db, DY)
    if kr > k1:
        rsq = rsq - DY[kr]**2
        coef = rsq/sw
    else:
        coef = big
    
    coefAndKnot = {}
    coefAndKnot['default'] = (coef, knot)
    
    KEEP = False
    EXECUTE = \
        not (lx[variable] == 3 or bfIndex >= nBasisFunctions or nnt <= me+mel)
    if EXECUTE:
        START = time.time()

        #find the best knot & coefficient
        coef, knot = \
            bestKnot(n, X, y, yb, w, sw, MM, partialEval, db,  
                     eps, big, me, mel, mn, nnt, kr, DY)#fvr, DY)
        KEEP = True
        coefAndKnot['proposed'] = (coef, knot)
        
        END = time.time()
        if 'marsgo2' in pymars.TIME.keys(): pymars.TIME['marsgo2'] += [(START, END)]

    return db, rsq, KEEP, DY, coefAndKnot 
def getPotentialKnots(n, X, MM, w, partialEval, me, mel, mn, nnt):
    """This function gets a list of pairs of indices starting at the
    top of the range and working backwards to the bottom.  It uses 
    the skipping algorithm described on page 24 of the Friedman paper. 
    Note at the very beginning the input data are sorted for this
    reason. Also this calculation could be performed right after
    the sort and stored for later use. """
    nnl = nnt
    nst = 0
    nnr = -1
    j = n
    potentialKnots = []
    while j > 1:
        j0 = j
        while j > 1:
            mj = MM[j]
            h = partialEval[mj]
            if w[mj] > 0.0 and h > 0.0:
                nst = nst+1
                nnl = nnl-1
                nnr = nnr+1
            if (X[MM[j-1]] >= X[MM[j]] or
                nst < mn or nnl < mel or nnr < me):
                j = j-1
            else:
                nst = 0
                break
        potentialKnots += [(j, j0)]
        j = j-1
    return potentialKnots
def bestKnot(n, X, y, yb, w, sw, MM, partialEval, db, 
             eps, big, me, mel, mn, nnt, kr, DY):
    """"This function is the controller that finds the best hockey
    stick in the specified variable. It relies on the partial response
    surface, that is the portion of the response surface built to this
    point in the calculation.  This is the array partialEval."""
    inputs = []
    trialKnot = 0.0
    coef = big
    knot = 0.0

    v, u, t, we, se, sy, dy = 7*[0.0]
    D = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    DB = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    txt = X[MM[1]] + X[MM[n]]
    xt = 0.5*txt #midpoint between the smallest & largest
    PARALLEL = False#pymars.MPI
    #define the maximum number of inputs
    #the limit in the following depends on how much memory is available
    #to store the inputs.
    MAXINPUTS = 2000
    potentialKnots = getPotentialKnots(n, X, MM, w, partialEval, me, mel, mn, nnt)  

    #this seems to be equation 52 in the paper
    for I,(j,j0) in enumerate(potentialKnots):
        FIRST_TRIAL = (j0 == n)
        LastPair = (I+1 == len(potentialKnots))
        if j > 1:
            lastTrialKnot = trialKnot
            trialKnot = X[MM[j]]
            input = \
                (FIRST_TRIAL, lastTrialKnot, 
                 trialKnot, X, y, yb, w, sw, MM[j:j0+1], xt, partialEval, kr, db)
            inputs += [input]
        COMPUTE_COEF = False    
        if len(inputs) >= MAXINPUTS or LastPair:
            outputs = []
            if PARALLEL:
                inputs = mpi.scatter(inputs)

            #else:
            for input in inputs:
                FIRST_TRIAL, lastTrialKnot, trialKnot = input[0:3]
                V, T, U, WE, SE, SY, DYY, ST, SU = computeCoef(input[2:])
                outputs += \
                        [(FIRST_TRIAL, lastTrialKnot, trialKnot, V, T, U, WE, SE, SY, DYY, ST, SU)]
            inputs = []
            if PARALLEL:
                try:
                    outputs = mpi.gather(outputs)
                except: 
                    mpi.abort()
            COMPUTE_COEF = True
        if COMPUTE_COEF:
            for (FIRST_TRIAL, lastTrialKnot, trialKnot, V, T, U, WE, SE, SY, DYY, ST, SU) in outputs:

                if not FIRST_TRIAL:#j0 != n:
                    dx = lastTrialKnot - trialKnot
                    dy = dy + dx*sy
                    we = we + dx*se
                    v = v + dx*(2.0*u - (lastTrialKnot + trialKnot - txt)*t)
                    D[1:kr+1] = D[1:kr+1] + dx*DB[1:kr+1]
                #accumulate the local contributions to the best knot & coefficient
                v += V
                t += T
                u += U
                we += WE
                se += SE
                sy += SY
                dy += DYY
                D[1:kr+1]  = D[1:kr+1]  + ST[1:kr+1]
                DB[1:kr+1] = DB[1:kr+1] + SU[1:kr+1]
                #perform the check that a new best knot & coefficient has been found
                newCoef, Coef = checkCoef(coef, eps, kr, sw, v, we, se, dy, D, DY)
                if newCoef:
                    coef = Coef
                    knot = trialKnot
 
    return coef, knot
def computeCoef((trialKnot, X, y, yb, w, sw, MM, xt, partialEval, kr, db)): 
                #v, u, t, we, se, sy, txt, xt, dy, partialEval, kr, db)):  
    """This function computes the local calculations involved in deciding
    a new knot and coefficient. It is called from bestKnot and once complete
    the accumulation of the local contributions is done before it is decided
    if the best one is found."""
    #print 'computeCoef'
    v, u, t, we, se, sy, dy = 7*[0.0]
    ST = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    SU = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    for mj in MM:
        h = partialEval[mj]
        if w[mj] > 0.0 and h > 0.0:
            xx = X[mj]
            xd = xx - trialKnot
            su = w[mj]*h
            st = su*xd
            yc = y[mj] - yb
            dy = dy + st*yc
            sy = sy + su*yc
            we = we + st
            se = se + su
            sj = w[mj]*h**2
            v = v + sj*xd**2
            t = t + sj
            u = u + sj*(xx - xt)
            ST[1:kr+1] = ST[1:kr+1] + st*db[mj,1:kr+1]
            SU[1:kr+1] = SU[1:kr+1] + su*db[mj,1:kr+1]
    return v, t, u, we, se, sy, dy, ST, SU
def checkCoef(coef, eps, kr, sw, v, we, se, dy, D, DY):
    """This function checks to see if the proposed coefficient is best."""
    dv = v-we**2/sw
    if dv > 0.0:
        a = numpy.inner(D[1:kr+1], DY[1:kr+1])
        b = numpy.inner(D[1:kr+1], D[1:kr+1])
        b = dv-b
        if b > eps*dv:
            b = -(dy-a)**2/b
            #this check is a very important one in the algorithm. if it
            #passes the new hockey stick might be accepted and the response
            #surface changes. In the original code there was no accounting for 
            #numerical issues. The addition of 1.e-10 was added only after
            #some painful debugging. The effect of adding this is to say 
            #that the new coefficient (b) is taken as the new one if it's
            #really smaller.
            if b + 1.e-10< coef:
                return True, b #new coefficient
    return False, None
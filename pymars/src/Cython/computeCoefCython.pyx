import numpy as np
cimport numpy as np

#def computeCoef((trialKnot, X, y, yb, w, sw, MM, xt, partialEval, kr, db)): 
                #v, u, t, we, se, sy, txt, xt, dy, partialEval, kr, db)):  
def computeCoef(double trialKnot, 
                np.ndarray[np.float64_t, ndim=1] X, 
                np.ndarray[np.float64_t, ndim=1] y, 
                double yb, 
                np.ndarray[np.float64_t, ndim=1] w, 
                double sw, 
                np.ndarray[np.int_t, ndim=1] MM, 
                double xt,  
                np.ndarray[np.float64_t, ndim=1] partialEval, 
                int kr, 
                np.ndarray[np.float64_t, ndim=2] db):

    """This function computes the local calculations involved in deciding
    a new knot and coefficient. It is called from bestKnot and once complete
    the accumulation of the local contributions is done before it is decided
    if the best one is found."""
    
    #returned parameters 
    cdef double v, u, t, we, se, sy, dy
    cdef np.ndarray[np.float64_t, ndim=1] ST, SU 
    
    #local parameters
    cdef int mj
    cdef double h, xx, xd, su, st, yc, sj
    
    #print 'computeCoef'
    #output = ComputeCoef.output(input.FIRST_TRIAL, input.lastTrialKnot, input.trialKnot, input.kr)
    v, u, t, we, se, sy, dy = 7*[0.0]
    ST = np.zeros(shape = kr+1, dtype = np.float)
    SU = np.zeros(shape = kr+1, dtype = np.float)
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
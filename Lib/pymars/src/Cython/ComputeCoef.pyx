#cimport numpy
import numpy as np
cimport numpy as np
#import numpy

import ComputeCoef
cimport ComputeCoef

cdef class input:
    cdef public long FIRST_TRIAL
    cdef public double lastTrialKnot 
    cdef public double trialKnot
    cdef public np.ndarray X
    cdef public np.ndarray y
    cdef public double yb
    cdef public np.ndarray w
    cdef public double sw 
    cdef public np.ndarray MM
    cdef public double xt 
    cdef public np.ndarray partialEval
    cdef public int kr
    cdef public np.ndarray db
    
    def  __init__(self, long FIRST_TRIAL, 
                        double lastTrialKnot, 
                        double trialKnot, 
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
        self.FIRST_TRIAL = FIRST_TRIAL
        self.lastTrialKnot = lastTrialKnot
        self.trialKnot = trialKnot
        self.X = X
        self.y = y
        self.yb = yb
        self.w = w
        self.sw  = sw
        self.MM = MM#[j:j0+1]
        self.xt = xt
        self.partialEval = partialEval
        self.kr = kr
        self.db = db
    
        return
cdef class output:
    cdef public int FIRST_TRIAL
    cdef public double lastTrialKnot 
    cdef public double trialKnot
    cdef public double v
    cdef public double t
    cdef public double u
    cdef public double we
    cdef public double se 
    cdef public double sy
    cdef public double dy 
    cdef public np.ndarray ST
    cdef public np.ndarray SU
    
    def  __init__(self, int FIRST_TRIAL, 
                        double lastTrialKnot, 
                        double trialKnot,
                        int kr):
        self.FIRST_TRIAL = FIRST_TRIAL
        self.lastTrialKnot = lastTrialKnot
        self.trialKnot = trialKnot
        self.v = 0.0
        self.t = 0.0
        self.u = 0.0
        self.we = 0.0
        self.se = 0.0
        self.sy = 0.0
        self.dy = 0.0
        self.ST = np.zeros(shape = kr+1, dtype=np.float)
        self.SU = np.zeros(shape = kr+1, dtype=np.float)
    
        return
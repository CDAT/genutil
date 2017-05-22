import numpy as np
import ctypes
from pymars import cmars, debug

def computeCoef(trialKnot, X, y, yb, w, sw, MM, xt, partialEval, kr, db, db_nrows, db_ncols): 
    class output(ctypes.Structure):
        _fields_ = [('v', ctypes.c_double),
                    ('u', ctypes.c_double),
                    ('t', ctypes.c_double),
                    ('we', ctypes.c_double),
                    ('se', ctypes.c_double),
                    ('sy', ctypes.c_double),
                    ('dy', ctypes.c_double),
                    ('ST', ctypes.POINTER(ctypes.c_double)),
                    ('SU', ctypes.POINTER(ctypes.c_double))]

 
    cmars.computeCoef.argtypes = [ctypes.c_double,
                                  np.ctypeslib.ndpointer(np.double),
                                  np.ctypeslib.ndpointer(np.double),
                                  ctypes.c_double,
                                  np.ctypeslib.ndpointer(np.double),
                                  ctypes.c_double,
                                  np.ctypeslib.ndpointer(np.int64),
                                  ctypes.c_int,
                                  ctypes.c_double,
                                  np.ctypeslib.ndpointer(np.double),
                                  ctypes.c_int,
                                  np.ctypeslib.ndpointer(np.double),
                                  ctypes.c_int,
                                  ctypes.c_int]
    
    cmars.computeCoef.restype = output
    OUT = cmars.computeCoef(trialKnot, X, y, yb, w, sw, MM, len(MM), xt, partialEval, kr, db, db_nrows, db_ncols)
    
    ST = np.ctypeslib.as_array(OUT.ST, shape=(kr+1,))
    SU = np.ctypeslib.as_array(OUT.SU, shape=(kr+1,))
    
    return OUT.v, OUT.t, OUT.u, OUT.we, OUT.se, OUT.sy, OUT.dy, ST, SU
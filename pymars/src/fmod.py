#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.fmrs import fmrs
from pymars.cmrs import cmrs

def fmod (m, n, nk, x, az, tb, cm, kp, kv, lp, lv, bz, tc):
    if m == 1:
        f = fmrs(n, x, nk, az, tb, cm)
    else:
        f = cmrs(n, x, cm, kp, kv, lp, lv, bz, tc)
    return f 
from .fmrs import fmrs
from .cmrs import cmrs

def fmod (m, n, nk, x, az, tb, cm, kp, kv, lp, lv, bz, tc):
    if m == 1:
        f = fmrs(n, x, nk, az, tb, cm)
    else:
        f = cmrs(n, x, cm, kp, kv, lp, lv, bz, tc)
    return f 
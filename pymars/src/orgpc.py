#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.scpc import scpc
from pymars.border import border
import numpy
def orgpc(xm, xs, lp, lv, tc):

    la = 1
    l1 = la

    while lp[1, l1] != 0:
        l = lp[1,l1]
        nt = lp[3,l1]
        lb = la + 5*l*nt - 1
        for j in range(1, nt+1):
            jv = numpy.array(lv[lp[2,l1]:lp[2,l1]+l])
            jv = border(jv) 
            TC = numpy.array(tc[la:la+5*l*nt])
            TC.shape = 5*l,nt
            TC = TC.T
            TC = border(TC)
            TC, b = scpc(xm, xs, j, l, nt, jv, TC, tc[lb+j])
            TC = TC[1:,1:]
            TC = TC.flatten()
            tc[la:la+5*l*nt] = TC
            tc[lb+j] = b
        la = lb + nt + 1
        l1 = l1 + 1
    return tc

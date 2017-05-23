#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from .que import que
from .border import border
from .icat import icat
from .basisFunction import categoricalBasisFunction
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
def cmrs(n, x, cm, kp, kv, lp, lv, bz, tc):
    #print 'cmrs'
    y = numpy.array([0.0] + n*[bz], dtype=FLOAT_DTYPE)
    sc = numpy.zeros(shape = (n+1,3), dtype=FLOAT_DTYPE)

    ll = 1
    la = ll
    l1 = la
    while kp[1, ll] >= 0:
        sc[1:n+1,1] = 1.0
    
        FLAG1 = (kp[1,ll] > 0)
        FLAG2 = (kp[3,ll] < 0)
        FLAG3 = (kp[3,ll] > 0)

        if FLAG1:
            jl = int(kp[1,ll])
            for il in range(1, jl+1):
                k = int(kp[2,ll]) + il - 1
                jj = int(kv[1,k])
                j = abs(jj)
                kk = int(kv[2,k])
                for i in range(1, n+1):
                    if sc[i,1] != 0.0:
                        if parameters['cmrs']['ifg'] == 0:
                            ic = icat(x[i,j], j, cm)
                        else: 
                            ic = int(x[i,j]+.1) 
                        if ic == 0:
                            sc[i,1] = 0.0
                        else:   
                            sc[i,1] = cm[ic+kk]
                        if jj < 0:
                            if sc[i,1] == 0.0:
                                sc[i,1] = 1.0
                            else:  
                                sc[i,1] = 0.0
    
        if FLAG1 and FLAG2:
            print 'cmrs'
            k = -kp[3,ll]
            (i,) = numpy.where(sc[i,1] !=0.0)
            i = i+1
            y[i] = y[i] + tc[k]
            #for i in range(1, n+1):
                #if sc[i,1] != 0.0:
                    #y[i] = y[i] + tc[k]
    
        if (FLAG1 and not FLAG2) or (not FLAG1 and FLAG3):
            kp3 = kp[3,ll]
            for m in range(1, kp3+1):
                l = lp[1,l1]
                nt = lp[3,l1]
                lb = la + 5*l*nt-1
                for j in range(1, nt+1):
                    sc[1:n+1,2] = sc[1:n+1,1]
                    pnt = int(lp[2,l1])
                    LV = lv[pnt: pnt+l]
                    LV = border(LV)
                    TC = numpy.array(tc[la: la+nt*5*l])
                    TC.shape = 5*l, nt
                    TC = TC.T
                    TC = border(TC)
                    sc[:,2] = que(j, l, nt, LV, n, x, TC, sc[:,2])
                    y[1:n+1] = y[1:n+1] + tc[lb+j]*sc[1:n+1,2]
                la = lb+nt+1
                l1 = l1+1
    
        ll = ll+1

    return y

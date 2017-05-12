#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.basisFunction import defaultBasisFunction,categoricalBasisFunction
def cvmod (i, n, x, y, w, nk, mk, tb, cm, sc, cv0, cv):
    #print 'cvmod'
    #very similar to phi
    for m in range(1, nk+1): 
        orient = tb[2,m]
        knot = tb[3,m]
        j = int(abs(orient)+.1)
        if cm[2*j] > 0.0: 
            u = categoricalBasisFunction(orient, knot, cm, int(x[i,j]+.1), x[i,j])
        else:
            u = defaultBasisFunction(tb[2,m], tb[3,m], 1.0, x[i,j])
        l = int(tb[4,m]+.1)
        if l > 0: 
            cv[m,2] = u*cv[l,2]
        else:
            cv[m,2] = u
            
    kp = nk+4
    cv0 = cv0 + w[i]*(y[i]-sc[4])**2
    for m in range(1, nk+1):
        kp = kp+1
        s = sc[kp]
        for l in range(1, nk+1):
            kp = kp+1
            if l <= mk:
                s = s + sc[kp]*cv[l,2]
        cv[m,1] = cv[m,1] + w[i]*(y[i]-s)**2
    return cv0, cv

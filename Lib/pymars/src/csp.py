#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE, debug
def csp(nc, n, X, y, w,  catKnots, variable,
        yb, kr, ntt, sw, me,  partialEval, db, DY):
    """ It appears to have similar functionality as marsgo2
        only for categorical variables. It is called in marsgo1.
    """
    #print "csp"
    eps, big = (1.e-4, 9.9e30)
    nop = 0
    if nc <= 1: 
        coef = big
        defaultKnot = X[1]
        CM = None
        return coef, defaultKnot, CM, nop
    n1 = nc + 1

    sp = numpy.zeros(shape = (kr+1, nc+2), dtype=FLOAT_DTYPE)
    sp1 = numpy.zeros(shape = (3+1, nc+2), dtype=FLOAT_DTYPE)
    mm = numpy.zeros(shape = (nc+1, 2+1), dtype=INT_DTYPE)
    D = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    
    for i in range(1, n+1):
        h = partialEval[i]
        if h > 0.0 and w[i] > 0.0: 
            wh = w[i]*h
            k = int(X[i]+.1)
            mm[k,2] = mm[k,2]+1   #count number of categories 
            sp1[3,k] = sp1[3,k] + wh
            sp1[2,k] = sp1[2,k] + wh*(y[i]-yb)
            sp1[1,k] = sp1[1,k] + wh*h
            sp[1:kr+1, k] = sp[1:kr+1, k] + wh*db[i,1:kr+1]

    mm[1:nc+1, 1] = range(1,nc+1) 
    bof0 = big
    ns = 0
    jj = nc
    nrt = 0
    k1 = 1
    EXECUTE = True
    while EXECUTE:
        bof1 = big
        js = 0
        for j in range(1, jj+1):
            k = mm[j,1]
            if mm[k,2] != 0: 
                nr = nrt + mm[k,2]
                if nr > me and ntt-nr > me: 
                    dy = sp1[2,n1] + sp1[2,k]
                    a  = sp1[3,n1] + sp1[3,k]
                    dv = sp1[1,n1] + sp1[1,k] - a**2/sw
                    if dv > 0.0:
                        D[1:kr+1] = sp[1:kr+1, n1] + sp[1:kr+1, k]
                        a = numpy.inner(D[1:kr+1], DY[1:kr+1])
                        b = numpy.inner(D[1:kr+1], D[1:kr+1])
                        b = dv - b
                        #debug.info('b, eps, dv' +repr((b, eps, dv)))
                        if b > eps*dv: 
                            nop = nop + 1
                            b = -(dy-a)**2/b
                            #debug.info('b, bof0, bof1' +repr((b, bof0, bof1)))
                            #these checks are very important ones in the algorithm. if it
                            #passes the new hockey stick might be accepted and the response
                            #surface changes. In the original code there was no accounting for 
                            #numerical issues. The addition of 1.e-10 was added only after
                            #some painful debugging. The effect of adding this is to say 
                            #that the new coefficient (b) is taken as the new one if it's
                            #really smaller. This check is similar for continuous variables.
                            if b + 1.e-10 < bof1:
                                bof1 = b
                                js = j
                            if b + 1.e-10 < bof0: 
                                bof0 = b
                                ns = jj
                    if nc == 2: 
                        break 
        #debug.info('js='+repr(js))
        if js != 0:
            k = mm[js,1]
            mm[js,1] = mm[jj,1]
            mm[jj,1] = k
            sp1[2,n1] = sp1[2,n1] + sp1[2,k]
            sp1[3,n1] = sp1[3,n1] + sp1[3,k]
            nrt = nrt + mm[k,2]
            sp1[1,n1] = sp1[1,n1] + sp1[1,k]
            sp[1:kr+1, n1] = sp[1:kr+1, n1] + sp[1:kr+1, k]
            jj = jj-1
        else:
            break
        EXECUTE = (jj > 2)

    coef = bof0
    defaultKnot = catKnots[variable] #kcp

    CM = numpy.zeros(nc+1)
    #for j in range(1, nc+1):
    #    cm[j+kcp] = 0.0
    if ns == 0:
        return coef, defaultKnot, CM, nop
    #print cm[kcp+1:kcp+nc+1]
    CM[mm[ns:nc+1,1]] = 1.0
    #for j in range(ns, nc+1):
    #    cm[mm[j,1]+kcp] = 1.0
    #debug.info('in csp, kcp, ns, nc, mm='+repr((kcp, ns, nc,mm[ns:nc+1, 1])))
    #debug.info('in csp, cm='+repr((cm[kcp+1:kcp+nc+1])))
    return coef, defaultKnot, CM, nop

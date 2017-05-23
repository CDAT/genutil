#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def jfvc(l, m, tb, cm, nv, jv, jp):

    ip = m
    nv = 0
    while ip > 0:
        j = int(abs(tb[2,ip])+.1)
        EXECUTE = True
        if l == 1:
            if cm[2*j] > 0.0: 
                ip = int(tb[4,ip]+.1)
                EXECUTE = False
        else:   
            if cm[2*j] == 0.0: 
                ip = int(tb[4,ip]+.1)
                EXECUTE = False
        if EXECUTE:
            nv = nv+1
            jv[nv] = j
            if l != 1 and tb[2,ip] < 0.0:
                jv[nv] = -j
            if l != 1: 
                jp[nv] = int(tb[3,ip]+.1)
            ip = int(tb[4,ip]+.1)
    
    if nv <= 1:
       return nv, jv, jp
    j = nv-1
    ll = 1
    while ll != 0:
        ll = 0
        for i in range(1, j+1):
            if abs(jv[i]) > abs(jv[i+1]): 
                ll = 1
                k = jv[i]
                jv[i] = jv[i+1]
                jv[i+1] = k
                if l != 1:
                    k = jp[i]
                    jp[i] = jp[i+1]
                    jp[i+1] = k
    return nv, jv, jp

#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def orgpl(xm, xs, nk, tb, cm):
      
    for m in range(1, nk+1):
        j = int(abs(tb[2,m])+.1)
        if cm[2*j] <= 0.0: 
            tb[3,m] = xm[j] + xs[j]*tb[3,m]
    for m in range(1, nk+1):
        if tb[1,m] != 0.0:
            scl = 1.0
            ip = m
            while ip > 0:
                j = int(abs(tb[2,ip])+.1)
                if cm[2*j] == 0.0:
                    scl = scl*xs[j]
                ip = int(tb[4,ip]+.1)
            tb[1,m] = tb[1,m]/scl
    return tb

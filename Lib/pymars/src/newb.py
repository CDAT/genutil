#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from .ieq import ieq
def newb(m, tb):
    #this looks like code from marsgo
    newb = 0
    for k in range(1, m):
        if ieq(tb[2,k],tb[2,m],1.0) != 0: 
            if ieq(tb[3,k],tb[3,m],1.0) != 0: 
                if ieq(tb[4,k],tb[4,m],1.0) != 0: 
                    newb = 1
                    break 
    return newb
def isNewHS(x, y):
    a,n=y.shape
    for i in range(n):
        if (abs(x-y[:,i]) <= 1.e-5).all():
            return False
    return True
    
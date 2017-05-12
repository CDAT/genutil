#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.pr import pr
def side(l, nt, jv, xe, x):
    #print 'side'
    l2 = l+l
    l3 = l2+l
    l4 = l3+l
    for k in range(1, l+1):
        xl = xe[1,jv[k]]
        xr = xe[2,jv[k]]
        for j in range(1, nt+1):
            z = x[j,k]
            if z <= xl: 
                x[j,k+l] = xl
                x[j,k+l2] = x[j,k+l]
                x[j,k+l3] = 0.0
                x[j,k+l4] = x[j,k+l3]
            else:   
                dl = z-xl
                dr = xr-z
                x1 = xl
                x2 = xr
                for m in range(1, nt+1):
                    a = x[m,k]
                    if a != z: 
                        dx = a-z
                        if dx < 0.0 and -dx < dl: 
                            dl = -dx
                            x1 = a
                        if dx > 0.0 and dx < dr: 
                            dr = dx
                            x2 = a
                x1 = 0.5*(x1+z)
                x2 = 0.5*(x2+z)
                if x[j,k+l] > 0.0: 
                    x[j,k+l] = x1
                    x[j,k+l2] = x2
                else: 
                    x[j,k+l] = x2
                    x[j,k+l2] = x1
                x[j,k+l3], x[j,k+l4] = pr(x[j,k+l], x[j,k], x[j,k+l2], x[j,k+l3], x[j,k+l4])

    return x

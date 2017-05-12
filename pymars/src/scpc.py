#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def scpc(xm, xs, jp, l, nt, jv, tc, b):

    l2 = l+l
    l3 = l2+l
    l4 = l3+l
    q = 1.0
    for k in range(1, l+1):
        j = jv[k]
        g = xm[j]
        h = xs[j]
        q = q*h
        tc[jp,k+l] = g+h*tc[jp,k+l]
        tc[jp,k] = g+h*tc[jp,k]
        tc[jp,k+l2] = g+h*tc[jp,k+l2]
        tc[jp,k+l3] = tc[jp,k+l3]/h
        tc[jp,k+l4] = tc[jp,k+l4]/h**2
    b = b/q
    return tc, b

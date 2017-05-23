#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def pr(um, u, up, p, r):
    #this is where the cubic is being spliced at the knot?
    s=1.0
    if um > up: 
    	s=-1.0
    p = s*(2.0*up + um - 3.0*u)/(up-um)**2
    r = s*(2.0*u - up - um)/(up-um)**3
    return p, r
#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
from pymars.blf import blf
def blf0(NEST, l, ja, n, x ,w, cm, sc):
    nnt=0
    bl = blf(l, n, sc)
    if ja > 0:
    	bl = NEST.cmpnst(ja, n, x, cm, bl)

    for i in range(1, n+1):
        if bl[i] > 0.0 and w[i] > 0.0:
      		nnt = nnt+1
    return nnt, bl

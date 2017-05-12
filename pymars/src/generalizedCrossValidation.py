#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
def generalizedCrossValidation(a,b,c):
    gcv = a/(1.0 - b/c)**2
    return gcv
def logisticGCV(n, ef, w, sw, y, Y):
    s = 0.0
    t = s
    for i in range(1, n+1):
        yh = 1.0/(1.0 + numpy.exp(-Y[i]))
        gcv = ef*(y[i]-yh)**2
        s = s + w[i]*gcv
        t = t + w[i]*yh*(1.0-yh)
    s = s/sw
    t = t/sw
    return s, t
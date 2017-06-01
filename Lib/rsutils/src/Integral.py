#from quadpack import *
from Numeric import *
from math import *

def GaussIntegral(a, b):
#compute the integral of the Gaussian distribution
#mean 0 & standard deviation 1
    epsabs = 1.0e-10
    epsrel = 1.0e-10
    key = 6
    limit = 400
    last = 0
    lenw = 10*limit

    iwork = array(limit*[0.], Int32)
    work = array(lenw*[0.], Float64)

    result, abserr, neval, ier, last = dqag(a, b, epsabs, epsrel, key, 
                                    limit, lenw, iwork, work)
                            
    return result[0]
    
def gauinv(p):
    #Algorithm AS 70 Applied Statistics, 1974, vol 23 no 1
    #gauinv finds the percentage points of the normal distribution
    #It computes the inverse of the Normal CDF
    
    zero = 0.0
    one = 1.0
    half =.5
    alimit = 1.0e-20
    
    p0 = -.322232431088
    p1 = -1.0
    p2 = -.342242088547
    p3 = -.204231210245e-1
    p4 = -.453642210148e-4
    q0 = .993484626060e-1
    q1 = .588581570495
    q2 = .531103462366
    q3 = .103537752850
    q4 = .38560700634e-2
    
    ifault = 1
    cutoff = zero
    ps=p
    if ps > half:
        ps = one - ps
    if ps < alimit:
        return (cutoff, ifault)
    ifault = 0
    if ps == half:
        return (cutoff, ifault)
    yi = sqrt(log(one/ps**2))
    cutoff = yi +  ((((yi*p4 + p3)*yi + p2)*yi + p1)*yi + p0) \
                 / ((((yi*q4 + q3)*yi + q2)*yi + q1)*yi + q0)

    if p < half:
        cutoff = -1*cutoff
    return (cutoff, ifault)

def alnorm(x, upper):
    #Algorithm 66 Applied Statistics, 1973, 22, no 3
    #evaluates the tail area of the standard normal curve
    #from x to infinity if upper is true or from
    #-infinity to x if upper is false
    
    ltone = 7.0
    utzero = 18.66
             
    zero = 0.0
    one = 1.0
    half =.5
    con  = 1.28
    
    up = upper
    z = x
    if not (z >= zero):
        up = not up
        z = -z
    if not (z <= ltone or (up and z <= utzero)):
        integral = zero
        if not up:
            integral = one - integral
        return integral
    
    y = half*z**2
    if not (z > con):
        integral = half - z*(.39894228044 - .399903438504 * y \
                        / (y + 5.75885480458 - 29.8213557808 \
                        / (y + 2.62433121679 + 48.6959930692 \
                        / (y + 5.92885724438))))
    else:
        integral = .398942280385*exp(-y) \
                 / (z - 3.8052e-8 + 1.00000615302 \
                 / (z + 3.98064794e-4 + 1.98615381364 \
                 / (z - .151679116635 + 5.29330324926 \
                 / (z + 4.8385912808 - 15.1508972451 \
                 / (z + .742380924027 + 30.789933034 \
                 / (z + 3.99019417011))))))
    if not up:
        integral = one -integral
     
    return integral
from VVutil.Statistics import *
import scipy, scipy.stats
from numpy import *

def computeEDF(population, samples):
#compute the empirical distribution function
    d = 1.0/len(samples)
    edf = [(-scipy.Infinity, 0.0)]
    total = 0.0
    for x in population:
        total = total + d*samples.count(x)
        edf = edf + [(x, total)]
    return edf
def computeDstat(f1, f2):
#compute the absolute difference of f2 relative to f1
    Dstat = -scipy.Infinity

    #form intervals of f2
    Intervals = []
    (low, y) = f2[0]
    for (x, y) in f2[1:]:
        Intervals = Intervals + [(low, x)]
        low = x
    Intervals = Intervals + [(x, scipy.Infinity)]

    lastx = -scipy.Infinity
    lasty = 0.0
    for (x, y) in f1[1:]:
        j = 0
        #search for interval containing x
        keepSearching = True
        while keepSearching:
            (low, high) = Intervals[j]
            keepSearching = not(low <= x and x < high)
            if keepSearching:
                j = j+1
       
        (x2, y2) = f2[j]
        d = y2 - y

        if d >= 0 and low <= lastx:
            d=y2 - lasty
        lastx = x
        lasty = y
        
        d = abs(d)
        Dstat = max(d, Dstat)
        #print x, y, y2, d

    return Dstat

edf = [(-scipy.Infinity, 0.0), (1.0, 10.), (10.0, 11.)]
cdf = [(-scipy.Infinity, 0.0), (2.0, 2.), (3., 3.), (4., 4.), (5.,5.), (11., 6.)]
print 'edf=', edf
print 'cdf=', cdf
Dstat = computeDstat(cdf, edf)
print 'Dstat=', Dstat

Dstat = computeDstat(edf, cdf)
print 'Dstat=', Dstat

print 'check normal samples'
nsamples = 40 
samples = scipy.stats.norm.rvs(size=nsamples).tolist()
samples.sort()
xedf = findPopulation(samples) 
edf = computeEDF(xedf, samples)

ycdf = scipy.stats.norm.cdf(array(xedf)) 
cdf=[(-scipy.Infinity, 0.0)]
for j in range(nsamples):
    cdf = cdf + [(xedf[j], ycdf[j])]
    
print 'edf=', edf
print 'cdf=', cdf
Dstat = computeDstat(cdf, edf)
ks = scipy.special.smirnov(nsamples, Dstat) 

#the confidence quoted here is the probability that the 
#two distributions are the same. See
#Pratical nonparametric statistics by Conover, p 428 
conf = 1. - (1. - ks)**2  #=Probability(D > Dstat) 
                             #where Dstat is the maximal absolute difference
                             #between the empirical distribution and the
                             #proposed distribution; in this case normal

print 'Dstat=', Dstat, 'KS=', ks, 'confidence=', conf
print 'mu=', scipy.mean(samples), 'std=', scipy.std(samples)
x = Dstat
x = 9./40
x=10./40
n = nsamples
nx = int(modf(n*x)[1])

num = 1.0
for i in range(nx):
    num = num*(n-i)
den = 1.0
for i in range(nx):
    den = den*(n+i+1)
print 1-(num/den)


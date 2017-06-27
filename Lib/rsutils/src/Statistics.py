from math import *
import scipy
def findPopulation(samples):
    population = []
    for sample in samples:
        if population.count(sample) == 0:
            population = population + [sample]
    return population
def computeStats(samples):
    nsamples = len(samples)
    sum = 0.0
    for sample in samples:
        sum = sum + sample
    mu = sum/nsamples
    
    sum = 0.0
    for sample in samples:
        sum = sum + (sample - mu)**2
    sd = sqrt(sum/(nsamples-1))
    
    return (mu, sd)
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

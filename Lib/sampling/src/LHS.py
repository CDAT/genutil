import random, sys
from math import *
from Roots import *

def genUniformLHSdata(npts, box):

    dim = len(box)
    samples =[]
    if npts == 0:
        return LHSdata
        
    for (low, high) in box:
        delta = (high - low)/npts
        randInd = randomIndices(npts)
        pnts = []
        for i in randInd:
            L = low + i*delta
            H = L + delta
            pnt = (H - L)*random.random() + L
            pnts = pnts + [pnt]
        samples = samples + [pnts]
            
    return samples  
def genGeometricLHSdata(npts, box, degree):
    
    dim = len(box)    
    samples = []
    
    for (low, high) in box:
        randInd = randomIndices(npts)
        pnts = []
        
        midpnt = (high + low)/2.
        delta = (high - low)/8
        middle = []
        midlow = midpnt  
        if 2*(npts/2) != npts:
            midlow = midpnt - delta
            midhigh = midpnt + delta
            middle = [(midlow, midhigh)]
            
        lowintervals = makeIntervals(npts/2, (low, midlow), degree)
                
        #reflect the low intervals about the midpoint
        highintervals = []
        for (a,b) in lowintervals:
            highintervals = [(2*midpnt - b, 2*midpnt - a)] + highintervals
        
        intervals = lowintervals + middle + highintervals

        for i in randInd:
            (L, H) = intervals[i]
            pnt = (H - L)*random.random() + L
            pnts = pnts + [pnt]
        samples = samples + [pnts]
            
    return samples 
def randomIndices(npts):

    r=[]
    for i in range(npts):
        r = r + [(i, random.random())]
        
    r.sort(lambda x, y: cmp(x[1], y[1]))
    
    indices = []
    for pair in r:
        indices = indices + [pair[0]]
    return indices             
def makeIntervals(nint, (low, high), degree):
    
    delta = (high - low)/(nint**degree)
    poly = [1.0 - (high - low)/delta]
    for i in range(nint-1):
        poly = poly + [1.0]
    epsilon = Newton.findRoot(poly, (low, high), .0001)

    intervals = []
    a=low
    width = delta
    for i in range(nint):
        b = a + width
        intervals = intervals + [(a,b)]
        a=b
        width = epsilon*width
        
    return intervals    
import random, sys
from math import *
def genMC(npts, box):
    
    MCdata =[]
    for (low, high) in box:
        pnts = []
        for i in range(npts):
            pnt = (high - low)*random.random() + low
            pnts = pnts + [pnt]
        MCdata = MCdata + [pnts]
            
    return MCdata           
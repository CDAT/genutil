#!/usr/local/bin/python

import Sampling
from MARStools.ResponseSurfaceTools import *
from math import *

class SamplingController(object):
    def __init__(self, npts=1000, dim=3, low=-2.0, high=2.0, type='GeometricLHS', degree=2, output='samples.dat'):
        self.__npts = npts
        self.__dim = dim
        self.__low = low
        self.__high = high
        self.__type = type
        self.__degree = degree
        self.__output = output
        self.__status = 0

    def getStatus(self):
        return self.__status
              
    def run(self):
        box = genBox(self.__dim, self.__low, self.__high)
        out_filename = self.__type + self.__output
        if self.__type == 'GeometricLHS':           
            s = Sampling.LHS.genGeometricLHSdata(self.__npts, box, self.__degree)
        if self.__type == 'UniformLHS':
            s = Sampling.LHS.genUniformLHSdata(self.__npts, box)
        if self.__type == 'MC':
            s = Sampling.MonteCarlo.genMC(self.__npts, box)
            
        try:
            f = open(out_filename, 'w')
            f.write(str(s))
            self.__status = 1
        except IOError:
            print "Error writing %s " % out_filename
            
        f.close()
        
################################################################################   

if __name__=="__main__":        
            
    sc = SamplingController()
    sc.run()
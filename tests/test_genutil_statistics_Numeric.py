#!/usr/bin/env python

import numpy,genutil
import unittest

class GENUTIL(unittest.TestCase):
    def assertArraysEqual(self,A,B):
        self.assertTrue(numpy.all(numpy.equal(A,B)))
    def testStatisticsNumpy(self):
        a=numpy.ones((15,25),'d')
        rk = [0.0, 91.66666666666667, 87.5, 83.33333333333333, 79.16666666666667, 75.0, 70.83333333333333, 66.66666666666667, 62.5, 58.333333333333336, 54.166666666666664, 95.83333333333333, 50.0, 41.666666666666664, 37.5, 33.333333333333336, 29.166666666666668, 25.0, 20.833333333333332, 16.666666666666668, 12.5, 8.333333333333334, 4.166666666666667, 45.833333333333336, 100.0]
        # rk will be copied over and over
        self.assertArraysEqual(genutil.statistics.rank(a,axis=1),rk)
        self.assertTrue(numpy.allclose(genutil.statistics.variance(a,axis=0),0.))

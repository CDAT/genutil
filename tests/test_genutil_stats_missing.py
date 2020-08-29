# Adapted for numpy/ma/cdms2 by convertcdms.py
import genutil, MV2
import unittest
import numpy

class GENUTIL(unittest.TestCase):
    def testError(self):
        error = genutil.StatisticsError("Custom error message")

        self.assertTrue(str(error) == "Custom error message")

    def testStatsMissing(self):
        a=MV2.array([1,2,3,4,5],mask=[0,0,1,0,0])
        self.assertTrue(numpy.allclose(genutil.statistics.std(a),1.58113883008))
        self.assertTrue(numpy.ma.masked_all(genutil.statistics.std(a,max_pct_missing=10.)))
        self.assertTrue(numpy.allclose(genutil.statistics.std(a,max_pct_missing=79),1.58113883008))

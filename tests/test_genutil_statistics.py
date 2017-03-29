#!/usr/bin/env python
# Adapted for numpy/ma/cdms2 by convertcdms.py

import cdms2,genutil,os,sys
import unittest
import cdat_info
import numpy

class GENUTIL(unittest.TestCase):
    def assertArraysEqual(self,A,B):
        self.assertTrue(numpy.all(numpy.equal(A,B)))

    def testStatisitcs(self):
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'clt.nc'))

        u=f('u',time=slice(0,1),level=slice(0,1),squeeze=1)
        v=f('v',time=slice(0,1),plev1=slice(0,1),squeeze=1)

        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),"genutil_statistics.nc"))
        print 'Lagged correlation'
        nm = "lagged_corr"
        self.assertArraysEqual(genutil.statistics.laggedcorrelation(u,v,axis=0), f( nm+"_1"))
        self.assertArraysEqual(genutil.statistics.laggedcorrelation(u,v,axis=0,lag=4), f( nm+"_2"))
        self.assertArraysEqual(genutil.statistics.laggedcorrelation(u,v,axis=0,lag=4,noloop=1), f( nm+"_3"))
        self.assertArraysEqual(genutil.statistics.laggedcorrelation(u,v,axis=0,lag=[4,8,10]), f(nm+"_4"))


        print 'Lagged covariance'
        nm = "lagged_cov"
        self.assertArraysEqual(genutil.statistics.laggedcovariance(u,v,axis=0), f(nm+"_1"))
        self.assertArraysEqual(genutil.statistics.laggedcovariance(u,v,axis=0,lag=4), f(nm+"_2"))
        self.assertArraysEqual(genutil.statistics.laggedcovariance(u,v,axis=0,lag=4,noloop=1), f(nm+"_3"))
        self.assertArraysEqual(genutil.statistics.laggedcovariance(u,v,axis=0,lag=[4,8,10]), f(nm+"_4"))

        print 'Auto correlation'
        nm = "auto_corr"
        self.assertArraysEqual(genutil.statistics.autocorrelation(u,axis=0), f(nm+"_1"))
        self.assertArraysEqual(genutil.statistics.autocorrelation(u,axis=0,lag=4), f(nm+"_2"))
        self.assertArraysEqual(genutil.statistics.autocorrelation(u,axis=0,lag=4,noloop=1), f(nm+"_3"))
        self.assertArraysEqual(genutil.statistics.autocorrelation(u,axis=0,lag=[4,8,10]), f(nm+"_4"))

        print 'Auto covariance'
        nm = "auto_cov"
        self.assertArraysEqual(genutil.statistics.autocovariance(u,axis=0), f(nm+"_1"))
        self.assertArraysEqual(genutil.statistics.autocovariance(u,axis=0,lag=4), f(nm+"_2"))
        self.assertArraysEqual(genutil.statistics.autocovariance(u,axis=0,lag=4,noloop=1), f(nm+"_3"))
        self.assertArraysEqual(genutil.statistics.autocovariance(u,axis=0,lag=[4,8,10]), f(nm+"_4"))
        f.close()

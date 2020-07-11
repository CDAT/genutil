from __future__ import print_function
import cdms2
import genutil
import os
import unittest
import cdat_info
import numpy

dump = False


class GENUTIL(unittest.TestCase):
    def assertArraysEqual(self, A, nm, skip=False, dump=dump):
        B = self.good(nm)
        print("A-B max difference", numpy.abs(A - B).max(), skip)
        print("A-B % difference", numpy.abs((A - B) / A).max(), skip)
        if dump:
            self.out.write(A, id=nm)
            self.out.sync()

        if not skip:
            self.assertTrue(numpy.ma.allclose(A, B))

    def testPercentiles(self):
        a = numpy.arange(1.0, 6.0)
        p = genutil.statistics.percentiles(a)

        assert (p == a).all()

        with self.assertRaises(Exception):
            p = genutil.statistics.percentiles(a, axis=1)

        aa = numpy.array([a for _ in range(10)])

        p = genutil.statistics.percentiles(aa)

        assert (p == aa).all()

        p = genutil.statistics.percentiles(aa, axis=1)

        assert (p == numpy.full((1, 10), 3.0)).all()

        p = genutil.statistics.percentiles(aa, [0, 50, 100], axis=1)

        assert (p == numpy.array([numpy.full((10,), x) for x in (1.0, 3.0, 5.0)])).all()

        a = numpy.arange(1.0, 6.0).reshape((5, 1))

        p = genutil.statistics.percentiles(a)

        assert (p == numpy.array([[3.0]])).all()

    def testStatisitcs(self):
        f = cdms2.open(os.path.join(cdat_info.get_sampledata_path(), "clt.nc"))

        u = f("u", time=slice(0, 1), level=slice(0, 1), squeeze=1)
        v = f("v", time=slice(0, 1), plev1=slice(0, 1), squeeze=1)

        self.good = cdms2.open(
            os.path.join(cdat_info.get_sampledata_path(), "genutil_statistics.nc")
        )
        if dump:
            self.out = cdms2.open("genutil_statistics_new.nc", "w")
        print("Lagged correlation")
        nm = "lagged_corr"
        self.assertArraysEqual(
            genutil.statistics.laggedcorrelation(u, v, axis=0), nm + "_1"
        )
        self.assertArraysEqual(
            genutil.statistics.laggedcorrelation(u, v, axis=0, lag=4), nm + "_2"
        )
        self.assertArraysEqual(
            genutil.statistics.laggedcorrelation(u, v, axis=0, lag=4, noloop=1),
            nm + "_3",
        )
        self.assertArraysEqual(
            genutil.statistics.laggedcorrelation(u, v, axis=0, lag=[4, 8, 10]),
            nm + "_4",
        )

        print("Lagged covariance")
        nm = "lagged_cov"
        self.assertArraysEqual(
            genutil.statistics.laggedcovariance(u, v, axis=0), nm + "_1"
        )
        self.assertArraysEqual(
            genutil.statistics.laggedcovariance(u, v, axis=0, lag=4), nm + "_2"
        )
        self.assertArraysEqual(
            genutil.statistics.laggedcovariance(u, v, axis=0, lag=4, noloop=1),
            nm + "_3",
        )
        self.assertArraysEqual(
            genutil.statistics.laggedcovariance(u, v, axis=0, lag=[4, 8, 10]), nm + "_4"
        )

        print("Auto correlation")
        nm = "auto_corr"
        self.assertArraysEqual(genutil.statistics.autocorrelation(u, axis=0), nm + "_1")
        self.assertArraysEqual(
            genutil.statistics.autocorrelation(u, axis=0, lag=4), nm + "_2"
        )
        self.assertArraysEqual(
            genutil.statistics.autocorrelation(u, axis=0, lag=4, noloop=1), nm + "_3"
        )
        self.assertArraysEqual(
            genutil.statistics.autocorrelation(u, axis=0, lag=[4, 8, 10]), nm + "_4"
        )

        print("Auto covariance")
        nm = "auto_cov"
        self.assertArraysEqual(genutil.statistics.autocovariance(u, axis=0), nm + "_1")
        self.assertArraysEqual(
            genutil.statistics.autocovariance(u, axis=0, lag=4), nm + "_2"
        )
        self.assertArraysEqual(
            genutil.statistics.autocovariance(u, axis=0, lag=4, noloop=1), nm + "_3"
        )
        self.assertArraysEqual(
            genutil.statistics.autocovariance(u, axis=0, lag=[4, 8, 10]), nm + "_4"
        )
        f.close()

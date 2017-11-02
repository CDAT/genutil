#!/usr/bin/env python

import cdms2,genutil,cdtime,os,sys
import unittest
import cdat_info
import numpy

cdms2.setAutoBounds('on')
class GENUTIL(unittest.TestCase):
    def testPicker(self):
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'ta_ncep_87-6-88-4.nc'))
        levels = [1000,700,800]
        with self.assertRaises(Exception):
            s=f('ta',slice(0,1),genutil.picker(level=levels,match=1))


        s=f('ta',slice(0,1),genutil.picker(level=levels,match=0))

        self.assertEqual(s.shape[1],3)
        self.assertFalse((s.getLevel()[:]!=levels).any(),"Error did not retrieve the right levels!")

        self.assertTrue(numpy.all(s[0,-1].mask))

        levels = [1000,700,850]
        s3=f('ta',genutil.picker(time=['1987-7','1988-1',cdtime.comptime(1988,3)],level=[1000,700,850]))

        self.assertEqual(s3.shape,(3, 3, 73, 144))
        t1= cdtime.componenttime(1987,7)
        t2= cdtime.componenttime(1988,1)
        t3= cdtime.componenttime(1988,3)
        tc = s3.getTime().asComponentTime()
        for i,g in [t1, t2, t3]:
            self.assertEqual(tc[i].cmp(g),0)
        test = s3.getLevel()[:]!=levels
        self.assertFalse(test.any())

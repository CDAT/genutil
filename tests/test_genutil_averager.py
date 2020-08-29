# Adapted for numpy/ma/cdms2 by convertcdms.py
"""This is a set of tests for averager - which can be used to calculate area weighted averages."""
import cdms2, MV2, numpy, numpy.ma, cdtime, os, sys, cdat_info
from genutil import averager, AveragerError
cdms2.setAutoBounds('on')
import unittest

class GENUTIL(unittest.TestCase):
    def testError(self):
        error = AveragerError("Custom error message")

        print(str(error))

        self.assertTrue(str(error) == "Custom error message")

    def testAverager(self):
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'tas_ukmo_con.nc'))
        x = f('tas')
        ans = averager(x, axis='yt', weights = ['weighted','equal'], combinewts=1)
        #


        ans2 = averager(x, axis='yx', weights = [numpy.ones(len(x.getLatitude()), numpy.float), 'generate'], combinewts=1)


        b = numpy.array([1.,2.3, 3.0, 4.3])
        b1 = numpy.ma.array(b)

        b2,w = numpy.ma.average(b1, weights=None, returned=1)
        self.assertTrue( numpy.ma.allclose(b2, 2.65))

        b2b,wb = averager(b1, weights=None, returned=1)
        self.assertTrue( numpy.ma.allclose(b2b, 2.65))
        self.assertTrue( numpy.ma.allclose(wb, w))


        ba = averager(b)
        self.assertTrue( numpy.ma.allclose(ba, 2.65))


        bs = averager(b, action='sum')
        self.assertTrue( numpy.ma.allclose(bs, 10.6))



        b = numpy.ones((3,4,5))
        c = averager(b, weights=None, axis=[0,2], action='sum')
        self.assertTrue( numpy.ma.allclose(c, [ 15., 15., 15., 15.,]))
        self.assertTrue( isinstance(c, numpy.ndarray))

        with self.assertRaises(AveragerError):
            c2 = averager(b, weights=numpy.ones(b.shape), axis=[0,2], action='sum')


        with self.assertRaises(AveragerError):
            c2 = averager(b, weights=numpy.ones((3,5,4)), axis=[0,2], action='sum')


        d = averager(b, weights='equal', axis=2, action='sum')

        d0MA = numpy.ma.average(b, weights=numpy.ma.array([1.,2.,3.,4.,5.]), axis=2)
        d0 = averager(b, weights=numpy.ma.array([1.,2.,3.,4.,5.]), axis=2) 
        self.assertTrue( numpy.ma.allclose(d0, d0MA))


        d1 = averager(b, weights=numpy.ma.array([1.,2.,3.,4.,5.]), axis=2, action='sum')

        d2 = averager(b, weights=[numpy.ma.array([1.,2.,3.]),'equal'], axis=[0,2], action='sum')


        b = numpy.array([1.,2.3, 3.0, 4.3])
        bs = averager(b, action='sum')
        self.assertTrue( numpy.ma.allclose(bs, 10.6)  )


        b = numpy.ones((3,4,5))
        c = averager(b, weights=None, axis=[0,2], action='sum')
        self.assertTrue( numpy.ma.allclose(c, [ 15., 15., 15., 15.,]))
        self.assertTrue( isinstance(c, numpy.ndarray))

        b1 = numpy.ma.array(b)
        c1 = averager(b1, weights=None, axis=[0,2])
        self.assertTrue( isinstance(c1,numpy.ndarray))

        averager(x, axis = 'xyt')

        ans = averager(x)

        averager(x, axis='2t')

        averager(x, axis='t(lat)')

        with self.assertRaises(AveragerError):
            averager(x, axis='t(lev)')


        with self.assertRaises(AveragerError):
            averager(x, axis='t3')

        averager(x, axis='t')

        averager(x, axis='tx', weight=['generate', 'generate']) 

        averager(x, axis='tx', weight=['equal', 'generate']) 

        with self.assertRaises(AveragerError):
            averager(x, axis='tx', weight='equal') 

        with self.assertRaises(AveragerError):
            a = numpy.array(list(range(10)), numpy.float)
            averager(x, axis='tx', weight=['equal', a]) 

        with self.assertRaises(AveragerError):
            b=numpy.ma.array(a)
            averager(x, axis='tx', weight=['equal', b]) 

        result = averager(x, axis='tx', weight=['equal', 'equal'])

        result = averager(x, axis='2t', weight=['generate', 'equal'])


        #**********************************************************************
        #
        # Create the area weights 
        #
        faw=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),"area_weights.nc"))
        aw = faw("aw")
        aw.setAxisList(x.getAxisList())
        
        #
        #
        #**********************************************************************


        result = averager(x, axis='x', weight=aw)


        result = averager(x, axis='xy', weight=aw) 


        #
        # Now I want the Longitude axis to be area weighted (including any missing data)
        # but the latitude axis to be equally weighted
        #
        result, newwts = averager(x, axis='x', weight=aw, returned=1) 
        new_result = averager(result, axis='y', weight='equal')
        result = averager(x, axis='21', weight=aw) 

        result = averager(x, axis=0) 

        #******************************************************************************************
        # Now a real world check!!
        #******************************************************************************************

        #
        # The values in this file were calculated using old routines in CDAT2.4
        #
        fcheck = cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'tas_gavg_rnl_ecm.nc'))
        start = cdtime.componenttime(1979, 1, 1)
        end =  cdtime.componenttime(1979, 12, 1)
        correct = fcheck('tas', time=(start, end))

        #
        # The source for the averages in the above file is the raw data in the file below.
        # The difference in units (degrees C and degrees K) should be the obvious answer.
        #
        f = cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'tas_ecm_1979.nc'))

        #
        # I can invoke averager by passing the f('tas') as the variable instead of first extracting
        # the variable...... Youcan obviously get more fancy with the kind of selected variables you pass.
        #
        result = averager(f('tas'), axis='xy', weight=['generate', 'generate'])


        #
        # For the purposes of this test, I am using the extracted variable as below
        #
        x = f('tas')


        #
        # I am going to calculate the answer in 2 ways.
        #
        # First method: Use the 'generate' options to average over x and y
        #
        result1 = averager(x, axis='xy', weight=['generate', 'generate'])

        #
        # Second Method: Create the area weights and
        #                convert it into an MV2 before passing to the averager.
        #
        aw2 = faw("aw2")
        aw2.setAxisList(x.getAxisList())
        result2 = averager(x, axis='x(lat)', weight=aw2)


        #
        # Third Method: Use the area weights from Charles' function created above
        #               but call the averaging only one 1 dimension at a time. Use the average &
        #               weights from the first step in the averaging at the second step!.
        #
        temp_step, temp_wts = averager(x, axis='x', weight=aw2, returned=1)
        result3 = averager(temp_step, axis='y', weight=temp_wts)

        #
        # Note that the above way of doing multiple steps is useful when you want the temp_wts
        # altered between steps.......

        #
        # Now check the 3 results....... they will be different by 273.0 (Centigrade to Kelvin)
        #
        diff1 = result1 - correct
        self.assertTrue( MV2.allclose(273.0, diff1))
        diff2 = result2 - correct
        self.assertTrue( MV2.allclose(273.0, diff1))
        diff3 = result3 - correct
        self.assertTrue( MV2.allclose(273.0, diff1))
        self.assertTrue( MV2.allclose(result1, result2))
        self.assertTrue( MV2.allclose(result2, result3))

        #
        # This test is to verify the action='sum' option
        #
        tasm_file = os.path.join(cdat_info.get_sampledata_path(),'tas_cru_1979.nc')

        ftasm = cdms2.open(tasm_file)
        xtasm = ftasm('tas')
        ywt = faw("ywt")
        ywt.setAxisList(xtasm.getAxisList())
        #
        # This is a good way to compute the area fraction that the data is non-missing
        #
        ysum2 = averager(ywt, axis='xy', weight=['equal', 'equal'], action='sum')

        #
        # Verification of combine weights for accuracy in the presence of missing data
        #
        xavg_1 = averager(xtasm, axis = 'xy', weights = ywt)
        xavg_2 = averager(xtasm, axis = 'xy', weights = ['generate', 'generate', 'generate'], combinewts=1)
        self.assertTrue( MV2.allclose(xavg_1, xavg_2))

        #
        # Real world Averager Test #2
        #
        newf = cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'clt.nc'))
        u = newf('u')
        u2 = averager(u, axis='1')
        u3 = averager(u, axis='(plev)')
        self.assertTrue( numpy.ma.allclose(u2, u3))

        uw = faw("uw")
        uw.setAxisList(u.getAxisList())
        u4 = averager(u, axis='1x', weight=uw)
        #
        # SUM and AVERAGE should return the same answer when dimensions are averaged
        # together or individually (in the same order)if the data has no missing values.
        # Test this!
        #
        clt = newf('clt')

        # First try the average
        clt_ave = averager(clt, axis='txy')
        clt1a = averager(clt, axis='t')
        clt2a = averager(clt1a, axis='x')
        clt3a = averager(clt2a, axis='y')
        self.assertTrue( numpy.ma.allclose(clt_ave, clt3a))

        # Now try the sum
        clt_sum = averager(clt, axis='txy', action='sum')
        clt1 = averager(clt, axis='t', action='sum')
        clt2 = averager(clt1, axis='x', action='sum')
        clt3 = averager(clt2, axis='y', action='sum')
        self.assertTrue( numpy.ma.allclose(clt_sum, clt3))
        print("ALL RAN!")
        f.close()
if __name__=="__main__":
    GENUTIL().testAverager()

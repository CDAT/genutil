# Adapted for numpy/ma/cdms2 by convertcdms.py
# Adapted for numpy/ma/cdms2 by convertcdms.py
import numpy
import genutil
import cdms2
import numpy.ma
import os
import sys
import unittest
import cdat_info

class GENUTIL(unittest.TestCase):
    ### EXTRACT TESTS
    def assertArraysEqual(self,A,B):
        self.assertTrue(numpy.all(numpy.equal(A,B)))

    def test_arrayindexing_1D(self):
        ### 1D
        A=numpy.array([6,7,8,9,2],'f')
        B=numpy.argsort(A).astype('i')
        C=genutil.arrayindexing.get(A,B)
        self.assertArraysEqual(C,numpy.array([ 2.,  6.,  7.,  8.,  9.]))
        D=numpy.array([6.5,7.5,8.5,9.5,2.5],'f')
        # Sets values of D into A at indices defined in B
        E=genutil.arrayindexing.set(A,B,D)
        self.assertArraysEqual(E,numpy.array([ 7.5,  8.5,  9.5,  2.5,  6.5]))

    def test_arrayindexing_2D(self):
        ## 2D
        A=numpy.array([[1,2],[3,4,],[5,6],[7,8]],numpy.float)
        B=numpy.array([3,2],numpy.int32)

        # Extract from A at indices specified by B
        # From C directly
        C=genutil.array_indexing.extract(A,B)
        self.assertEqual(C.dtype.char,'d')
        self.assertArraysEqual(C,[7.,6.])
        # From python interface
        C=genutil.arrayindexing.get(A,B)
        self.assertEqual(C.dtype.char,'d')
        self.assertArraysEqual(C,[7.,6.])

        #### Set tests
        V=numpy.array([1345,34],A.dtype.char)
        B=numpy.array([-3,2],numpy.int)
        # Checks setting negative indices
        A=genutil.arrayindexing.set(A,B,V)
        self.assertArraysEqual(A,[[  1.,   2.],
             [  1345,   4.],
              [  5.,   34],
               [  7.,   8.]])

        A=numpy.array([[1,2],[3,4,],[5,6],[7,8]],numpy.float)
        B=numpy.array([[1,2],[3,0,],[1,2],[0,3]],numpy.int)
        V=numpy.array([[10.,21.],[13,.4,],[1.5,6.4],[77.7,9.8]],numpy.float)
        C=genutil.arrayindexing.set(A,B,V)
        self.assertArraysEqual(C,[[ 77.7,   0.4],
             [  1.5,   4. ],
              [  5. ,   6.4],
               [ 13. ,   9.8]])

        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'clt.nc'))
        clt=f('clt')
        ## clt=cdms2.MV2.average(clt,2)

        M=numpy.ma.maximum.reduce(clt,axis=0)
        marg=numpy.ma.argmax(clt,axis=0)
        M2=genutil.arrayindexing.get(clt,marg)
        self.assertArraysEqual(M2,M)

        M=numpy.ma.maximum.reduce(clt,axis=1)
        marg=numpy.ma.argmax(clt,axis=1)
        marg=cdms2.MV2.array(marg)
        marg.setAxis(0,clt.getAxis(0))
        marg.setAxis(1,clt.getAxis(2))
        M2=genutil.arrayindexing.get(clt,marg,axis=1)
        self.assertArraysEqual(M2,M)

        clt=cdms2.MV2.masked_greater(clt,80)
        M=numpy.ma.maximum.reduce(clt,axis=1)
        marg=numpy.ma.argmax(clt,axis=1)
        marg=cdms2.MV2.array(marg)
        marg.setAxis(0,clt.getAxis(0))
        marg.setAxis(1,clt.getAxis(2))
        M2=genutil.arrayindexing.get(clt,marg,axis=1)
        self.assertArraysEqual(M2,M)

    def test_arrayindexing_3D(self):
        ## 3D
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'clt.nc'))
        clt=f('clt')

        # Checks we can extract, not happy with random though...
        I=numpy.random.random(clt.shape)*clt.shape[0]
        I=I.astype('i') # integers required
        M2=genutil.arrayindexing.get(clt,I)




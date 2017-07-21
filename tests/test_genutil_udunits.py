from __future__ import print_function
import genutil,numpy
import unittest

class GENUTIL(unittest.TestCase):
    def testUdunits(self):

        m=genutil.udunits(5,'m')

        cm=genutil.udunits(7,'cm')

        i=genutil.udunits(3,'in')

        f=genutil.udunits(5,'feet')


        dC=genutil.udunits(33,'degC')
        dF=genutil.udunits(83,'degF')
        dK=genutil.udunits(300,'K')

        o = m+cm
        print(m,'+',cm,'=',o)
        self.assertEqual(o.units,m.units)
        self.assertTrue(numpy.allclose(o.value,5.07))


        o=cm+i
        print(cm,'+',i,'=',cm+i)
        self.assertEqual(o.units,cm.units)
        self.assertTrue(numpy.allclose(o.value,14.62))
        o=i+cm
        print(i,'+',cm,'=',i+cm)
        self.assertEqual(o.units,i.units)
        self.assertTrue(numpy.allclose(o.value,5.75590551181))

        m.units='km'
        print(m)
        self.assertTrue(numpy.allclose(m.value,5E-3))

        # Make sure incompatible units can't be added
        with self.assertRaises(Exception):
            m2=m+dK

        #Just to be sure test against another system
        o=dK-dC
        print(dK,'-',dC,'=',dK-dC)
        self.assertEqual(o.units,dK.units)
        self.assertTrue(numpy.allclose(o.value,-6.15))

        o = dF/dC
        print(dF,'/',dC,'=',dF/dC)
        #division by compatible units result in pur float (unitless)
        self.assertTrue(isinstance(o,float))
        self.assertTrue(numpy.allclose(o,0.908096280088))
        o=dC/dF
        print(dC,'/',dF,'=',dC/dF)
        self.assertTrue(isinstance(o,float))
        self.assertTrue(numpy.allclose(o,1.16470588235))
        o=dC/m
        print(dC,'/',m,'=',dC/m)
        self.assertEqual(o.units,"degC/km")
        self.assertTrue(numpy.allclose(o.value,6600.0))


        o=dF*dC
        print(dF,'*',dC,'=',dF*dC)
        self.assertEqual(o.units,"degF*degF")
        self.assertTrue(numpy.allclose(o.value,7586.2))
        o=dC*dF
        print(dC,'*',dF,'=',dC*dF)
        self.assertEqual(o.units,"degC*degC")
        self.assertTrue(numpy.allclose(o.value,935.0))
        o=dC*m
        print(dC,'*',m,'=',dC*m)
        self.assertEqual(o.units,"degC*km")
        self.assertTrue(numpy.allclose(o.value,0.165))

        #Test simple conversion
        o=dC.to("K")
        print(dC,'to K',dC.to('K'))
        self.assertEqual(o.units,"K")
        self.assertTrue(numpy.allclose(o.value,306.15))

        print('how to go from',dF.units,'to K',dF.how('K'))
        ## dF.show()
        ## print dF.available()

        """
        kn=dF.known_units()
        for k in kn.keys():
            print k,kn[k]
        kn=dF.known_units(bytype=1)

        for k in kn.keys():
            print 'UNITS OF',k
            ln=kn[k]
            for l in ln:
                print '\t',l
        """

        # Test remaining ufuncs
        o=m**2
        print(m**2)
        self.assertEqual(o.units,"km**2")
        self.assertTrue(numpy.allclose(o.value,2.5e-05))
        o=5+m
        print(5+m)
        self.assertEqual(o.units,"km")
        self.assertTrue(numpy.allclose(o.value,5.005))
        o=5-m
        print(5-m)
        self.assertEqual(o.units,"km")
        self.assertTrue(numpy.allclose(o.value,4.995))
        o=5*m
        print(5*m)
        self.assertEqual(o.units,"km")
        self.assertTrue(numpy.allclose(o.value,0.025))
        o=m*5
        print(m*5)
        self.assertEqual(o.units,"km")
        self.assertTrue(numpy.allclose(o.value,0.025))
        o=5/m
        print(5/m)
        self.assertEqual(o.units,"1/(km)")
        self.assertTrue(numpy.allclose(o.value,1000.))
        o=m/5
        print(m/5)
        self.assertEqual(o.units,"km")
        self.assertTrue(numpy.allclose(o.value,.001))

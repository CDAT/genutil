import genutil
import numpy
import unittest
class GENUTIL(unittest.TestCase):
    def testUdunits2Additions(self):
        ## Test for udunits2 feature wrapped
        ## Author: Charles Doutriaux
        ## The units created bellow
        ## have absolutely no meaning
        ## They are just used to test that the features exposed
        ## work properly

        genutil.udunits_wrap.init()

        ## Create a new baseunits called eq
        genutil.addBaseUnit("eq")

        ## Create a dimensionless units named dimless
        genutil.addDimensionlessUnit("dimless")

        ## Created a scaled units for dimless
        genutil.addScaledUnit("pourcent",.01,"dimless")

        ## Create an offsetted units for dimless
        genutil.addOffsettedUnit("fakeCelsius",273.15,"dimless")

        ## Create mutliplied units
        genutil.addMultipliedUnits("efC","eq","fakeCelsius")
        genutil.addMultipliedUnits("efP","eq","pourcent")

        ## Create divded units
        genutil.addDividedUnits("defC","eq","fakeCelsius")
        genutil.addDividedUnits("defP","eq","pourcent")

        ## Create inverted unit
        genutil.addInvertedUnit("iefC","defC")

        ## Test scaled
        p = genutil.udunits(1,"pourcent")
        ## Test new base unit
        eq = genutil.udunits(1,"eq")

        o=p.to("dimless")
        print o
        self.assertEqual(o.units,"dimless")
        self.assertEqual(o.value , 0.01)

        fC = genutil.udunits(2,"fakeCelsius")
        o=fC.to("dimless")
        print o
        self.assertEqual(o.units,"dimless")
        self.assertTrue(numpy.allclose(o.value,275.15))
        o=fC.to("pourcent")
        print o
        self.assertEqual(o.units,"pourcent")
        self.assertTrue(numpy.allclose(o.value,27515))

        #Trying multiplied
        o=eq*fC
        print o
        self.assertEqual(o.units,"eq*fakeCelsius")
        self.assertEqual(o.value , 2)
        o=o.to("efC")
        print o
        self.assertEqual(o.units,"efC")
        self.assertEqual(o.value , 2)
        o=o.to("efP")
        print o
        self.assertEqual(o.units,"efP")
        self.assertEqual(o.value , 200)

        #Trying divided
        o=eq/fC
        print o
        self.assertEqual(o.units,"eq/fakeCelsius")
        self.assertEqual(o.value , .5)
        o=o.to("defC")
        print o
        self.assertEqual(o.units,"defC")
        self.assertEqual(o.value , .5)
        o=o.to("defP")
        print o
        self.assertEqual(o.units,"defP")
        self.assertEqual(o.value , 5E-3)

        ## Trying inverted
        o=fC/eq
        o=o.to("iefC")
        print o
        self.assertEqual(o.units,"iefC")
        self.assertEqual(o.value , 2)

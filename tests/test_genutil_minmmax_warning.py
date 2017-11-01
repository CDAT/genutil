import unittest
import genutil
import cdat_info
import os
import cdms2
import warnings

class GENUTIL(unittest.TestCase):
    def testMinMax(self):
        path=cdat_info.get_sampledata_path()
        f=cdms2.open(os.path.join(path,"clt.nc"))
        clt=f('clt', time=slice(0,1), squeeze=1)
        with warnings.catch_warnings(record=True) as w:
            mn,mx = genutil.minmax(clt,clt*100,clt/3.)
            self.assertEqual(len(w),0)
        self.assertEqual(mn,0.)
        self.assertEqual(mx,10000.0)

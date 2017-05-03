# Adapted for numpy/ma/cdms2 by convertcdms.py
import cdms2,sys,os
import cdat_info
import unittest
import genutil

class GENUTIL(unittest.TestCase):
    def testGrower(self):
        f=cdms2.open(os.path.join(cdat_info.get_sampledata_path(),'clt.nc'))
        s=f('clt')
        works = genutil.statistics.variance(s,axis='t',weights=['equal'])

import genutil,sys,cdat_info
import unittest

class GENUTIL(unittest.TestCase):
    def testReadCol(self):
        ## Test reading of ascii file organized in columns
        print 'Testing genutil.ASCII.read_col',cdat_info.get_sampledata_path()+'/test_col.asc'
        vars = genutil.ASCII.read_col(cdat_info.get_sampledata_path()+'/test_col.asc',header=4,cskip=1,idrow=True,axis=True)
        self.assertEqual(len(vars),3, 'genutil.ASCII: Error should have returned 3 variables, returned: %s\nCheck cskip option or axis=True option' % len(vars))

        ids_ok = ['TEST','TEST2','TEST4']
        for i in range(3):
            v=vars[i]
            ax = v.getAxis(0)
            self.assertEqual(v.id,ids_ok[i],'genutil.ASCII: Error axis on variables should be named axis, looks like the axis setup failed!\nCheck idrow=True option')
            self.assertEqual(ax.id,'axis','genutil.ASCII: Error axis on variables should be named axis, looks like the axis setup failed!\nCheck axis=True option')
            self.assertTrue((ax[:]==[1.,2.,2.13,4.]).any(), 'genutil.ASCII: Error axis on variables have wrong values, looks like the axis setup failed!\nCheck axis=True option')


           

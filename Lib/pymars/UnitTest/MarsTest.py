import unittest, sys, os, subprocess, genutil.pymars, genutil.mars
from subprocess import Popen, call, PIPE, STDOUT
from genutil.mars import logit as F_logit

sep = os.sep
pardir = os.pardir
logdir = '..'+sep+'logfiles'+sep
baselinedir = '..'+sep+'baselineLogFiles'+sep

def suite():

    s = unittest.TestSuite()
    for methodName in dir(MarsTest):
        if methodName.startswith("test"):
            s.addTest(MarsTest(methodName))        
    return s

class MarsTest(unittest.TestCase):
    #There are 2 tests defined in the class: continuous & categorical
    #The implementation of each test uses a subprocess to perform
    #the test.  The reason is that in python externally generated
    #modules such as mars (generated using f2py) can not be reloaded.
    #What this means if something changes in one execution of mars
    #then it remains changed in the next execution of mars. It was
    #df that brought this problem to light.  Specifically when 
    #testing cross validation df is changed and was picked up in the
    #next call to mars. There is a sparse amount of information on the
    #web regarding this problem, which indicates why it hasn't been
    #solved.  The only solution that seams to work is to start a new
    #process for each execution of mars.
    def __init__(self, testName):
        self.fns = []
        self.cases =[(0,0), (0,2), (1,0), (1,2)]
        unittest.TestCase.__init__(self, testName)
    def SetUp(self):
        pass
    def testContinuous(self):
        #test continuous data
        self.cases = [(0,0)]
        output =[]
        for logisticRegression, crossValidation in self.cases:
            #ARG = 'srun -n2 -p pdebug pyMPI RunMarsAndPymars.py ' 
            ARG = 'python RunMarsAndPymars.py ' 
            ARG = ARG + 'ContinuousTest '
            ARG = ARG + '5 2 1 '
            ids=[('il=',logisticRegression), ('ix=',crossValidation)]
            logfile = filename('pymars_cont',ids) + '.log'
            self.fns = self.fns + [logfile]
            Flogfile = filename('mars_cont',ids) + '.log'
            ARG = ARG + str(logisticRegression) + ' ' + str(crossValidation) + ' '
            ARG = ARG + logdir+logfile + ' ' + logdir+Flogfile + ' '
            ARG = ARG + '[1,1]'
            ARG = ARG.split()
            
            run = Popen(ARG, stderr=PIPE, stdout=PIPE)
            stdout, stderr = run.communicate()

            output += [logfile+'\n', 'stderr=\n', stderr,'stdout=\n',stdout, '\n']
            if len(stderr) > 0:
                print 'stderr=\n', stderr
            print 'stdout=\n',stdout

        #dumpDebug('pymars_cont'+'.debug', output)

    def xtestCategorical(self):
        #test categorical data
        #self.cases = [(0,0)]
        output =[]
        for logisticRegression, crossValidation in self.cases:
            #ARG = 'srun -n3 -p pdebug pyMPI RunMarsAndPymars.py '
            ARG = 'python RunMarsAndPymars.py '
            ARG = ARG + 'CategoricalTest '
            ARG = ARG + '5 1 1 '
            ids=[('il=',logisticRegression), ('ix=',crossValidation)]
            logfile = filename('pymars_cat',ids) + '.log'
            self.fns = self.fns + [logfile]
            Flogfile = filename('mars_cat',ids) + '.log'
            ARG = ARG + str(logisticRegression) + ' ' + str(crossValidation) + ' '
            ARG = ARG + logdir+logfile + ' ' + logdir+Flogfile + ' '
            ARG = ARG + '[1,-1,1]'
            ARG = ARG.split() 
            
            run = Popen(ARG, stderr=PIPE, stdout=PIPE)
            stdout, stderr = run.communicate()
            output += [logfile+'\n', 'stderr=\n', stderr, 'stdout=\n',stdout, '\n']
            if len(stderr) > 0:
                print 'stderr=\n', stderr
            #print 'stdout=\n',stdout

        #dumpDebug('pymars_cat'+'.debug', output)
        
    def xtestMixedContinuousCategorical(self):
        #test categorical data
        self.cases = [(0,0, '[1,-1,-1]'), (0,0,'[1,1,-1]'), (0,2,'[1,-1,-1]'), (0,2,'[1,1,-1]')]
        output =[]
        for logisticRegression, crossValidation, lx in self.cases:
            ARG = 'srun -n4 -p pdebug pyMPI RunMarsAndPymars.py '
            ARG = 'python RunMarsAndPymars.py '
            ARG = ARG + 'MixedContinuousCategorical '
            ARG = ARG + '5 1 1 '
            lxStr = lx[1:-1:] #remove brackets
            lxStr = lxStr.replace (',','') #remove commas
            ids=[('il=',logisticRegression), ('ix=',crossValidation),('lx=', lxStr)]
            logfile = filename('pymars_mixed',ids) + '.log'
            self.fns = self.fns + [logfile]
            Flogfile = filename('mars_mixed',ids) + '.log'
            ARG = ARG + str(logisticRegression) + ' ' + str(crossValidation) + ' '
            ARG = ARG + logdir+logfile + ' ' + logdir+Flogfile + ' '
            ARG = ARG + lx
            ARG = ARG.split() 
            
            run = Popen(ARG, stderr=PIPE, stdout=PIPE)
            stdout, stderr = run.communicate()
            output += [logfile+'\n', 'stderr=\n', stderr, 'stdout=\n',stdout, '\n']
            
            if len(stderr) > 0:
                print 'stderr=\n', stderr
            #print 'stdout=\n',stdout
        #dumpDebug('pymars_mixed'+'.debug', output)
    def xtestMultiCategorical(self):
        #test categorical data
        self.cases = [(0,0, '[-1,-1,-1]'), (0,2,'[-1,-1,-1]')]
        output =[]
        for logisticRegression, crossValidation, lx in self.cases:
            ARG = 'srun -n4 -p pdebug pyMPI RunMarsAndPymars.py '
            ARG = 'python RunMarsAndPymars.py '
            ARG = ARG + 'MultiCategorical '
            ARG = ARG + '5 1 1 '
            lxStr = lx[1:-1:] #remove brackets
            lxStr = lxStr.replace (',','') #remove commas
            ids=[('il=',logisticRegression), ('ix=',crossValidation),('lx=', lxStr)]
            logfile = filename('pymars_multicat',ids) + '.log'
            self.fns = self.fns + [logfile]
            Flogfile = filename('mars_multicat',ids) + '.log'
            ARG = ARG + str(logisticRegression) + ' ' + str(crossValidation) + ' '
            ARG = ARG + logdir+logfile + ' ' + logdir+Flogfile + ' '
            ARG = ARG + lx
            ARG = ARG.split() 
            
            run = Popen(ARG, stderr=PIPE, stdout=PIPE)
            stdout, stderr = run.communicate()
            output += [logfile+'\n', 'stderr=\n', stderr, 'stdout=\n',stdout, '\n']
            
            if len(stderr) > 0:
                print 'stderr=\n', stderr
            #print 'stdout=\n',stdout
        dumpDebug('pymars_multicat'+'.debug', output)
    def xtestDebug(self):   
        logisticRegression, crossValidation = 0,2
        ARG = 'python RunMarsAndPymars.py '
        ARG = ARG + 'ContinuousTest '
        ARG = ARG + '5 2 1 '
        ids=[('il=',logisticRegression), ('ix=',crossValidation)]
        logfile = filename('pymars_cont',ids) + '.log'
        self.fns = self.fns + [logfile]
        Flogfile = filename('mars_cont',ids) + '.log'
        ARG = ARG + str(logisticRegression) + ' ' + str(crossValidation) + ' '
        ARG = ARG + logdir+logfile + ' ' + logdir+Flogfile + ' '
        #ARG = ARG + '[1,-1,1]'
        ARG = ARG.split()
        print ARG
        #run = subprocess.Popen(ARG)
        #rcode = run.wait()
        run = Popen(ARG, stderr=PIPE, stdout=PIPE)
        stderr, stdout = run.communicate()
        #print stderr, stdout
        #print 'returned code=', rcode      
        return
    def tearDown(self):
        #compare the output of pymars with the baseline output 
        for fn in self.fns:
            ok = fileCompare(logdir+fn, baselinedir+fn)
            print 'passed=',ok, fn
        self.fns = []
        sys.stdout.flush()
        return
def filename(root, ids):   
    fn = root
    for id,val in ids:
        fn = fn + '_' + id + str(val)
    return fn
def fileCompare(f1, f2):
    ok = False
    try:
        in1 = open(f1, 'r')
        in2 = open(f2, 'r')  
        ok = (in1.readlines() == in2.readlines())
        in1.close()
        in2.close()
    except:
        print 'file problems with ', f1, f2
    return ok
def dumpDebug(fn, output):
    f=open(fn, 'w')
    for o in output:
        f.write(o)
    f.close()
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MarsTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
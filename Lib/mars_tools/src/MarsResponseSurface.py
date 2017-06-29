from ParseMarsLogFile import *
from LinearSpline import *
from MarsBasisFunction import *
from ParseMarsMemory import *
from types import *
import numpy
from numpy import shape

from genutil.mars import mars as _mars
from genutil.mars import fmod as _fmod

class MarsResponseSurface:
    def __init__(self, arg, SIGN=+1, fm=None, im=None):
        self.const = 0.
        self.functions = []
        self.error = 'none'
        self.inputs = {}
        self.outputs = {}
        self.fm = fm
        self.im = im
        self.SIGN=SIGN

        
        if type(arg) is FloatType:
            self.const = arg
        elif type(arg) is TupleType:
            self.fm, self.im = arg
            mem = ParseMarsMemory(arg)
            BFs = mem.basisFunctions
            self.inputs = mem.inputs
            self.outputs = mem.outputs
            self.createResponseSurface(BFs)   
            #self.createResponseSurfaceFromMarsLogFile(BFs)         
        elif type(arg) is StringType:
            try:
                log = ParseMarsLogFile(arg)
                self.inputs = log.inputParms()
                self.outputs = log.outputParms()
                BFs = log.getBasisFunctions()
                self.createResponseSurfaceFromMarsLogFile(BFs)
            except:
                self.error = 'error reading log file:' + arg
    def __str__(self):
        s = self.error + '\n'
        s = s + 'constant=' + str(self.const) + '\n'
        i = 0
        for i,fnc in enumerate(self.functions):
            s = s + 'function' + str(i) + ' = (' + str(fnc) + ')\n'
        for input in self.inputs:
            val = self.inputs[input]
            s = s + str(input)+ ' = ' + str(val) + '\n'
        for output in self.outputs:
            val = self.outputs[output]
            s = s + str(output)+ ' = ' + str(val) + '\n'
        s = s + 'SIGN=' + str(self.SIGN) + '\n'
        return s
    def __call__(self, xinput):
        if not isinstance(xinput, numpy.ndarray):
            print 'MarsResponseSurface requires numpy (nsamps,p) input array'
            return
        if not isinstance(self.fm, numpy.ndarray):
            print 'MarsResponseSurface requires numpy fm array'
            return
        if not isinstance(xinput, numpy.ndarray):
            print 'MarsResponseSurface requires numpy im array'
            return
        (n,p) = xinput.shape
        if 'dim' in self.inputs.keys() and not p == self.inputs['dim']:
            print 'input array does not have correct number of predictor variables'
            return        
        return _fmod(1,n,p,xinput,self.fm, self.im)
    def __len__(self):
        return len(self.functions)
    def __cmp__(self, x):
        same = (self.const == x.const)
        if same:
            try:
                for s in self.functions:
                    same = same and (x.functions.index(s) >= 0)
                for s in x.functions:
                    same = same and (self.functions.index(s) >= 0)
            except:
                same = False
        compare =int(same)-1
        return compare         
    def addBasisFunction(self,bf):
        self.functions = self.functions + [bf]
    def createResponseSurfaceFromMarsLogFile(self, BFs):
        """ This method does what the name suggests.  Given the basis function
        previously read from a log file, make the linear splines,
        make the basis functions in a computable form and then make the 
        response function in a computable form. """
        
        self.const = BFs[0]
        self.functions = []

        for i in range(1, len(BFs)):
            coef = BFs[i]['coef']
            fnc = MarsBasisFunction(coef)
            parent = i
            while parent != 0:
                #print parent
                curbf = BFs[parent]
                ls = LinearSpline(curbf['var'], curbf['knot'], curbf['sign'])
                fnc.addSpline(ls)
                parent = curbf['parent']
            if coef != 0.0:
                self.addBasisFunction(fnc)
    def createResponseSurface(self, BFs):
        """ This method is a fix to eliminate any basis function
            whose coefficients are zero."""
        
        self.const = BFs[0]
        self.functions = []

        for i in range(1, len(BFs)):
            coef = BFs[i]['coef']
            fnc = MarsBasisFunction(coef)
            parent = i
            while parent != 0:
                #print parent
                curbf = BFs[parent]
                ls = LinearSpline(curbf['var'], curbf['knot'], curbf['sign'])
                fnc.addSpline(ls)
                parent = curbf['parent']
            if coef != 0.0:
                self.addBasisFunction(fnc)
if __name__ == '__main__':   
    from MarsBasisFunction import *
    from LinearSpline import *
    
    ls=LinearSpline(1, .5, -1)
    bf=MarsBasisFunction(2.3)
    bf.addSpline(ls)
    ls1=LinearSpline(2, .5, 1)
    bf.addSpline(ls1)
    rs = MarsResponseSurface(4.5)
    rs.addBasisFunction(bf)
    print str(rs)
    print rs((-1,.75,3))
    mrs = MarsResponseSurface('mars.log', SIGN=-1)
    print mrs
    #log = ParseMarsLogFile('mars.log')
    #BFs = log.getBasisFunctions()
    #mrs.createResponseSurfaceFromMarsLogFile(BFs)
    #print mrs
    print mrs((.5, .5))
    mrs1 = MarsResponseSurface(2.0)
    print (mrs == mrs1)

    
    
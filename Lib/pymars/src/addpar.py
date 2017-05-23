#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy, genutil.pymars

class ADDPAR():
    """ 
    This class implements the Parent Priority Queue as 
    described in Friedman's Fast MARS paper, section 3.
    """
    def __init__(self):
        from genutil.pymars import FLOAT_DTYPE, INT_DTYPE
        self.maxdph = 2000
        self.mpr = 10
        self.mtr = 5
        self.beta = 1.0
        self.ktr, self.last, self.kp, self.itr = 0,0,0,0
        self.que = numpy.zeros(shape=(2+1, self.maxdph+1), dtype = FLOAT_DTYPE)
        self.bfIndex = numpy.zeros(shape=self.maxdph+1, dtype = INT_DTYPE)
        self.m = numpy.zeros(shape=self.maxdph+1, dtype = INT_DTYPE)
        self.familyTree = numpy.zeros(shape=(2+1, self.maxdph+1), dtype = INT_DTYPE)
        return
    def addpar(self,bfIndex):
        #print 'addpar'
        if bfIndex == 0: 
            self.last = 1
            self.que[1,1], self.que[2,1] = numpy.inf, 0.0
            self.m[1] = 1
            self.kp, self.itr, self.bfIndex[1], self.familyTree[[1,2],1] = 0,0,0,0
            self.ktr = int(0.5*(self.mpr-1)+.1)
            return
        
        #keep certain points
        (I,) = numpy.where(self.que[1, 1:self.last+1] >= -0.5)
        I = I+1
        N = len(I)
        self.bfIndex[1:N+1] = self.bfIndex[I]
        self.familyTree[1,1:N+1] = self.familyTree[1,I] 
        self.familyTree[2,1:N+1] = self.familyTree[2,I]
        self.que[1,1:N+1] = self.que[1,I]
        self.que[2,1:N+1] = self.que[2,I]
        self.last = N

        self.last = self.last+1
        if self.last > self.maxdph:
            raise 'increase maxdph in addpar to %10i or larger' %(self.last)
        
        #initialize new entry
        self.que[1,self.last] = numpy.inf
        self.que[2,self.last] = 0.0
        self.bfIndex[self.last] = bfIndex
        self.familyTree[[1,2],self.last] = 0
        self.m[1:self.last+1] = numpy.arange(1,self.last+1)
        
        temp = self.que[1,0:self.last+1]
        temp[0] = -numpy.inf
        #psort(self.sp, self.m, 1, self.last)
        self.m[0:self.last+1] = temp.argsort() 

        J = self.m[1:self.last+1]
        I = numpy.arange(1,self.last+1)
        temp[J] = I + self.beta*(self.itr - self.que[2,J])
        #genutil.pymars.debug.info('itr='+repr((self.itr,self.itr - self.que[2,J])))
        #genutil.pymars.debug.info('temp='+repr(temp[J]))
        self.m[0:self.last+1] = temp.argsort()
        self.kp = max(0, self.last-self.mpr)
        #genutil.pymars.debug.info('kp, mpr='+repr((self.kp, self.mpr)))
        return
    
    def nxtpar(self, parent, child):
        #print 'nxtpar'
        self.kp = self.kp+1
        if self.kp > self.last:
            l = -1
            return -1, child
        parent = self.bfIndex[self.m[self.kp]]
        if self.itr-self.familyTree[2,self.m[self.kp]] > self.mtr or \
            self.itr <= self.ktr:
            self.familyTree[1,self.m[self.kp]] = 0
        child = self.familyTree[1, self.m[self.kp]]
        return parent, child
    def updpar(self, child, val):
        #print 'updpar'
        self.que[1,self.m[self.kp]] = val
        self.que[2,self.m[self.kp]] = self.itr
        if self.familyTree[1,self.m[self.kp]] == 0:
            self.familyTree[1,self.m[self.kp]] = child
            self.familyTree[2,self.m[self.kp]] = self.itr
        return
    def selpar(self, bfIndex):
        #print 'selpar'
        LIST = range(1, self.last+1)
        if len(LIST) > 1:
            LIST.reverse()
        for i in LIST:
        #do 13 i=last,1,-1
            if self.bfIndex[i] == bfIndex: 
                self.familyTree[1,i] = 0
                break 
        return
    def itrpar(self, iarg):
        self.itr = iarg
        return
    def setmpr(self, iarg):
        self.mpr = iarg
        return
    def setbta(self, arg):
        self.beta = arg
        return
    def setfrq(self, arg):
        self.mtr = 1.0/int(max(arg,0.01)+.1)
        return

import numpy
class SETINT():
    def __init__(self):
        from genutil.pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
        from .border import border
        self.mlist = 1001
        self.il = 0
        self.m = numpy.zeros(shape = (2, self.mlist), dtype = INT_DTYPE)
        self.m = border(self.m)
        return
    def setint(self, i,j,k):
        print 'setint'
    	if 0 in [i,j]:
        	self.il = 0
         	return
        if i==j: 
            return
    
        m1 = min(i,j)
        m2 = max(i,j)
        if k == 0:
            #the only way the following loop executes is
            #if there are multiple calls to this method.
            #also there are NO calls to this in the fortran
            #so it appears as though the user must do it in the app
            l = 1
            while l <= self.il:
                if m1 == self.m[1,l] and m2 == self.m[2,l]: 
                    return
                l = l+1
            self.il = self.il+1
            if self.il > self.mlist:
                print " increase parameter mlist in subroutine setint \
                to greater than %i5 and recompile." %(self.il)
                raise 'stop'
            self.m[1,self.il] = m1
            self.m[2,self.il] = m2
            return
       
        ig=0
        l=1
        while m1 != self.m[1,l] or m2 != self.m[2,l]:
            l=l+1
            if l > self.il: 
                break
        else:
            ig=1   
        if ig == 0: 
            return
    
        self.il = self.il-1
        ll=1
        while ll <= self.il:
            ll = ll+1
            self.m[1,ll] = self.m[1,ll+1]
            self.m[2,ll] = self.m[2,ll+1]      
        return
       
    def intlst(self):
        from genutil.pymars import LOG, logger
        #print 'intlst'
        if LOG:
            if self.il > 0:
                logger.info(" interactions prohibited between: ")
                for l in range(1, self.il + 1):
                    logger.info('    var %i  and  var %i '%(self.m[1,l], self.m[2,l]))
        return
    
    def intalw(self, i, j):
        #print 'intalw'
        #this method does get called in the fortran
        #none of the attributes change
        k=1
        m1 = min(i,j)
        m2 = max(i,j)
        l=1
        while m1 != self.m[1,l] or m2 != self.m[2,l]:
            l=l+1
            if l > self.il:
                break 
        else:
            k=0
        return k

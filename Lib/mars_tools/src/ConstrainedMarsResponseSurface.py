from MARStools.MarsResponseSurface import *
import scipy
class ConstrainedMarsResponseSurface(MarsResponseSurface):
    def __init__(self, args, **keywords):
        (initparms, box) = args
        MarsResponseSurface.__init__(self, initparms, **keywords)
        self.box = box
        #print keywords
    def __call__(self, x):
        flag=True
        i=0
        for (a,b) in self.box: 
            if a > x[i] or x[i] > b:
                flag=False
            i=i+1
        
        if flag:        
            return MarsResponseSurface.__call__(self, x)
        else:
            return scipy.Infinity*self.SIGN
                                     
if __name__ == '__main__':   
    box = [(0,1.), (0.,1)]
    cmrs = ConstrainedMarsResponseSurface(('mars.log', box))
    print cmrs
    
    mrs=ConstrainedMarsResponseSurface(('mars.log', box), SIGN=-1)
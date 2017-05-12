#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
#from pymars import LOG, logger
from pymars.border import border
from pymars.nord import nord
#from pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
class NEST():
    def __init__(self):
        from pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
        self.mlist = 201
        self.nlist = 2001
        self.il = 0
        self.jl = 0
        self.m  = numpy.zeros(shape = (4, self.mlist), dtype = INT_DTYPE)
        self.m = border(self.m)
        self.vm = numpy.zeros(shape = self.nlist, dtype = INT_DTYPE)
        self.vm = border(self.vm)
        return
    def nest(self, n,i,j,nv,vals):
        if 0 in [i,j]:
            self.il = 0
            self.jl = self.il
            return
        if i == j:
            return

        print 'nest'    
        ig = 0
        if nv > 0:
            k = 1
            while k <= self.il:
                if m[1,k] in [i, j]:
                    return
                k = k+1
            self.il = self.il+1
            if self.il > self.mlist:
                print "increase parameter mlist in subroutine nest \
                to greater than %i5 and recompile." %(self.il)
                raise 'stop'
            self.m[1,self.il] = i
            self.m[2,self.il] = j
            self.m[3,self.il] = nv
            self.m[4,self.il] = self.jl
            if self.jl+nv > self.nlist:
                print "increase parameter nlist in subroutine nest \
                to greater than %i5 and recompile." %(self.jl+nv)
                stop
            for k in range(1, nv+1):
                self.jl = self.jl+1
                self.vm[self.jl] = vals(k)
            return
      
        k = 1
        while self.m[1,k]!= i or self.m[2,k] != j:
            k = k+1
            if k > self.il:
                break
        else:  
            ig = 1
        if ig == 0: 
            return
    
        self.il = self.il-1
        ll = k
        while ll <= self.il:
            for i in [1,2,3,4]:
                self.m[l,ll] = self.m[l,ll+1]
            ll = ll+1
        return  
    def nstlst(self):
        from pymars import LOG, logger
        if not LOG:
            return
        if self.il == 0:
            return
        print 'nstlst'
        logger.info('\n variable nesting:')
        for k in range(1, self.il+1):
            if self.m[3,k] > 5:
                logger.info( "%3i: var(%3i) exists for var(%i3) = " %(k,self.m[1,k],self.m[2,k]))
                for VM in self.vm[self.m[4,k]+1: self.m[4,k]+self.m[3,k]+1]:
                    logger.info( '%7.1f'%(VM))
            else:
                logger.info( "%3i : var(%3i) exists for var(%3i) = )" %(k, self.m[1,k], self.m[2,k]))
                for VM in self.vm[self.m[4,k]+1: self.m[4,k]+self.m[3,k]+1]:
                    logger.info( '%7.1f'%(VM))
        return      
    def oknest(self, p,lx,cm):
        from pymars import LOG, logger
        if not LOG:
            return
        #print 'oknest'
        l = 1
        while l <= self.il:
            j1 = self.m[1,l]
            jn = self.m[2,l]
            jv = self.m[3,l]
            jp = self.m[4,l]
            if j1 < 1 or j1 > p:
                logger.info(' nesting entry %3i invalid variable %3i to be nested.'%(l,j1))
            if jn < 1 or jn > p: 
                logger.info(' nesting entry %3i invalid nesting variable %3i'%(l,jn))
            if lx[jn] >= 0: 
                logger.info(' nesting entry %3i lx(%3i = %2i  must be < 0.' %(l,jn,lx[jn]))
            k1 = int(cm[2*jn] + .1)
            k2 = int(cm[2*jn+1] + .1)
            
            for k in range(jp+1, jp+jv+1):
                ig = 0
                for kk in range(k1, k2+1):
                    if self.vm[k] == cm[kk]:
                        ig = 1
                        break 
                if ig == 0: 
                    logger.info(' nesting entry %i3, categorical value %12.4g,/,' +  \
                    'not among the data values for variable,%i3,' %(l,self.vm[k],jn))
            l = l+1
        return
    def isnstr(self, j):
        #print 'isnstr'
        jb = 0
        k = 1
        while self.m[2,k] != j:
            if k > self.il: 
                return jb
            k = k+1
            jb = self.m[1,k]
        return jb
    def isfac(self, lm, j, mk, tb, cm):
        #print 'isfac'
        ja = 0
        (ind, ) = numpy.where(j == self.m[1,1:self.il+1])
        if len(ind) == 0:
            return ja
        l = ind[0]+1
        jn = self.m[2,l]
        if cm[2*jn] == 0.0:
            return ja
        jv = self.m[3,l]
        jp = self.m[4,l]
        ig = 0
        ip = lm
        while ip > 0:
            j1 = int(abs(tb[2,ip])+.1)
            if j1 == jn: 
                ig = 1
                break 
            ip = int(tb[4,ip] + .1)
          
        if ig != 0: 
            nc = int(cm[2*jn+1] - cm[2*jn] + 1.1)
            t = tb[2,ip]
            kp = int(tb[3,ip] + .1)
            for l in range(1, nc+1):
                lp = l+kp
                #complicated logic; it must be simpler
                flag1 = (t <= 0)
                flag2 = (cm[lp] == 0.0)
                EXECUTE = (flag1 and flag2) or (not flag1 and not flag2)
                if EXECUTE:
                    ex = cm[int(cm[2*jn]+.1) + l - 1]
                    ig = 0
                    for k in range(jp+1, jp+jv+1):
                        if ex == self.vm[k]:
                            ig = 1
                            break 
                    if ig == 0:
                        ja = -1
                        return ja
            return ja
        ja = l
        norm = nord(lm, tb) + 1
        nc = int(cm[2*jn+1] - cm[2*jn] + 1.1)
        for lk in range(1, mk+1):
            if nord(lk,tb) == norm: 
                jg = 0
                ip = lk
                while ip > 0:
                    t1 = tb[2,ip]
                    j1 = int(abs(t1)+.1)
                    if j1 == jn:
                        kp = int(tb[3,ip] + .1)
                        kg = 0
                        for l in range(1, nc+1):
                            lp = l+kp
                            lon = int(cm[lp] + .1)
                            if t1 < 0.0:
                                if lon == 0:
                                    lon = 1
                                else: 
                                    lon = 0
                            ex = cm[int(cm[2*jn]+.1) + l - 1]
                            ig = 0
                            for k in range(jp+1, jp+jv+1):
                                if ex == vm[k]:
                                    ig = 1
                                    break 
                            if lon == 1 and ig == 0: 
                                kg = 1
                                break 
                            elif lon == 0 and ig == 1: 
                                kg = 1
                                break 
       
                        if kg == 0: 
                            jg = 1
                            break 
                    ip = int(tb[4,ip] + .1)
                if jg != 0:
                    if ieqbf(lk,lm,tb,cm) == 1: 
                        ja = -1
                        return ja
        return ja
    def cmpnst(self, ja, n, x, cm, bl):
        print 'cmpnst'
        jn = self.m[2,ja]
        jv = self.m[3,ja]
        jp = self.m[4,ja]
        for i in range(1, n+1):
            kx = int(x(l,jn) + .1)
            ex = cm[int(cm[2*jn]+.1) + kx-1]
            ig = 0
            for k in range(jp+1, jp+jv+1):
                if ex == self.vm[k]:
                    ig = 1
                    break
            if ig != 1: 
                bl[l] = 0.0
        return bl  
    def getnst(self, ja, cm, j, nv, vals):
        print 'getnst'
        j  = self.m[2,ja]
        jv = self.m[3,ja]
        jp = self.m[4,ja]
        nv = int(cm[2*j+1] - cm[2*j] + 1.1)
        vals[1:nv+1] = 0.0
        #for k in range(1, nv+1):
            #vals[k] = 0.0
        for l in range(jp+1, jp+jv+1):
            k = icat(self.vm[l],j,cm)
            if k > 0: 
                vals[k] = 1.0
        return j, vals

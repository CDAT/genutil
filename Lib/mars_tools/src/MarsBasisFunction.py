class MarsBasisFunction:
    def __doc__(self):
        return """MARS basis function is coefficient times a product
                of linear splines in distinct variables"""  
    def __init__(self, coef):
        self.coef = coef
        self.splines = []
        self.error = ''
    def __str__(self):
        s = 'coef=' + str(self.coef) + '\n'
        i = 0
        for ls in self.splines:
            s = s + 'spline' + str(i) + ' = (' + str(ls) + ')\n'
            i=i+1
        return s
    def __call__(self, x):        
        f = self.coef
        for s in self.splines:
            i = s.var-1
            f = f*s(x[i])
        return f
    def __cmp__(self, x):
        same = (self.coef == x.coef)
        if same:
            try:
                for s in self.splines:
                    same = same and (x.splines.index(s) >= 0)
                for s in x.splines:
                    same = same and (self.splines.index(s) >= 0)
            except:
                same = False
        compare =int(same)-1
        return compare         
    def addSpline(self,ls):
        add = True
        for s in self.splines:
            add = add and (ls.var != s.var)
        if add:
            self.splines = self.splines + [ls]
        else:
            self.error = 'spline has the same variable as one in the list'
        return add

if __name__ == '__main__':    
    from LinearSpline import *       
    ls=LinearSpline(1, .5, -1)
    bf=MarsBasisFunction(2.3)
    bf.addSpline(ls)
    ls1=LinearSpline(2, .5, 1)
    bf.addSpline(ls1)
    bf1 = MarsBasisFunction(2.3)
    bf1.addSpline(ls)
    print bf1.splines
    print bf1.error
    print cmp(bf1, bf)
    print (bf1 == bf)
    print bf.splines.index(ls)
    print str(bf)
    print bf((-1,.75,3))
    print bf.__doc__()

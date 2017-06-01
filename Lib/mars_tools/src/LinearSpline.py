class LinearSpline:
    def __init__(self, var, knot, sign):
        self.var = var
        self.knot = knot
        self.sign = sign
    def __str__(self):
        s = 'var=' + str(self.var) + ' '
        s = s + 'knot=' + str(self.knot) + ' '
        s = s + 'sign=' + str(self.sign) 
        return s
    def __type__(self):
        return 'LinearSpline'
    def __cmp__(self, x):
        compare = 1
        if self.var == x.var and self.knot == x.knot and self.sign == x.sign:
            compare = 0
        if self.var == x.var and self.knot < x.knot:
            compare = -1
        return compare
    def __call__(self, x):
        if self.sign > 0:
            if x <= self.knot:
                y = 0.
            else:
                y = x-self.knot
        else:
            if x >= self.knot:
                y = 0.
            else:
                y = self.knot - x            
        return y
    def dump(self):
        print self.var, self.knot, self.sign

if __name__ == '__main__':    
    ls=LinearSpline(1, .5, -1)
    ls.dump()
    print str(ls)
    print ls(.4)
    ls1=LinearSpline(2, .4, -1)
    print str(ls1)
    print ls1 <= ls
    #rint cmp(ls, ls1)
    ls.dump()
    print ls(.6)
    print ls.var
    print type(ls)
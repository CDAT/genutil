class NextVariable():
    def __init__(self, variable, nextVariable, hockeyStick, compare):
        self.coefficients = []
        self.hockeyStick = hockeyStick
        self.nextVariable = nextVariable
        self.compare = compare
        return
    def decide(self, tb, cm, fvr, txl, tx1, txi, variable, nextVar, index, 
                           newBF, ja, jas, COMPARE, CAT_CHECK= (False, 0, 0, 0, 0)):
        (CHECK, ict, kcp0, kcp, nc) = CAT_CHECK
        newbf = isNewBasisFunction(tb[2:5, index], tb[2:5, 1:index])
        """
        This function does some sort of check of the hockey stick function
        and decides what the next variable is for the current basis function.
        COMPARE is either <= in findBestKnot or < in keepHS
        """
        s = "%e %s %e" %(fvr*tb[1, index], COMPARE, txl)
        if eval(s) and  newbf:
            txl = fvr*tb[1, index]
            tx1 = tb[1, index]
            nextVar = variable
        s = "%e %s %e" %(fvr*tb[1, index], COMPARE, txi)
        if eval(s) and  newbf:
            txi = fvr*tb[1, index]
            newBF[1] = tb[1, index]
            newBF[2:5] = tb[2:5, index]
            jas = ja
            if CHECK and ict != 0:
                cm[kcp0+1:kcp0+nc+1] = cm[kcp+1:kcp+nc+1]
                kcp = kcp0 + nc
                newBF[3] = kcp0
        return tb, cm, newBF, txl, tx1, txi, jas, kcp, nextVar
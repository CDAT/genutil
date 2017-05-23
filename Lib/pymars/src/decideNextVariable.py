import numpy
from genutil.pymars import debug, FLOAT_DTYPE
from .newb import isNewHS
from .border import border
def decideNextVariable(hockeyStick, CM, cm, fvr, txl, tx1, txi, variable, nextVar, 
                       newHS, ja, COMPARE, CAT_CHECK= (False, 0, 0, 0, 0)):
    
    """
    This function does some sort of check of the hockey stick function
    and decides what the next variable is for the current basis function.
    COMPARE is either <= in findBestKnot or < in keepHS
    """
    
    (CHECK, ict, kcp0, kcp, nc) = CAT_CHECK
    jas = -1
    coef = hockeyStick[0]
    #debug.info('hockeyStick='+repr((hockeyStick)))
    #the checks below are very important one in the algorithm. if they
    #pass the new hockey stick might be accepted and the response
    #surface changes. In the original code there was no accounting for 
    #numerical issues. The addition of 1.e-10 was added only after
    #some painful debugging. The effect of adding this is to say 
    #that the new coefficient (b) is taken as the new one if it's
    #really smaller.
    s = "%e %s %e" %(fvr*coef + 1.e-10, COMPARE, txl)
    txlEval = eval(s)
    print 'txlEval=', txlEval
    if eval(s):
        #debug.info('pass1')
        txl = fvr*coef
        tx1 = coef #used to compute the score in update
        nextVar = variable


    #determine if the new hockey stick is best
    s = "%e %s %e" %(fvr*coef + 1.e-10, COMPARE, txi)
    txiEval = eval(s)
    print 'txiEval=', txiEval
    if eval(s):
        txi = fvr*coef
        newHS[1:5] = hockeyStick
        jas = ja
        #adjustments for categorical hockey stick
        if CHECK and ict != 0:
            print 'new Knot=', newHS[2:5]
            cm[kcp0+1:kcp0+nc+1] = CM[1:]#cm[kcp+1:kcp+nc+1]#CM[1:]#
            kcp = kcp0 + nc
            newHS[3] = kcp0
            debug.info('in decideNextVariable:kcp, kcp0, nc ='+repr((kcp, kcp0, nc )))
            debug.info('cm= '+repr((cm[kcp0+1:kcp0+nc+1])))
    return nextVar, newHS, txl, tx1, txi, jas, cm, kcp
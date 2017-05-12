import numpy
from pymars import debug, FLOAT_DTYPE
from pymars.newb import isNewHS
from pymars.border import border
def decideWinner(states, txl, tx1, txi, cm, kcp0, nextVar, newHS):

    """
    This function makes the final decision for the winning hockey stick
    and associated parameters. 
    """

    jas = -1  
    for state in states:    
        if state.makeDecision:
            coef = state.HS[0]
            #debug.info('hockeyStick='+repr((hockeyStick)))
            #the checks below are very important one in the algorithm. if they
            #pass the new hockey stick might be accepted and the response
            #surface changes. In the original code there was no accounting for 
            #numerical issues. The addition of 1.e-10 was added only after
            #some painful debugging. The effect of adding this is to say 
            #that the new coefficient (b) is taken as the new one if it's
            #really smaller.
            
            if state.COMPARE(state.fvr*coef + 1.e-10,  txl):
                txl = state.fvr*coef
                tx1 = coef #used to compute the score in update
                nextVar = state.variable
        
            #determine if the new hockey stick is best
            if state.COMPARE(state.fvr*coef + 1.e-10,  txi):
                txi = state.fvr*coef
                newHS[1:5] = state.HS
                jas = state.ja
                #adjustments for categorical hockey stick
                if state.ict != 0:
                    #cm[kcp0+1:kcp0+nc+1] = cm[kcp+1:kcp+nc+1]
                    #kcp = kcp0 + nc
                    cm[kcp0+1:kcp0+state.nc+1] = state.CM[1:]
                    newHS[3] = kcp0
    return nextVar, newHS, tx1, txi, jas, cm
class hsState():
    def __init__(self, makeDecision, variable, HS, COMPARE, fvr, ja, ict, nc, CM):
        self.makeDecision = makeDecision
        self.variable = variable
        self.HS = HS
        self.COMPARE = COMPARE
        self.fvr = fvr
        self.ja = ja
        self.ict = ict
        self.nc = nc
        self.CM = CM
        return
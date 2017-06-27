import numpy
def mnspan(ms, alf, nep, nnt):
    """ Equations 42-45 in Friedman paper."""
    al2=0.693147
    al25=1.732868
    allf = -numpy.log(1.0-alf)
    fmn = -numpy.log(allf/(nep*nnt))/al25 #equation 43 in natural log form
    fme = -numpy.log(alf*0.125/nep)/al2   #equation 45 in natural log form

    if ms > 0:
        me = int(ms*fme/fmn+0.5)
        mn = ms
    else:
        me = int(fme + 0.5)
        mn = int(fmn + 0.5)
    me = max(me,mn,2)
    nst = nnt-2*me-1
    nnr = int(float(nst)/mn)
    nnl = nst-nnr*mn
    nnr = (nnr+1)*mn-nst
    nst = min(nnl,nnr)
    if nnl <= nnr: 
        nnl = 1
    else:
        nnl = -1
    nnr = int(float(nst)/2)
    mel = me
    me = me+nnl*nnr
    mel = mel+nnl*nnr
    if numpy.fmod(nst,2) != 0:
        mel = mel+nnl
    return mn, me, mel

from genutil.pymars import debug
def holl (jp, cm, t):

    j1 = int(cm[2*jp]+.1)
    j2 = int(cm[2*jp+1]+.1)
    j2 = j2-j1+1
    if j2 > 28:
        h = '   cat. factor > 28 values  '
        return h
    h = '                            '
    j1 = int(float(28-j2)/2)
    j2 = j1+j2-1
    k = int(t+.1)
    #debug.info('t, k, j1, j2 ='+repr((t, k, j1, j2)))
    for j in range(j1, j2+1):
        #debug.info('j ='+repr((j, cm[k+j-j1+1])))
        if cm[k+j-j1+1] > 0.0:
            h = h[:j]+'1'+h[j+1:]
        else:
            h = h[:j]+'0'+h[j+1:]
    #debug.info(h)
    return h

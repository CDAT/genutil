def icat(x, j, cm):
    #print 'icat'
    j0 = int(cm[2*j]+.1) 
    j1 = j0
    j2 = int(cm[2*j+1]+.1) 
    while j2 != j1+1:
        k = int(float(j1+j2)/2)
        if cm[k] == x:
            icat = k-j0+1
            return icat
        if cm[k] < x: 
            j1 = k
        else:
            j2 = k

    if x == cm[j1]:
        icat = j1-j0+1
    elif x == cm[j2]: 
        icat = j2-j0+1
    else:   
        icat = 0
    return icat
    

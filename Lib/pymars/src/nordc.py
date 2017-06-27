def nordc(l, m, tb, cm):
    ip = m
    nordc = 0
    while ip > 0:
        j = int(abs(tb[2,ip])+.1)
        if l == 1: 
            if cm[2*j] == 0.0:
                nordc = nordc+1
        else:    
            if cm[2*j] > 0.0:
                nordc = nordc+1
        ip = int(tb[4,ip]+.1)
    return nordc
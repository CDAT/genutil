#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def lcm(p,nk,tb,cm):
    #my hunch is that lcm means the last relevant entry of the cm array
    #cuurently it is not used
    
    ix = 0
    jj = 0
    for m in range(1, nk+1):
        j = int(abs(tb[2,m])+.1) # added int()
        if cm[2*j] != 0.0: 
            if int(tb[3,m]+.1) > ix: 
                ix = int(tb[3,m]+.1) #added int()
                jj = j
    if ix > 0: 
        lcm = ix + int(cm[2*jj+1]+.1) - int(cm[2*jj]+.1) + 1
    else:
        lcm = 2*p+1
        for j in range(1, p+1):
            if cm[2*j] != 0.0:
                lcm = lcm + int(cm[2*j+1]+.1) - int(cm[2*j]+.1) + 1
    return lcm
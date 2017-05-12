#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def exch(nk, m, k, d, a, b):
    #print 'exch'
    l = k+1
    km1 = k-1
    kp2 = k+2
    r = a[k]
    a[k] = a[l]
    a[l] = r
    r = b[k]
    b[k] = b[l]
    b[l] = r
    for j in [1,2]:
        i = nk+j
        t = d[k,i]
        d[k,i] = d[l,i]
        d[l,i] = t
    t = d[k,k]
    d[k,k] = d[l,l]
    d[l,l] = t
    
    j = 1
    while j <= km1:
        t = d[j,k]
        d[j,k] = d[j,l]
        d[j,l] = t
        j = j+1 
      
    j = kp2
    while j <= m:
        t = d[k,j]
        d[k,j] = d[l,j]
        d[l,j] = t
        j = j+1 
    return d, a, b
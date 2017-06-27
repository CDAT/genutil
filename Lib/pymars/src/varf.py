def varf(nk, d, a, sw, k1, k2):
    #print 'varf'

    s = 0.0
    for i in range(k1, k2+1):
        t = 0.0
        for j in range(k1, k2+1):
            if j <= i:
                u = d[j,i]
            else:
                u = d[i,j]
            t = t + a[j]*u
        s = s + a[i]*t
    varf = s/sw
    return varf

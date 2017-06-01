from numpy import *
def norms(x, y):
    "Computes the L1, L2, Linfinity & RMS difference of two arrays"
    
    n=len(x)    
    if n <= 0:
        return "no data in array"
    
    if len(y) != n:
        return "mismatch of array length"

    #L1 norm
    L1 = 0.0
    L2 = 0.0
    Linfinity = 0.0
    for i in range(n):
         diff = abs(x[i] - y[i])
         L1 = L1 + diff
         L2 = L2 + diff**2
         Linfinity = max(Linfinity, diff)
         
    RMS = sqrt(L2/n)
    L1 = L1/n
    L2 = sqrt(L2/n)    
    return [L1, RMS, Linfinity]
    
def compare(d, request, model):
    
    nvars =d['nvars']-1
    npts = d['npts']
    
    ind=zeros((nvars,npts))
    for i in range(nvars):
        x = d['data'][i]
        for j in range(npts):
            k=0
            go = 1
            while go:
                if request[i][k] <= x[j] and x[j] < request[i][k+1]:
                    go = 0
                else:
                    k = k+1
            ind[i, j] = k
    y=resize(model, (npts, npts))
#    print y
    for i in range(npts):
#        k = ind[0,i] + ind[1,i]*(npts-1)
        print d['data'][2][i], y[ind[0,i], ind[1,i]], ind[0,i], ind[1,i], d['data'][0][i],d['data'][1][i]
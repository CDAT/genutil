#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def jfv(m, tb, jv):

    ip = m
    j = 0
    while ip > 0:
        j = j+1
        jv[j] = int(abs(tb[2,ip])+.1)
        ip = int(tb[4,ip]+.1)

    if j == 1:
        return jv

    j = j-1
    l = 1
    while l != 0:
        l = 0
        for i in range(1, j+1):
            if jv[i] > jv[i+1]: 
                k = jv[i]
                jv[i] = jv[i+1]
                jv[i+1] = k
                l = 1
    return jv

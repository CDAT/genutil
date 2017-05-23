#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def jf(m, j, tb):   
    ip = m
    jf = 0
    while ip > 0:
        jp = int(abs(tb[2,ip])+.1)
        if jp == j:
            jf = 1
        ip = int(tb[4,ip]+.1)
    return jf

#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def nnord(NEST, m, tb):
    ip = m
    nnord = 0
    while ip > 0:
        if NEST.isnstr(int(abs(tb[2,ip])+.1)) == 0:
            nnord = nnord+1
        ip = int(tb[4,ip]+.1)
    return nnord
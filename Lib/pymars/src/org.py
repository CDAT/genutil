#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def org(m1, m2, tb, cm, xs, a):
    k = 0
    for m in range(m1, m2+1):
        k = k+1
        if tb[1,m] == 0.0:
            a[k] = 0.0
        else:
            s = 1.0
            ip = m
            while ip > 0:
                j = int(abs(tb[2,ip])+.1)
                if cm[2*j] == 0.0:
                    s = s*xs[j]
                ip = int(tb[4,ip]+.1)
            a[k] = tb[1,m]/s
    return a
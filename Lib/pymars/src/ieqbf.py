#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def ieqbf(lk, lm, tb, cm):
    print "ieqbf"
    ipo = lm
    lg = 0
    ko = 0
    nc = 0
    while ipo > 0:
        to = tb[2,ipo]
        jo = int(abs(to)+.1)
        jg = 0
        if cm[2*jo] == 0.0: 
            t = tb[3,ipo]
            ic = 0
        else:
            ko = int(tb[3,ipo]+.1)
            nc = int(cm[2*jo+1] - cm[2*jo] + 1.1)
            ic = 1
        ip = lk
        while ip > 0:
            t1 = tb[2,ip]
            j1 = int(abs(t1)+.1)
            if j1 == jo: 
                if ic == 0:
                    if to*t1 > 0.0:
                        if ieq(t,tb[3,ip],1.0) == 1: 
                            jg = 1
                            break 
                else:
                    kp = int(tb[3,ip]+.1)
                    kg = 0
                    for l in range(1, nc+1):
                        lo = l+ko
                        lp = l+kp
                        lon = int(cm[lo]+.1)
                        lop = int(cm[lp]+.1)
                        if to < 0.0:
                            if lon == 0:
                                lon = 1
                            else:
                                lon = 0
                        if t1 < 0.0: 
                            if lop == 0: 
                                lop = 1
                            else:
                                lop = 0
                        if lon != lop:
                            kg = 1
                            break 
                    if kg == 0: 
                        jg = 1
                        break 
                    else:
                        ip = int(tb[4,ip]+.1)
        if jg == 0:
            lg = 1
            break 
        else:
            ipo = int(tb[4,ipo]+.1)

    if lg == 0: 
        ieqbf = 1
    else:
        ieqbf = 0
    return ieqbf
#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
def ibfext(m, tb, cm):
    print 'ibfext'
    ibfext = 0
    norm = nord(m,tb)
    for l in range(1, m):
        if nord(l,tb) == norm: 
            if ieqbf(l, m, tb, cm) != 0: 
                ibfext = 1
                return ibfext
    return ibfext
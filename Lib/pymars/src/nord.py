#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
def nord(m, tb):
    #print 'nord'
    #x = numpy.array(tb[4,m:]+.1)
    #(j,) = numpy.where(x <= 0)
    #nord = 1 + len(j)
    ip = m
    nord = 0
    while ip > 0:
        nord = nord+1
        ip = int(tb[4,ip]+.1)

    return nord
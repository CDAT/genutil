import numpy
def init(w, y):
    sw = w[1:].sum()
    wn = numpy.inner(w[1:], w[1:])
    yb = numpy.inner(w[1:], y[1:])
    yb = yb/sw
    wn = sw**2/wn
    s = numpy.inner(w[1:], (y[1:] - yb)**2)
    return sw, wn, yb, s
# Emulate array_indexing

from . import array_indexing
import numpy


def extract(a, b):
    # -1 means all missing so 0 is fine
    b = numpy.where(numpy.equal(b, -1), 0, b)
    return b.choose(a)


rank = array_indexing.rank
set = array_indexing.set

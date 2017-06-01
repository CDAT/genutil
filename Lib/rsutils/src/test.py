from VVutil.Geometry import *
from mars import *
import numpy, sets

x = makeTupleSet([[1, 2, 3, 0],[4, 5, 6,4]])
y = makeTupleSet([[6, 7],[5, 6]])
print x.union(y)
print x.intersection(y)
print project(x.union(y), 0)
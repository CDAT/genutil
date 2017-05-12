#Copyright (c) 2010, LLNS, LLC.  See "Copyright" for full copyright notice.
import numpy
from pymars import ARRAY_SIZE, FLOAT_DTYPE, INT_DTYPE
class CPTB:
	""" make a copy of the tb array"""
	def __init__(self, nk):
		self.ub = numpy.zeros(shape=(5+1, nk+1), dtype=FLOAT_DTYPE)
		return
	def cptb(self, nk, tb):
		self.ub[1:,1:] = tb[1:,1:]
		return
	def setz(self, l):
	    self.ub[1,l] = 0.0
	    return

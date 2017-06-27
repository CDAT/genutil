from VVutil.Metrics import *
from Numeric import *
from VVutil.ColumnIO import *

d=initDict()
d['npts']=3
d['nvars']=3
d['data']=[[1.5, 2.4, 2.7], [4.5, 5.5, 4.8], [7., 8., 9.]]

request=array([[1.,2.,3.], [4., 5., 6.]])

model = array([-5, -4, -3, -2., -1., 0.0, 1.0, 2., 3.])

#print d
#print request
#print model

compare(d, request, model)
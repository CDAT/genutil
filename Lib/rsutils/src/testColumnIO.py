import os
from VVutil.ColumnIO import *

# Read a simple column IO data file
simplePath = 'C:' + os.sep + 'projects' + os.sep + 'workspace' + os.sep + 'MarsController' + os.sep + 'marsctl' + os.sep + 'MarsTrnDatTest002.tpl'
             
simpleData = readData(simplePath)

print "simple data read"

print simpleData['filename']
print simpleData['nvars']
print simpleData['npts']
print simpleData['vars']
print simpleData['data']


# Read an 'omar' column IO data file
omarPath = 'C:' + os.sep + 'projects' + os.sep + 'workspace' + os.sep + 'VVutil' + os.sep + 'omar.dat'

omarData = readData(omarPath, "omar")

print "omar data read"

print omarData['filename']
print omarData['nvars']
print omarData['npts']
print omarData['vars']
print omarData['data']

# Read a 'ptab' generated data file

ptabPath = 'C:' + os.sep + 'projects' + os.sep + 'workspace' + os.sep + 'ResponseSurfaceAnalysisToolkit' + os.sep + 'atp001a_0.tab'
ptabData = readData(ptabPath, "ptab")

print "ptab data read"

print ptabData['filename']
print ptabData['nvars']
print ptabData['npts']
print ptabData['vars']
print ptabData['data']
from string import *
import re
import numpy
from types import *

def readData(fn, headerType=None):
    "read data from file fn, return a dictionary containing the keys:\r\n\
    vars, nvars, filename, npts, and data.  vars is a list of strings\r\n\
    data is a list of numeric arrays that corresponds to the vars. \r\n\
    Example:  d = readData('indata')"
    
    data = {}
    d=[]
    
    p_ptab_header = re.compile('#') #will match a '#'
    p_integer = re.compile('[-+]?\d+')
    p_float = re.compile('[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?')

    f = open(fn, 'r')
    
    for i, line in enumerate(f):
        s=line.strip().split()
        if i == 0:
            #Process header line
            #Special processing for special headers
            if headerType == 'omar':
                s[0] = "edge#"
                s.append("inters1d")
                s.append("inters0d")
        
            if headerType == "ptab":
                #find first element with non-alphanumeric character
                #remove all elements before that element, and that element itself
                for j, elem in enumerate(s):
                    if p_ptab_header.match(s[j]) != None:
                        for k in range(j+1):
                            s.remove(s[0])
                        break

                #remove trailing comment in header line
                #find first element with non-alphanumeric character
                #remove that element itself, and all elements after it
                for j, elem in enumerate(s):        
                    if p_ptab_header.match(s[j]) != None:
                        for k in range(j-1,len(s)):
                            s.remove(s[-1])
                            
            data['vars']= s
            nvars = len(s)
            
            for j in range(nvars):
                d.append([])

            data['nvars']=nvars
            data['filename']=fn
        
        else:
            for j, elem in enumerate(s):
                if '.' not in elem:
                    if headerType == "ptab":
                        try:
                            felem = float(elem)
                            d[j].append(felem)
                        except TypeError:
                            print "VVUtil.ColumnIO: error converting %s to float" % elem
                            try:
                                ielem = int(elem)
                                d[j].append(ielem)
                            except TypeError:
                                print "VVUtil.ColumnIO: error converting %s to int" % elem
                    elif headerType == "omar":
                        try:
                            ielem = int(elem)
                            d[j].append(ielem)
                        except TypeError:                               
                            try:
                                felem = float(elem)
                                d[j].append(felem)
                            except TypeError:
                                print "VVUtil.ColumnIO: error converting %s to float" % elem
                    else:
                        try:
                            felem = float(elem)
                            d[j].append(felem)
                        except TypeError:
                            print "VVUtil.ColumnIO: error converting %s to float" % elem
                            try:
                                ielem = int(elem)
                                d[j].append(ielem)
                            except TypeError:
                                print "VVUtil.ColumnIO: error converting %s to int" % elem
                else:
                    try:
                        felem = float(elem)
                        d[j].append(felem)
                    except TypeError:
                        print "VVUtil.ColumnIO: error converting %s to float" % elem
            #print "d[%d]: " % i
            #print d[i]               
    f.close()
    data['npts']=i
    
    din = []
    for i in range(nvars):
        typelist = [type(elem) for elem in d[i]]
        if FloatType in typelist:
            din.append(numpy.array(d[i], dtype=numpy.float64))
        else:
            din.append(numpy.array(d[i], dtype=numpy.int32))
        
    data['data']= din
    
    return data

def writeData(fn, dout):
    "write data from a dictionary to file fn in column format. The keys:\r\n\
    of the dictionary are: vars, nvars, filename, npts, and data, as described\r\n\
    in readData. Example: writeData('outdata', dout)\r\n"
    
    f=open(fn, 'w')
    cr = "\n"
    blank = "  "
    
    s= ''
    for i in range(dout['nvars']):
        s=s+dout['vars'][i] + blank
    s=s+cr
    f.write(s)
        
    for i in range(dout['npts']):
        s=''
        for j in range(dout['nvars']):
            s = s + str(dout['data'][j][i])+ blank
        s=s+cr
        f.write(s)
             
    f.close()
    return
def initDict():
    "initialize a default dictionary."
    
    d={}
    d['vars']= []
    d['nvars']=0
    d['filename']=" "
    d['npts']=0
    d['data']= []
    
    return d
            
           

from VVutil.ParseMarsLogFile import *
import string
def findLine(lines, searchString):
    nlines = len(lines)
    i = -1
    line = ''
    found = False
    while not found and i < nlines:
        found = (string.find(line, searchString) >= 0)
        if not found:
            i =i+1
            if i != nlines:
                line = lines[i]
    return (found, i)
def getBasisFunctions(file):
   
    bf_parms = []
    coefs = []

    (found, loc) = findLine(file, 'basfn(s)')
    if not found:
        return ['string not found'] 

    bf_lines = []
    next = loc+1
    line = file[next]
    while string.strip(line) != '':
        bf_lines = bf_lines + [line]
        next = next + 1
        line = file[next]
            
    for line in bf_lines:
        vals = string.split(line)
        bf = []
        for val in vals:
            if string.find(val, '.') >= 0:
                val = float(val)
            else:
                val = int(val)    
            bf = bf + [val]
        if type(bf[0]) == type(1) and type(bf[1]) == type(1):
            bf_parms = bf_parms + [bf[1:]]
            bf.remove(bf[1])
        bf_parms = bf_parms + [bf]

    testString = 'coef'
    coef_lines = []
    for line in file[next+1:]:
        loc = string.find(line, testString)
        if loc >= 0:
            coef_lines = coef_lines + [line]
            
    for line in coef_lines:
        c = string.split(line)
        for v in c[1:]:
            coefs = coefs + [float(v)]
                       
    BFs = [coefs[0]]
    nBFs = len(coefs)
    for i in range(1,nBFs):
        coef = coefs[i]
        bf_parm = bf_parms[i]
        var = int(bf_parm[-3])
        knot = bf_parm[-2]
        parent = int(bf_parm[-1])
        bf = {'var':var, 'knot':knot, 'coef':coef, 'parent': parent}
        BFs = BFs + [bf]
        
    return BFs
    
file = readLogFile('mars.log')
i = findLine(file, 'basfn(s)')
for bf in getBasisFunctions(file):
    print bf
import string
class ParseMarsLogFile:
    def __init__(self, fn):
        self.fn = fn
        self.file = []
        f=open(fn, 'r')

        line = f.readline()
        for line in f.readlines():
            self.file.append(line)
        f.close()
    def nBasisFunctionsUsed(self):    
        testString = "anova decomposition on"
        nBFU = -1
        for line in self.file:
            loc = string.find(line, testString)
            if loc >= 0:
                loc = loc+len(testString)
                nBFU = string.atoi(line[loc:loc+3])
                break
        return nBFU
    def inputParms(self):
        testString = "input parameters"
        for line in self.file:
            loc = string.find(line, testString)
            if loc >= 0:
                i = self.file.index(line)
                names = self.file[i+1]
                names = string.split(names)
                values  = self.file[i+2]
                values = string.split(values)
            
                inputs = {}
                for i in range(len(names)):
                    val = values[i]
                    if string.find(val, '.') >= 0:
                        val = float(val)
                    else:
                        val = int(val)                    
                    inputs[names[i]] = val
        return inputs
    def outputParms(self):
        testString = "piecewise cubic fit on"
        for line in self.file:
            loc = string.find(line, testString)
            if loc >= 0:            
                outputs = {}
                i = line.find('n')+1
                j = line.find('basis')
                outputs['nBFused'] = int(line[i:j])
                
                i = line.find('=')+1
                outputs['gcv'] = float(line[i:])
        return outputs
    def findLine(self, lines, searchString):
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
    def getBasisFunctions(self):
   
        bf_parms = []
        coefs = []
        
        #get the list of mins and maxs from output
        (found, loc) = self.findLine(self.file, 'var ')
        next = loc+1
        line = self.file[next]
        minsmaxs = []
        while string.strip(line) != '':
            vals = string.split(line)
            minsmaxs = minsmaxs + [[float(vals[1]), float(vals[-1])]]
            next = next+1
            line = self.file[next]
        #get a list of basis functions
        (found, loc) = self.findLine(self.file, 'basfn(s)')
        if not found:
            return ['string not found'] 

        bf_lines = []
        next = loc+1
        line = self.file[next]
        while string.strip(line) != '':
            bf_lines = bf_lines + [line]
            next = next + 1
            line = self.file[next]
            
        for line in bf_lines:
            vals = string.split(line)
            bf = []
            for val in vals:
                if string.find(val, '.') >= 0:
                    val = float(val)
                elif string.find(val, '*') >= 0:
                    val = 1.e30
                else:
                    val = int(val)    
                bf = bf + [val]
            #add an integer that will be changed to a sign(1 or -1)
            bf =bf + [0]
            if type(bf[0]) == type(1) and type(bf[1]) == type(1):
                bf[-1]= 1
                bf_parms = bf_parms + [bf[1:]]
                bf[-1]=-1
                bf.remove(bf[1])
            else:
                var = int(bf[-4])
                knot = float(bf[-3])
                minmax = minsmaxs[var-1]
                midpnt = (minmax[0]+minmax[1])/2.
                if knot < midpnt:
                    bf[-1]=1
                else:
                    bf[-1]=-1
            bf_parms = bf_parms + [bf]

        #get the coefficients
        testString = 'coef'
        coef_lines = []
        for line in self.file[next+1:]:
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
            var = int(bf_parm[-4])
            knot = float(bf_parm[-3])
            parent = int(bf_parm[-2])
            sign = bf_parm[-1]
            bf = {'var':var, 'knot':knot, 'coef':coef, 'parent': parent, 'sign':sign}
            BFs = BFs + [bf]
        
        return BFs
if __name__ == '__main__':  
    log = ParseMarsLogFile('mars.err')
    for bf in log.getBasisFunctions():
        print bf
    b=log.getBasisFunctions()[1]
    print type(b['var'])
    print type(b['sign'])
    print type(b['parent'])
    print type(b['knot'])
    print type(b['coef'])
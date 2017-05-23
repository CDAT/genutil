import numpy
class ResponseSurfaceOutput():
    def __init__(self, fn):
        self.file= open(fn, 'w')       
    def output(self, header, az, tb):
        """dump the response surface to a file"""
    
        lines = [header+'\n']
        i,nk=tb.shape
        (i,) = numpy.where(tb[1,1:] != 0.0)
        lines = lines + [str(1+len(i))+'\n']
        lines = lines + [str(1) + ' ' + str(az)+'\n']
        for i in range(1,nk):
            coef = tb[1,i]
            cnt = 1
            out = str(coef)+' '
            if coef != 0.0:
                j=i
                while j>0:
                    cnt+=3
                    var = int(tb[2,j])
                    knot = tb[3,j]
                    orient = 1.0
                    if var < 0:
                        orient = -1.0
                    j = int(tb[4,j])
                    out = out + str(abs(var)) + ' ' + str(knot) + ' ' + str(orient) + ' '
                lines = lines + [str(cnt) + ' ' + out + '\n']
        lines = lines + ['\n']
        self.file.writelines(lines)
        return   
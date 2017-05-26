from RunMarsAndPymars import *
from ioTools import *

logisticRegression, crossValidation = 0,0
logdir = './'
#ARG = 'python RunMarsAndPymars.py ' 
ARG = ''
ARG = ARG + 'ContinuousTest '
ARG = ARG + '5 2 1 '
ids=[('il=',logisticRegression), ('ix=',crossValidation)]
logfile = 'junk_pymars_cont.log'
#self.fns = self.fns + [logfile]
Flogfile = 'junk_mars_cont.log'
ARG = ARG + str(logisticRegression) + ' ' + str(crossValidation) + ' '
ARG = ARG + logdir+logfile + ' ' + logdir+Flogfile + ' '
ARG = ARG + '[1,1]'
ARG = ARG.split()

marsParameters()
runMarsAndPymars = RunMarsAndPymars()
args = ARG
print args
testName = args[0]
input, output = getData(testName)
nk = int(args[1])
mi = int(args[2])
m  = int(args[3])
logisticRegression = int(args[4])
crossValidation = int(args[5])
logfile = args[6]
Flogfile = args[7]
lx = eval(args[8])  #do this to convert string to a list

rMaPy = RunMarsAndPymars()
rMaPy.execute(testName, input, output, nk, mi, m, lx,
                         logisticRegression, crossValidation,
                         logfile, Flogfile)
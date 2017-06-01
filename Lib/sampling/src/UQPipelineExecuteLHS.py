#!/usr/local/bin/python

# ---------------------------------------------------------------------------
#                           University of California
#                                Copyright 2007
#
# UQPipelineExecuteLHS - Interface used by the UQ Pipeline Frameworks for executing LHS.
#
# Defined Interfaces:
#  ParameterIs(self,lowvalue,highvalue): Define a new parameter is sample.  
#    Provide low value and high value.
#    EXCEPTION Thrown: InvalidParameter
#
#   
#  Modifier            Date     Description
#  --------          --------  -----------
#  David Domyancic   04-12-07  Initial version
#  David Domyancic   07-23-07  Moved to the UQpipeline project. Modified to 
#                               read the lhs.in file to read in the number of
#                               sample points/parameters.
#
# ---------------------------------------------------------------------------

import string 

import LHS

class UQPipelineExecuteLHS:
    
    def __init__(self):
        
        self.__Parameters = []
        self.__NumberOfPoints = 0
    
    #Execptions defined for 'ParameterIs'
    InvalidParameter = "lowvalue is greater than highvalue"
    
    def ParameterIs(self,lowvalue,highvalue):
        
        if highvalue < lowvalue:
           raise InvalidParameter
        
        NewParameter = (lowvalue,highvalue)
        self.__Parameters.append(NewParameter)
    
    #Execptions defined for 'NumberOfPointsIs'    
    InvalidNumberOfPoints = "NumberOfPoints cannot be less than 0"
    
    def NumberOfPointsIs(self,NumberOfPoints):
        
        if NumberOfPoints < 0:
            raise InvalidNumberOfPoints
        
        self.__NumberOfPoints = NumberOfPoints
    #end 'NumberOfPointsIs'
    
    def NumberOfPoints(self):
        return self.__NumberOfPoints
    #end 'ParameterIs'
    
    ExecuteFailed = "Execution of UQPipelineExecuteLHS failed"
    def Execute(self):
        
        r = open("LHS.in",'r')
        
        #read in the number of sample points
        line = r.readline()
        line = line[:-1]
        self.NumberOfPointsIs(int(line))
        
        #skip this line: posesses names of variables
        line = r.readline()  
        
        #read in the parameters
        for line in r:
                        
            (lowerboundary,upperboundary) = line.split()
            self.ParameterIs(float(lowerboundary), float(upperboundary))
        #end for
            
        listofvalues = LHS.genUniformLHSdata(self.__NumberOfPoints,self.__Parameters)
        
        r.close()
        
        f = open("LHS.out", 'w')
        
        NumberOfParameters = len(listofvalues)
                
        for i in range(self.__NumberOfPoints):
            
            s = ""
            for j in range(NumberOfParameters):
            
                s += str(listofvalues[j][i]) + " "
            
            s += "\n"
            f.write(s)
        
        f.close()
    #end 'Execute'

if __name__ == "__main__":
    
    lhs = UQPipelineExecuteLHS()
    #lhs.NumberOfPointsIs(1000)
    #lhs.ParameterIs(0.0,1.5)
    #lhs.ParameterIs(0.0,1.5)
    lhs.Execute()

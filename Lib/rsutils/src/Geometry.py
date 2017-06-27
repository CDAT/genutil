#make a set of tuples from a list of lists
def makeTupleSet(input):
    import sets
    output = sets.Set()
    
    ncoord=len(input)
    npts=len(input[0])
    
    for i in range(npts):
        tup=()
        for j in range(ncoord):
            tup=tup + (input[j][i],)
        output.add(tup)
    return output
#compute the projection of a set of tuples onto the ith coordinate
def project(data, coord):
    import numpy
    
    output = []
    npts = len(data)
    for x in data:
        output = output + [x[coord]]
  
    return numpy.array(output)  
        
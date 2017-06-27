def genUniform1d(npts, box):  
    (low, high) = box[0]
    delta = (float(high) - low)/(npts-1)
    pnts = []
    for i in range(npts):
        pnt = i*delta + low
        pnts = pnts + [pnt]          
    return [pnts]  
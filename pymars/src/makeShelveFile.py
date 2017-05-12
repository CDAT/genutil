def makeShelveFile(fn, localdict, data):
    import shelve
    
    db = shelve.open(fn)

    for d in data:
        for k in localdict.keys():
            if id(localdict[k]) == id(d):
                db[k] = d
                break
    db.close()
if __name__ == '__main__':
    import numpy, shelve
    x = 2
    y = numpy.array([1,2,3])
    fn = 'shtest'
    makeShelveFile(fn, locals(), [x,y])
    
    db = shelve.open(fn)
    print db
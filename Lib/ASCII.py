# Adapted for numpy/ma/cdms2 by convertcdms.py
import MV2,numpy,cdms2,genutil

def get_parenthesis_content(code):
  opened = 0
  closed = 0
  for i,s in enumerate(code):
    if s=="(":
      opened+=1
    elif s==")":
      closed+=1
    if opened!=0 and closed==opened:
      return code[code.find("(")+1:i]
  return ""

def make_var(lap,id=None,shape=None):
    lap=MV2.array(lap)
    if shape is not None:
        lap=MV2.reshape(lap,shape)
    if id is not None:
        lap.id=id
    return lap

def readAscii( text_file ,header=0, ids=None, shape=None, next='------',separators=[';',',',':']):
    """
    Reads data from an ascii file to generate a list of transient(s)/varable(s)

    :Example:

        .. doctest:: genutil_ASCII_readascii

            >>> vars=genutil.ASCII.readAscii("vars.txt") # use default params

    :param text_file: A string, containing the path to an ASCII File to read from.
    :type text_file: str

    :param header: Number of header lines, these lines will be skipped.
    :type header: int

    :param ids: List of values to use as variable ids (1 per variable returned)
    :type ids: list

    :param shape: use the tuple/list in this list as the final shape of the variable read.
    :type shape: tuple or list

    :param next: character string marking separation between variables (i.e. '------')
    :type next: str

    :param separators: List of characters recognized as column separators (i.e. [';',',',':'])
    :type separators: list

    :returns: List containing transient(s) variable(s) possibly named after ids and reshaped from the 'shape' option.
    :rtype: list
    """
    sep=[]
    if isinstance(separators,str):
        separators=separators.split()
    for s in separators:
        sep.append(s)
        
    f=open( text_file )
    lst = f.readlines( )
    f.close( )
    if not isinstance(ids,(tuple,list)):
        ids=[ids]
    vars=[]
    lap=[]
    for l in lst[header:]:
        for s in sep:
            l=l.replace(s,' ')
        for s in l.split():
            if s==next:
                if len(vars)>len(ids)-1:
                    Id=None
                else:
                    Id=ids[len(vars)]
                if shape is None or len(vars)>len(shape):
                    sh = None
                else:
                    sh =shape[len(vars)]
                if lap!=[]:
                    vars.append(make_var(lap,shape=sh,id=Id))
                    lap=[]
            else:
                if s!='':
                    lap.append(float(s))
    if len(vars)>len(ids)-1:
        Id=None
    else:
        Id=ids[len(vars)]
    if lap!=[]:
        vars.append(make_var(lap,shape=shape,id=Id))
    if len(vars)>1:
        return vars
    else:
        return vars[0]


def read_col( text_file ,header=0, cskip=0, cskip_type='columns', axis=False, ids=None, idrow=0, separators=[';',',', ':']):
    """
    Reads column-stored data from ASCII files

    :Example:

        .. doctest:: genutil_ASCII_read_col

            >>> vars = genutil.ASCII.read_col("vars.txt") # use default params


    :param text_file: ASCII File to read from.
    :type text_file:

    :param header: Number of header lines, these lines will be skipped.
    :type header: int

    :param cskip: Number of 'column'/'character' to skip (dummy column)
    :type cskip: int

    :param cskip_type: One of 'columns' or 'characters'. Specifies which should be skipped.
    :type cskip_type: str

    :param axis: Boolean flag indicating whether to use as the values for the first column as
        variable axis (x values in y(x)).
    :type axis: bool

    :param idrow: Is the first row representing the ids of var generated.
    :type idrow:

    :param ids: (None) use the values in this list as variable ids (1 per column returned)
    :type ids:

    :param separators: ([';',',', ':']) List of character recognized as column separator
    :type separators:

    :returns: List containing 1 transient variable per column in the files.
                Variable ids are optionaly determined by first row.
                Variable axis may be the first column.
    :rtype: list
    """

    sep=[]
    if isinstance(separators,str):
        separators=separators.split()
    for s in separators:
        sep.append(s)
        
    f=open( text_file )
    lst = f.readlines( )
    f.close( )
    lst=lst[header:]
    if not isinstance(ids,(tuple,list)):
        ids=[ids]        
    vars=None
    for l in lst:
        if cskip_type=='characters':
            l=l[cskip:]
        for s in sep:
            l=l.replace(s,' ')
        sp=l.split()
        if cskip_type=='columns':
            sp=sp[cskip:]
        if vars is None:
            nvars=len(sp)
            vars=[]
            for i in range(nvars):
                vars.append([])
        if idrow:
            ids=sp
            idrow=0
        else:
            for i in range(nvars):
                vars[i].append(float(sp[i]))
    for i in range(nvars):
        vars[i]=MV2.array(vars[i])
        if ids!=[None]:
            vars[i].id=ids[i]
    if axis:
        id=vars[0].id
        axval=vars.pop(0).filled()
        axindices=numpy.argsort(axval)
        ax=cdms2.createAxis(numpy.sort(axval))
        ax.id=id
        for i in range(nvars-1):
            tmp=MV2.array(genutil.arrayindexing.get(vars[i],axindices))
            tmp.id=vars[i].id
            tmp.setAxis(0,ax)
            vars[i]=tmp
    if len(vars)>1:
        return vars
    else:
        return vars[0]
    
readAsciiCols = read_col

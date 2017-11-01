# Adapted for numpy/ma/cdms2 by convertcdms.py
import numpy
from . import statistics
import numpy.ma
import cdms2
import genutil


def get(Array, Indices, axis=0):
    """
    Arrayrrayindexing returns Array[Indices], indices are taken along dimension given with axis

    :Example:

        .. doctest:: arrayindexing_get

            >>> import numpy as np
            >>> Array=np.array([2,3,1,0,1,2,3])
            >>> Indices=[0,-1, len(Array)-2] # get the first, last, and second-to-last indices of the array
            >>> C=get(Array,Indices,axis=0) # i.e. C=Array[Indices]

    :param Array: A cdms2 variable, or numpy array, to access the indices of
    :type Array: cdms.tvariable.TransientVariable or numpy.array

    :param Indices: List of integers specifying the indices of Array to access and return

        .. note::

            Negative index value will access indices starting from the end of the array.
            i.e. -1 will be the last item.
    :type Indices: list

    :param axis: Axis of a cdms variable
    :type axis: int or str
    """
    # First some checks

    isma = numpy.ma.isMA(Array)
    if isinstance(Indices, int):
        return Array[Indices]
    if Indices.dtype not in [numpy.int, numpy.int32, numpy.int16]:
        raise "Error indices array must be made of integers (try: Indices=Indices.astype('l') first)"

    if cdms2.isVariable(Array):
        xatt = Array.attributes
        id = Array.id

    if len(Array.shape) != len(Indices.shape):
        Array, Indices, weights, axis, ax = statistics.__checker(
            Array, Indices, None, axis, smally=1)
        if isinstance(Indices, int):
            return Array[Indices]
        if Indices.shape != Array.shape[1:]:
            raise "Error incompatible shapes: " + \
                str(Array.shape) + " and " + str(Indices.shape)
    else:
        Array, Indices, weights, axis, ax = statistics.__checker(
            Array, Indices, None, axis)
        if Indices.shape != Array.shape:
            raise "Error incompatible shapes: " + \
                str(Array.shape) + " and " + str(Indices.shape)

    m = Array.mask
    if not isinstance(Indices, int):
        # Sometihng happened with masking of y by x mask
        Indices = Indices.data.astype('i')
    # print Array.data.dtype.char,Indices.dtype.char
    C = genutil.array_indexing.extract(Array.data, Indices)
    if m is not numpy.ma.nomask:
        M = genutil.array_indexing.extract(m.astype('i'), Indices)
        C = numpy.ma.masked_where(M, C, copy=0)
    elif isma:
        C = numpy.ma.array(C, copy=0, mask=None)
    if ax is not None:
        C = cdms2.createVariable(C, axes=ax, id=id, copy=0)
        for at in list(xatt.keys()):
            setattr(C, at, xatt[at])
    return C


def set(Array, Indices, Values, axis=0):
    """
    Arrayrrayindexing set Array[Indices] with Values, indices are taken along dimension given with axis

    :Example:

        .. doctest:: arrayindexing_set

            >>> import numpy as np
            >>> Array=np.array([2,3,1,0,1,2,3])
            >>> Indices=[0,-1, len(Array)-2] # get the first, last, and second-to-last indices of the array
            >>> Values = [5, 7, 9]
            >>> Array=set(Array,Indices,Values,axis=0) # i.e. Array[Indices]=Values


    :param Array: A cdms2 variable, or numpy array, to set the indices of
    :type Array: cdms.tvariable.TransientVariable or numpy.array

    :param Indices: List of integers specifying the indices of Array to access and set.

        .. note::

            Negative index value will access indices starting from the end of the array.
            i.e. -1 will be the last item.
    :type Indices: list

    :param axis: Axis of a cdms variable
    :type axis: int or str
    """
# if Indices.ndim==0:
# Array[Indices]=Values
    # First some checks
    # isma=numpy.ma.isMA(Array)
    if Indices.dtype not in [numpy.int, numpy.int32, numpy.int16]:
        raise "Error indices array must be made of integers (try: Indices=Indices.astype('l') first)"

    if cdms2.isVariable(Array):
        xatt = Array.attributes
        id = Array.id
    if len(Array.shape) != len(Indices.shape):
        crap, Indices, crap, axis, ax = statistics.__checker(
            Array, Indices, None, axis, smally=1)
        Array, Values, crap, axis, ax = statistics.__checker(
            Array, Values, None, axis, smally=1)
        if Indices.shape != Array.shape[1:]:
            raise "Error uncompatible shapes: " + \
                str(Array.shape) + " and " + str(Indices.shape)
    else:
        Array, Indices, Values, axis, ax = statistics.__checker(
            Array, Indices, Values, axis)
        if Indices.shape != Array.shape:
            raise "Error uncompatible shapes: " + \
                str(Array.shape) + " and " + str(Indices.shape)

    m = numpy.ma.getmask(Array)
    mv = numpy.ma.getmask(Values)
    if Indices.ndim > 0:
        Indices = Indices.data  # Something happened with masking of y by x mask
        Values = Values.data
    genutil.array_indexing_emulate.set(Array.data, Indices.astype('i'), Values)
    if m is not numpy.ma.nomask:
        if mv is not numpy.ma.nomask:
            genutil.array_indexing_emulate.set(m, Indices, mv)
    elif mv is not numpy.ma.nomask:
        m = numpy.zeros(mv.shape, mv.typcode())
        genutil.array_indexing_emulate.set(m, Indices, mv)
        if not numpy.ma.allequal(m, 0):
            Array = numpy.ma.masked_where(m, Array, copy=0)
    if ax is not None:
        Array = cdms2.createVariable(Array, axes=ax, id=id, copy=0)
        for at in list(xatt.keys()):
            setattr(Array, at, xatt[at])
    return Array

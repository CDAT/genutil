# Adapted for numpy/ma/cdms2 by convertcdms.py
# Adapted for numpy/ma/cdms2 by convertcdms.py
import numpy
import cdms2
import MV2
import genutil


def custom1D(x, filter, axis=0):
    """
    Apply a custom 1 dimensional filter to an array over a specified axis
    filter can be a list of numbers or a 1D array

    :param x: A CDMS TransientVariable
    :type x: cdms.tvariable.TransientVariable

    :param filter: numpy or MV array
    :type filter: array

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: str or int
    """
    isMV2 = cdms2.isVariable(x)
    if isMV2:
        xatt = x.attributes
    filter = MV2.array(filter)
    newx = MV2.array(x)
    initialorder = newx.getOrder(ids=1)
    n = len(filter)
    newx = newx(order=str(axis) + '...')
    sh = list(newx.shape)
    sh[0] = sh[0] - n + 1
    out = numpy.ma.zeros(sh, dtype=newx.dtype.char)
    ax = []
    bnds = []
    nax = newx.getAxis(0)
    for i in range(sh[0]):
        sub = newx[i:i + n]
        if i == 0:
            filter.setAxis(0, sub.getAxis(0))
            filter, sub = genutil.grower(filter, sub)
        out[i] = numpy.ma.average(sub, weights=filter, axis=0)
        if isMV2:
            a = nax.subAxis(i, i + n)
            try:
                b = a.getBounds()
                b1 = b[0][0]
                b2 = b[-1][1]
                ax.append((b1 + b2) / 2.)
                bnds.append([b1, b2])
            except BaseException:  # No bounds on this axis
                bnds = None
                ax.append(float(numpy.ma.average(a[:], axis=0)))
    out = MV2.array(out, id=newx.id)
    if isMV2:
        for k in list(xatt.keys()):
            setattr(out, k, xatt[k])
        for i in range(1, len(sh)):
            out.setAxis(i, newx.getAxis(i))
        if bnds is not None:
            bnds = numpy.ma.array(bnds)
        ax = cdms2.createAxis(ax, bounds=bnds)
        a = newx.getAxis(0)
        attr = a.attributes
        ax.id = a.id
        for k in list(attr.keys()):
            setattr(ax, k, attr[k])
        out.setAxis(0, ax)

    out = out(order=initialorder)
    if not isMV2:
        out = numpy.ma.array(out)
    return out


def smooth121(x, axis=0):
    """
    Apply a 121 filter to an array over a specified axis

    :Example:

        .. doctest:: genutil_filters_smooth121

            >>> import cdms2, vcs
            >>> vcs.download_sample_data()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> unfiltered=f('clt')
            >>> filtered = smooth121(unfiltered) # filter over axis at index 0

    :param x: A CDMS TransientVariable
    :type x: cdms.tvariable.TransientVariable

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: str or int
    """
    return custom1D(x, [1., 2., 1.], axis=axis)


def runningaverage(x, N, axis=0):
    """
    Apply a running average of length N to an array over a specified axis

    :Example:

        .. doctest:: genutil_filters_runningaverage

            >>> import cdms2, vcs
            >>> vcs.download_sample_data()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('clt')
            >>> smooth = runningaverage(x,12)

    :param x: A CDMS TransientVariable
    :type x: cdms.tvariable.TransientVariable

    :param N: length of the running average
    :type N: int

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: str or int
    """
    filter = numpy.ma.ones((N,), dtype='f')
    return custom1D(x, filter, axis=axis)

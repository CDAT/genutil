import numpy
# Adapted for numpy/ma/cdms2 by convertcdms.py


def minmax(*data):
    """
    Returns the minimum and maximum of a series of arrays/lists/tuples (or a combination of these)
    You can combine list/tuples/... pretty much any combination is allowed.

    :Example:

        .. doctest:: genutil_minmax

            >>> import genutil
            >>> s = range(7)
            >>> genutil.minmax(s)
            (0.0, 6.0)
            >>> genutil.minmax([s,s])
            (0.0, 6.0)
            >>> genutil.minmax([[s,s*2],4.,[6.,7.,s]],[5.,-7.,8,(6.,1.)])
            (-7.0, 8.0)
    """
    mx = numpy.finfo(numpy.float).min
    mn = numpy.finfo(numpy.float).max
    if len(data) == 1:
        data = data[0]
    global myfunction

    def myfunction(d, mx, mn):
        from numpy.ma import maximum, minimum, count
        if isinstance(d, (int, float)):
            return maximum(d, mx), minimum(d, mn)
        try:
            if count(d) == 0:
                return mx, mn
            mx = float(maximum(mx, maximum.reduce(d, axis=None)))
            mn = float(minimum(mn, minimum.reduce(d, axis=None)))
        except BaseException:
            for i in d:
                mx, mn = myfunction(i, mx, mn)
        return mx, mn
    mx, mn = myfunction(data, mx, mn)
    if mn == 1.E500 and mx == -1.E500:
        mn = mx = 1.E500
    return mn, mx

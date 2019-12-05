# Adapted for numpy/ma/cdms2 by convertcdms.py
import MV2
import numpy.ma
import cdms2
from .grower import grower
import numpy
import cdat_info
from . import arrayindexing
from . import array_indexing_emulate as array_indexing
from .stats_checker import __checker, StatisticsError


def __gammln1(x):
    cof = [76.18009172947146, -86.50532032941677,
           24.01409824083091, -1.231739572450155,
           .1208650973866179E-2, -.5395239384953E-5]
    stp = 2.5066282746310005
    y = x * 1.
    tmp = x + 5.5
    tmp = (x + 0.5) * numpy.ma.log(tmp) - tmp
    ser = 1.000000000190015
    for j in range(6):
        y = y + 1.
        ser = ser + cof[j] / y
    return tmp + numpy.ma.log(stp * ser / x)


def __betacf1(a, b, x):
    MAXIT = 100
    EPS = 3.E-7
    FPMIN = 1.E-30
    qab = a + b
    qap = a + 1.
    qam = a - 1.
    c = 1.
    d = 1. - qab * x / qap
    d = numpy.ma.where(numpy.ma.less(numpy.ma.absolute(d), FPMIN), FPMIN, d)
    d = 1. / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1. + aa * d
        d = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(d),
                FPMIN),
            FPMIN,
            d)
        c = 1. + aa / c
        c = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(c),
                FPMIN),
            FPMIN,
            c)
        d = 1. / d
        h = h * d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1. + aa * d
        d = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(d),
                FPMIN),
            FPMIN,
            d)
        c = 1. + aa / c
        c = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(c),
                FPMIN),
            FPMIN,
            c)
        d = 1. / d
        delet = d * c
        h = h * delet
        if numpy.ma.allclose(delet, numpy.ones(
                delet.shape), atol=EPS, rtol=0.):
            break
    h = numpy.ma.masked_where(
        numpy.ma.greater(
            numpy.ma.absolute(
                delet - 1.), EPS), h)
    return h


def __betai1(a, b, x):
    bt = numpy.ma.logical_or(numpy.ma.equal(x, 0.), numpy.ma.equal(x, 1.))
    bt = numpy.ma.where(bt, 0., numpy.ma.exp(
        __gammln1(a + b) - __gammln1(a) - __gammln1(b) +
        a * numpy.ma.log(x) + b * numpy.ma.log(1. - x)
    )
    )
    return numpy.ma.where(numpy.ma.less(x, (a + 1.) / (a + b + 2.)),
                          bt * __betacf1(a, b, x) / a,
                          1. - bt * __betacf1(b, a, 1. - x) / b)


def __probnd1(x):
    """
    FUNCTION PROBND1.

    Calculates the area under a normal curve (mean=0.0, variance=1.0)
    to the right of x. The accuracy is better than 7.5 * 10.**-8.

    REFERENCE:

    M. Abramowitz and I.A. Stegun.
    Handbook of Mathematical Functions.
    Dover, 1970, pages 931-932 (26.2.1 and 26.2.17).
"""
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    p = 0.2316419
    t = 1.0 / (1.0 + (p * x))
    term1 = ((((b1 * t) + (b2 * (t ** 2))) +
              (b3 * (t ** 3))) + (b4 * (t ** 4))) + \
        (b5 * (t ** 5))
    z = (1.0 / numpy.ma.sqrt(2.0 * numpy.pi)) * numpy.ma.exp(- ((x * x) / 2.0))
    return numpy.ma.where(numpy.ma.greater(x, 7.), 0., z * term1)


def __probf1(y, n1, n2, id):
    """
    FUNCTION PROBF1.

    The output is either the one- or two-tailed test area: i.e., the
    area under an F-curve (with N1 and N2 degrees of freedom) to the
    right of X if X exceeds 1.0 (one-tailed test) or twice this area
    (two-tailed test).

    Note: if X is less than 1.0, this function gives the area to the
    right of 1/X with reversed order for the degrees of freedom. This
    ensures the accuracy of the numerical algorithm.

    REFERENCE:

    M. Abramowitz and I.A. Stegun.
    Handbook of Mathematical Functions.
    Dover, 1970, page 947 (26.6.15).

    ** INPUT **
    real y            Calculated F-value
    real x            Inverse of Y if Y is less than 1.0
    integer n1, n2    Degrees of freedom
    integer id        Identifier for one- or two-tailed test

    ** OUTPUT **
    real probf1       Significance level (p-value) for F-value

    EXTERNALS:

    function PROBND1 - Calculates the area under a normal curve.
    """
    ly = numpy.ma.less(y, 1.)
    x = numpy.ma.where(ly, 1. / numpy.ma.array(y), y)
    n = numpy.ma.where(ly, n1, n2)
    n1 = numpy.ma.where(ly, n2, n1)
    n2 = numpy.ma.where(ly, n, n2)
    term1 = 2.0 / (9.0 * n1)
    term2 = 2.0 / (9.0 * n2)
    term3 = ((x ** (1.0 / 3.0)) * (1.0 - term2)) - (1.0 - term1)
    term4 = numpy.ma.sqrt(term1 + ((x ** (2.0 / 3.0)) * term2))
    term5 = term3 / term4
    probf1 = id * __probnd1(term5)

    #     The numerical algorithm can have problems when the F-value is
    #     close to 1.0 and the degrees of freedom are small. Therefore,
    #     insure that the probabilities returned cannot exceed 1.0.

    return numpy.ma.where(numpy.ma.greater(probf1, 1.), 1., probf1)


def __geometricmean(x):
    """
    Function: __geom
    Description:
       Returns the geometric mean (on first dimension)
    Usage:
    geo=__geometricmean(x)
    """
    g = numpy.ma.exp(
        numpy.ma.sum(
            numpy.ma.log(x) /
            numpy.ma.count(
                x,
                axis=0),
            0))
    # Now check for negative values
    t = numpy.ma.sum(numpy.ma.less(x, 0.), axis=0)
    g = numpy.ma.masked_where(numpy.ma.greater(t, 0.), g)
    t = numpy.ma.sum(numpy.ma.equal(x, 0), axis=0)
    g = numpy.ma.where(numpy.ma.greater(t, 0), 0., g)
    return g


def _treat_missing(out, x, max_pct_missing=100.):
    """
    # Max_pct_missing to missing data specified, calculate new mask
    # If fraction of data that is missing is more than the max_pct_missing allows mask point out
    """
    xmask = numpy.ma.less_equal((max_pct_missing / 100.), (numpy.ma.size(x, axis=0) -
                                                           numpy.ma.count(x, axis=0)) / float(numpy.ma.size(x, axis=0)))
    # numpy.ma.sum ignores missing_data, mask out results where amount of missing_data
    # is over the specified max_pct_missing
    return numpy.ma.masked_where(xmask, out)


def __covariance(x, y, weights=None, centered=1, biased=1):
    """
    Function: __covariance

    Description of function:
        Does the main computation for returning covariance. See documentation
        of covariance() for details.
    """
    if weights is not None and biased != 1:
        raise StatisticsError(
            'Error in covariance, you cannot have weights and unbiased together')

    if centered == 1:
        xmean = numpy.ma.average(x, weights=weights, axis=0)
        ymean = numpy.ma.average(y, weights=weights, axis=0)
        x = x - xmean
        y = y - ymean
        del(xmean)
        del(ymean)
    #
    if weights is None:
        weights = numpy.ma.ones(x.shape, dtype=x.dtype.char)
    if not ((x.mask is None) or (x.mask is MV2.nomask)):
        weights = numpy.ma.masked_where(x.mask, weights)
    if biased == 1:
        cov = numpy.ma.sum(x * y * weights, axis=0) / \
            numpy.ma.sum(weights, axis=0)
    else:
        cov = numpy.ma.sum(x * y, axis=0) / (numpy.ma.count(x * y, axis=0) - 1)

    return cov


def __variance(x, weights=None, centered=1, biased=1):
    """
    Function: __variance

    Description of function:
        Does the main computation for returning variance. See documentation
        of variance() for details.
    """
    return __covariance(x, x, weights=weights,
                        centered=centered, biased=biased)


def __std(x, weights=None, centered=1, biased=1):
    """
    Function: __std

    Description of function:
        Does the main computation for returning standard deviation. See
        documentation of std() for details.
    """
    return numpy.ma.sqrt(__variance(x, weights=weights,
                                    centered=centered, biased=biased))


def __correlation(x, y, weights=None, centered=1, biased=1):
    """
    Function: __correlation

    Description of function:
        Does the main computation for returning correlation. See documentation
        of correlation() for details.
    """
    cov = __covariance(x, y, weights=weights, centered=centered, biased=biased)
    sx = __std(x, weights=weights, centered=centered, biased=biased)
    sy = __std(y, weights=weights, centered=centered, biased=biased)
    return cov / (sx * sy)


def __rms(x, y, weights=None, centered=0, biased=1):
    """
    Function: __rms

    Description of function:
        Does the main computation for returning rms. See documentation
        of rms() for details.
    """

    return __std(x - y, centered=centered, biased=biased, weights=weights)


def __laggedcovariance(x, y, lag=1, centered=1, partial=1):
    """
    Function: __laggedcovariance

    Description of function:
        Does the main computation for returning lagged covariance. See
        documentation of laggedcovariance() for details.
    """
    if lag == 0:
        return __covariance(x, y, centered=centered)

    if partial == 1:
        if lag > 0:
            x = x[lag:]
            y = y[:-lag]
        else:
            x = x[:lag]
            y = y[-lag:]

    if centered == 1:
        xmean = numpy.ma.average(x, axis=0)
        ymean = numpy.ma.average(y, axis=0)
    else:
        xmean = 0.
        ymean = 0.
    x = x - xmean
    y = y - ymean
    del(xmean)
    del(ymean)

    if partial == 1:
        tmp = x * y
    else:
        if lag > 0:
            tmp = x[lag:] * y[:-lag]
        else:
            tmp = x[:-lag] * y[lag:]
    return numpy.ma.sum(tmp, axis=0) / numpy.ma.count(x * y, axis=0)


def __laggedcorrelation(x, y, lag, centered=1, partial=1, biased=1):
    """
    Function: __laggedcorrelation

    Description of function:
        Does the main computation for returning lagged correlation. See
        documentation of laggedcorrelation() for details.
    """
    cov = __laggedcovariance(x, y, lag, centered=centered, partial=partial)
    if partial == 1 and lag != 0:
        if lag > 0:
            sx = __std(x[lag:], centered=centered, biased=biased)
            sy = __std(y[:-lag], centered=centered, biased=biased)
        else:
            sx = __std(x[:lag], centered=centered, biased=biased)
            sy = __std(y[-lag:], centered=centered, biased=biased)

    else:
        sx = __std(x, centered=centered, biased=biased)
        sy = __std(y, centered=centered, biased=biased)

    return cov / (sx * sy)


def __autocovariance(x, lag, centered=1, partial=1):
    """
    Function: __autocovariance

    Description of function:
        Does the main computation for returning autocovariance. See
        documentation of autocovariance() for details.
    """
    return __laggedcovariance(x, x, lag, centered=centered, partial=partial)


def __autocorrelation(x, lag, centered=1, partial=1):
    """
    Function: __autocorrelation

    Description of function:
        Does the main computation for returning autocorrelation. See
        documentation of autocorrelation() for details.
    """
    if partial == 1 and centered == 1 and lag != 0:
        mean1 = numpy.ma.average(x[:-lag], axis=0)
        mean2 = numpy.ma.average(x[lag:], axis=0)
        x1 = x[:-lag] - mean1
        x2 = x[lag:] - mean2
        num = numpy.ma.sum(x1 * x2, axis=0)
        den = numpy.ma.sum(numpy.ma.power(x1, 2), axis=0) * \
            numpy.ma.sum(numpy.ma.power(x2, 2), axis=0)
        return num / numpy.ma.sqrt(den)
    else:
        return __autocovariance(x, lag, centered=centered, partial=partial) / \
            __autocovariance(x, 0, centered=centered, partial=partial)


def __meanabsdiff(x, y, weights=None, centered=1):
    """
    Function: __meanabsdiff

    Description of function:
        Does the main computation for returning mean absolute difference. See
        documentation of meanabsdiff() for details.
    """
    if centered == 1:
        xmean = numpy.ma.average(x, weights=weights, axis=0)
        ymean = numpy.ma.average(y, weights=weights, axis=0)
    else:
        xmean = 0.
        ymean = 0.
    x = x - xmean
    y = y - ymean
    del(xmean)
    del(ymean)
    if weights is None:
        weights = numpy.ma.ones(x.shape, dtype=x.dtype.char)
    if not ((x.mask is None) or (x.mask is MV2.nomask)):
        weights = numpy.ma.masked_where(x.mask, weights)
    if not ((y.mask is None) or (y.mask is MV2.nomask)):
        weights = numpy.ma.masked_where(y.mask, weights)

    return numpy.ma.sum(numpy.ma.absolute(x - y) * weights,
                        axis=0) / numpy.ma.sum(weights, axis=0)


def __linearregression(y, x, error=None, probability=None,
                       noslope=None, nointercept=None):
    """
    returns slope/intercept for linear regression of dependant var y and indep var x
    also possibly returns error and P values.
    """
    if (not (noslope is None or noslope == 0)) and (
            not (nointercept is None or nointercept == 0)):
        raise StatisticsError(
            'Error in __linearregression, at least one of the following argument as to be None:' +
            'noslope or nointercept, you are requesting nothing back !')
    if (probability is not None) and (error is None):
        raise StatisticsError(
            'Error in __linearregression, error must not be None if probability is defined, probability is:' +
            str(probability))
    if error is not None:
        if error > 3:
            raise StatisticsError(
                "Error in __linearregression, error must be None (0), 1, ,2 or 3")
    xmean = numpy.ma.average(x, axis=0)
    ymean = numpy.ma.average(y, axis=0)
    x = x - xmean
    y = y - ymean
    xy = numpy.ma.sum(y * x, axis=0)
    xx = numpy.ma.sum(x * x, axis=0)
    slope = xy / xx
    intercept = ymean - slope * xmean
    V = []
    if noslope is None or noslope == 0:
        V.append(slope)
    if nointercept is None or nointercept == 0:
        V.append(intercept)
    if error is None or error == 0:
        return V
    elif error == 1:
        E = []
        n1 = numpy.ma.count(y, axis=0)
        # Unadjusted errors
        res = (y + ymean) - (intercept + (x + xmean) *
                             numpy.ma.resize(slope, numpy.ma.shape(y)))  # x2
        ssd = numpy.ma.sum(res * res, axis=0)
        amsd1 = ssd / (n1 - 2.)  # amsd1=ssd/idfd1
        if noslope is None or noslope == 0:
            E.append(numpy.ma.sqrt(amsd1 / xx))
        if nointercept is None or nointercept == 0:
            s1 = xmean * xmean / xx + 1. / n1
            E.append(numpy.ma.sqrt(amsd1 * s1))
        if probability is None or probability == 0:
            return V, E
        else:
            Pt1 = []
            Pt2 = []
            Pf1 = []
            Pf2 = []
            f = numpy.ma.sum(y * y, axis=0) - ssd  # ssr
            f = f / amsd1
            aa1 = n1 / 2.0
            if noslope is None or noslope == 0:
                tb1 = slope / E[0]
                xx3 = n1 / (n1 + (tb1 * tb1))
                Pt1.append(__betai1(aa1, .5, xx3))
                Pt2.append(None)
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            if nointercept is None or nointercept == 0:
                tb1 = V[-1] / E[-1]
                xx3 = n1 / (n1 + (tb1 * tb1))
                Pt1.append(__betai1(aa1, .5, xx3))
                Pt2.append(None)
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            return V, E, Pt1, Pt2, Pf1, Pf2
    else:
        E = []
        # Adjusted error from residual
        n1 = numpy.ma.count(y, axis=0)
        res = (y + ymean) - (intercept + numpy.ma.resize(slope,
                                                         numpy.ma.shape(y)) * (x + xmean))  # x2
        ssd = numpy.ma.sum(res * res, axis=0)
        if error == 2:
            ac = __autocorrelation(res, 1, centered=1, partial=0)
            rdfd2 = 1.0 + ac
            rdfd2 = (1.0 - ac) / rdfd2
        elif error == 3:
            ac = __autocorrelation(y + ymean, 1, centered=1, partial=0)
            rdfd2 = 1.0 + ac
            rdfd2 = (1.0 - ac) / rdfd2
        rneff = n1 * rdfd2  # rneff
        amsd2 = ssd / (rneff - 2.)   # ssd/rdfd2
        if noslope is None or noslope == 0:
            E.append(numpy.ma.sqrt(amsd2 / xx))
        if nointercept is None or nointercept == 0:
            s1 = xmean * xmean / xx + 1. / n1
            E.append(numpy.ma.sqrt(amsd2 * s1))
        if probability is None or probability == 0:
            return V, E
        else:
            Pt1 = []
            Pt2 = []
            Pf1 = []
            Pf2 = []
            f = numpy.ma.sum(y * y, axis=0) - ssd  # ssr = amsr
            amsd1 = ssd / (n1 - 2.)  # amsd1=ssd/idfd1
            f = f / amsd1  # amsr/ssd
            aa1 = n1 / 2.0
            aa2 = rneff / 2.0
            if noslope is None or noslope == 0:
                tb2 = slope / E[0]
                xx1 = n1 / (n1 + (tb2 * tb2))
                xx2 = rneff / (rneff + (tb2 * tb2))
                Pt1.append(__betai1(aa1, .5, xx1))
                Pt2.append(__betai1(aa2, .5, xx2))
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            if nointercept is None or nointercept == 0:
                tb2 = V[-1] / E[-1]
                xx1 = n1 / (n1 + (tb2 * tb2))
                xx2 = rneff / (rneff + (tb2 * tb2))
                Pt1.append(__betai1(aa1, .5, xx1))
                Pt2.append(__betai1(aa2, .5, xx2))
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            return V, E, Pt1, Pt2, Pf1, Pf2


def covariance(x, y, weights=None, axis=0, centered=1,
               biased=1, max_pct_missing=100.):
    """
    Returns the covariance between 2 slabs. By default on the first dimension,
    centered and biased by default.

    :Example:

        .. doctest:: statistics_covariance

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> y=f('v')
            >>> cov = covariance(x, y) # default covariance of x and y

    :param x: The first slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param y: The second slab
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param weights: list of weights for assessing the weighted covariance

        .. note::

            Weighted covariance is inherently biased. Passing a value for weights but
            specifying an unbiased variance will cause an error

    :type weights: list
    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str
    :param centered: Integer flag for whether the covariance should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int
    :param biased: Flag indicating whether covariance should be calculated with bias. A value of 1 indicates covariance
        should be biased, while any other value indicates that it should not.
    :type biased: int
    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.covariance")
    if cdms2.isVariable(x):
        xatt = x.attributes
    if cdms2.isVariable(y):
        yatt = y.attributes
    x, y, weights, axis, ax = __checker(x, y, weights, axis)

    cov = __covariance(x, y, weights=weights, centered=centered, biased=biased)
    cov = _treat_missing(cov, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        cov = cdms2.createVariable(cov, axes=ax, id='covariance', copy=0)
        if 'units' in list(xatt.keys()) and 'units' in list(yatt.keys()):
            cov.units = xatt['units'] + '*' + yatt['units']
    return cov


def variance(x, weights=None, axis=0, centered=1,
             biased=1, max_pct_missing=100.):
    """
    Returns the variance from a slab. By default  on first dimension,
    centered, and biased.

    :Example:

        .. doctest:: statistics_variance

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> result = variance(x) # default variance of x

    :param x: Slab to compute variance of.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param weights: list of weights(floats) for assessing the weighted variance

        .. note::

            Weighted covariance is inherently biased. Passing a value for weights but
            specifying an unbiased variance will cause an error

    :type weights: list

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the variance should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param biased: Flag indicating whether variance should be calculated with bias. A value of 1 indicates variance
        should be biased, while any other value indicates that it should not.
    :type biased: int

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.variance")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, dum, weights, axis, ax = __checker(x, None, weights, axis)

    var = __variance(x, weights=weights, centered=centered, biased=biased)
    var = _treat_missing(var, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        var = cdms2.createVariable(var, axes=ax, id='variance', copy=0)
        if 'units' in list(xatt.keys()):
            var.units = xatt['units'] + '*' + xatt['units']
    return var


def checker(x, weights=None, axis=0, centered=1):
    x, dum, weights, axis, ax = __checker(x, None, weights, axis)
    return x, weights, axis, ax


def std(x, weights=None, axis=0, centered=1, biased=1, max_pct_missing=100.):
    """
    Returns the standard deviation from a slab. By default  on first
    dimension, centered, and biased.

    :Example:

        .. doctest:: statistics_std

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> result = std(x) # default standard deviation of x

    :param x: Slab to compute std of.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param weights: list of weights(floats) for assessing the weighted calculation

        .. note::

            Use of weights is inherently biased. Passing a value for weights but
            specifying a non-biased calculation will cause an error.

    :type weights: list

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param biased: Flag indicating whether std should be calculated with bias. A value of 1 indicates calculation
        should be biased, while any other value indicates that it should not.
    :type biased: int

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.std")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, dum, weights, axis, ax = __checker(x, None, weights, axis)
    std = __std(x, weights=weights, centered=centered, biased=biased)
    std = _treat_missing(std, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        std = cdms2.createVariable(
            std, axes=ax, id='standard_deviation', copy=0)
        if 'units' in list(xatt.keys()):
            std.units = xatt['units']
    return std


def correlation(x, y, weights=None, axis=0, centered=1,
                biased=1, max_pct_missing=100.):
    """
    Returns the correlation between 2 slabs. By default on the first
    dimension, centered and biased by default.

    :Example:

        .. doctest:: statistics_correlation

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> y=f('v')
            >>> result = correlation(x, y) # calculate the default correlation
    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param y: Second slab.
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param weights: list of weights(floats) for assessing the weighted calculation

        .. note::

            Use of weights is inherently biased. Passing a value for weights but
            specifying a non-biased calculation will cause an error.

    :type weights: list

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param biased: Flag indicating whether calculation should be biased. A value of 1 indicates calculation
        should be biased, while any other value indicates that it should not.
    :type biased: int

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.corelation")
    x, y, weights, axis, ax = __checker(x, y, weights, axis)

    cor = __correlation(
        x,
        y,
        weights=weights,
        centered=centered,
        biased=biased)
    cor = _treat_missing(cor, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        cor = cdms2.createVariable(cor, axes=ax, id='correlation', copy=0)
        cor.units = '-'
    return cor


def rms(x, y, weights=None, axis=0, centered=0,
        biased=1, max_pct_missing=100.):
    """
    Returns the root mean square difference between 2 slabs. By default from
    a slab (on first dimension) "uncentered" and "biased" by default.

    :Example:

        .. doctest:: statistics_rms

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> y=f('v')
            >>> result = rms(x, y) # default rms calculation
    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param y: Second slab.
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param weights: list of weights(floats) for assessing the weighted calculation

        .. note::

            Use of weights is inherently biased. Passing a value for weights but
            specifying a non-biased calculation will cause an error.

    :type weights: list

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param biased: Flag indicating whether calculation should be biased. A value of 1 indicates calculation
        should be biased, while any other value indicates that it should not.
    :type biased: int

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.rms")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, y, weights, axis, ax = __checker(x, y, weights, axis)
    rmsans = __rms(x, y, weights=weights, centered=centered, biased=biased)
    rmsans = _treat_missing(rmsans, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        rmsans = cdms2.createVariable(
            rmsans, axes=ax, id='RMS_difference', copy=0)
        if 'units' in list(xatt.keys()):
            rms.units = xatt['units']

    return rmsans


def laggedcovariance(x, y, lag=None, axis=0, centered=1,
                     partial=1, noloop=0, max_pct_missing=100.):
    """
    Returns the covariance between 2 slabs at lag k centered and partial by
    default.

    :Example:

        .. doctest:: statistics_laggedcovariance

            >>> import vcs, cdms2, os, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + '/clt.nc')
            >>> x=f('u')
            >>> y=f('v')
            >>> result = laggedcovariance(x, y) # default laggedcovariance calculation, with max. number of lags

    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param y: Second slab.
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param lag: Integer, list of integers, tuple of integers, or None. If None, maximum possible lags for specified axis
            are used.
    :type lag: int or list or tuple or None

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int
    :param partial: Integer flag. If 1, uses only common time for means.
    :type partial: int or None
    :param noloop: None | 0 | 1
        default value = 0 computes statistic at all lags upto 'lag'. If you
        set noloop=1 statistic is computed at lag only (not up to lag).
    :type noloop: int or None

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.laggedcovariance")
    if cdms2.isVariable(x):
        xatt = x.attributes
    if cdms2.isVariable(y):
        yatt = y.attributes
    x, y, w, axis, ax = __checker(x, y, None, axis)
    if lag is None:
        lags = list(range(x.shape[0]))
    elif isinstance(lag, int):
        if not noloop:
            lags = list(range(lag + 1))
        else:
            lags = [lag]
    elif not isinstance(lag, (list, tuple)):
        raise StatisticsError('lags type must be None, integer, list or tuple')
    else:
        lags = lag

    for k in lags:
        lcov = numpy.ma.array(
            __laggedcovariance(
                x,
                y,
                lag=k,
                centered=centered,
                partial=partial))
        lcov = _treat_missing(lcov, x, max_pct_missing=max_pct_missing)
        sh = list(lcov.shape)
        sh.insert(0, 1)
        lcov = numpy.ma.resize(lcov, sh)
        if k == lags[0]:
            lcovs = lcov
        else:
            lcovs = numpy.ma.concatenate((lcovs, lcov), 0)

    # lcovs=_treat_missing(lcovs,x,max_pct_missing=max_pct_missing)
    if ax is not None:
        newax = cdms2.createAxis(lags)
        newax.id = 'lag'
        ax.insert(0, newax)
        lcovs = cdms2.createVariable(
            lcovs,
            axes=ax,
            id='lagged_covariance' +
            str(lag),
            copy=0)
        if 'units' in list(xatt.keys()) and 'units' in list(yatt.keys()):
            lcovs.units = xatt['units'] + '*' + yatt['units']
    return lcovs


def laggedcorrelation(x, y, lag=None, axis=0, centered=1,
                      partial=1, biased=1, noloop=0, max_pct_missing=100.):
    """

    Returns the correlation between 2 slabs at lag k centered, partial and
    "biased" by default.

    :Example:

        .. doctest:: statistics_laggedcorrelation

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> y=f('v')
            >>> result = laggedcorrelation(x, y) # default laggedcorrelation, with max. number of lags


    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param y: Second slab.
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param lag: Integer, list of integers, tuple of integers, or None.
        If None, maximum possible lags for specified axis are used.
    :type lag: int or list or tuple or None

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int
    :param partial: Integer flag. If 1, uses only common time for means.
    :type partial: int or None
    :param noloop: None | 0 | 1
        default value = 0 computes statistic at all lags upto 'lag'. If you
        set noloop=1 statistic is computed at lag only (not up to lag).
    :type noloop: int or None

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.laggedcorrelation")
    x, y, w, axis, ax = __checker(x, y, None, axis)
    if lag is None:
        lags = list(range(x.shape[0]))
    elif isinstance(lag, int):
        if not noloop:
            lags = list(range(lag + 1))
        else:
            lags = [lag]
    elif not isinstance(lag, (list, tuple)):
        raise StatisticsError('lags type must be None, integer, list or tuple')
    else:
        lags = lag

    for k in lags:
        lcor = numpy.ma.array(
            __laggedcorrelation(
                x,
                y,
                lag=k,
                centered=centered,
                partial=partial,
                biased=biased))
        lcor = _treat_missing(lcor, x, max_pct_missing=max_pct_missing)
        sh = list(lcor.shape)
        sh.insert(0, 1)
        lcor = numpy.ma.resize(lcor, sh)
        if k == lags[0]:
            lcors = lcor
        else:
            lcors = numpy.ma.concatenate((lcors, lcor), 0)
    if ax is not None:
        newax = cdms2.createAxis(lags)
        newax.id = 'lag'
        ax.insert(0, newax)
        lcors = cdms2.createVariable(
            lcors,
            axes=ax,
            id='lagged_correlation' +
            str(lag),
            copy=0)
        lcors.units = '-'

    return lcors


def autocovariance(x, lag=None, axis=0, centered=1,
                   partial=1, noloop=0, max_pct_missing=100.):
    """
    Returns the autocovariance of a slab. By default over the first dimension,  centered, and partial.

    :Example:

        .. doctest:: statistics_autocovariance

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> result = autocovariance(x) # default autocovariance

    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param lag: Integer, list of integers, tuple of integers, or None.
        If None, maximum possible lags for specified axis are used.
    :type lag: int or list or tuple or None

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param partial: Integer flag. If 1, uses only common time for means.
    :type partial: int or None

    :param noloop: None | 0 | 1
        default value = 0 computes statistic at all lags upto 'lag'. If you
        set noloop=1 statistic is computed at lag only (not up to lag).
    :type noloop: int or None

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.autocovariance")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, dum, dum, axis, ax = __checker(x, None, None, axis)
    if lag is None:
        lags = list(range(x.shape[0]))
    elif isinstance(lag, int):
        if not noloop:
            lags = list(range(lag + 1))
        else:
            lags = [lag]
    elif not isinstance(lag, (list, tuple)):
        raise StatisticsError('lags type must be None, integer, list or tuple')
    else:
        lags = lag

    for k in lags:
        acov = numpy.ma.array(
            __autocovariance(
                x,
                lag=k,
                centered=centered,
                partial=partial))
        acov = _treat_missing(acov, x, max_pct_missing=max_pct_missing)
        sh = list(acov.shape)
        sh.insert(0, 1)
        acov = numpy.ma.resize(acov, sh)
        if k == lags[0]:
            acovs = acov
        else:
            acovs = numpy.ma.concatenate((acovs, acov), 0)
    if ax is not None:
        newax = cdms2.createAxis(lags)
        newax.id = 'lag'
        ax.insert(0, newax)
        acovs = cdms2.createVariable(
            acovs, axes=ax, id='autocovariance' + str(lag), copy=0)
        if 'units' in list(xatt.keys()):
            acovs.units = xatt['units'] + '*' + xatt['units']
    return acovs


def autocorrelation(x, lag=None, axis=0, centered=1, partial=1,
                    biased=1, noloop=0, max_pct_missing=100.):
    """
    Returns the autocorrelation of a slab at lag k centered,partial and
    "biased" by default

    :Example:

        .. doctest:: statistics_autocorrelation

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> result = autocorrelation(x)
    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param lag: Integer, list of integers, tuple of integers, or None.
        If None, maximum possible lags for specified axis are used.
    :type lag: int or list or tuple or None

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param partial: Integer flag. If 1, uses only common time for means.
    :type partial: int or None

    :param noloop: None | 0 | 1
        default value = 0 computes statistic at all lags upto 'lag'. If you
        set noloop=1 statistic is computed at lag only (not up to lag).
    :type noloop: int or None

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.autocorrelation")
    x, dum, dum, axis, ax = __checker(x, None, None, axis)
    if lag is None:
        lags = list(range(x.shape[0]))
    elif isinstance(lag, int):
        if not noloop:
            lags = list(range(lag + 1))
        else:
            lags = [lag]
    elif not isinstance(lag, (list, tuple)):
        raise StatisticsError('lags type must be None, integer, list or tuple')
    else:
        lags = lag

    for k in lags:
        acov = numpy.ma.array(
            __autocorrelation(
                x,
                lag=k,
                centered=centered,
                partial=partial))
        acov = _treat_missing(acov, x, max_pct_missing=max_pct_missing)
        sh = list(acov.shape)
        sh.insert(0, 1)
        acov = numpy.ma.resize(acov, sh)
        if k == lags[0]:
            acovs = acov
        else:
            acovs = numpy.ma.concatenate((acovs, acov), 0)
    if ax is not None:
        newax = cdms2.createAxis(lags)
        newax.id = 'lag'
        ax.insert(0, newax)
        acovs = cdms2.createVariable(
            acovs, axes=ax, id='autocorrelation' + str(lag), copy=0)
        acovs.units = '-'

    return acovs


def meanabsdiff(x, y, weights=None, axis=0, centered=1, max_pct_missing=100.):
    """
    Returns the mean absolute difference between 2 slabs x and y. By default
    on the first dimension and centered

    :Example:

        .. doctest:: statistics_meanabsdiff

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('u')
            >>> y=f('v')
            >>> result = meanabsdiff(x, y) # default mean absolute difference calculation

    :param x: First slab.
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param y: Second slab.
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param weights: list of weights for assessing the weighted calculation
    :type weights: list

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param centered: Integer flag for whether the calculation should be centered.
        0 or None means don't center, 1 means center.
    :type centered: int

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.meanabsdiff")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, y, weights, axis, ax = __checker(x, y, weights, axis)
    mad = __meanabsdiff(x, y, weights=weights, centered=centered)
    mad = _treat_missing(mad, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        mad = cdms2.createVariable(
            mad, axes=ax, id='mean_absolute_difference', copy=0)
        if 'units' in list(xatt.keys()):
            mad.units = xatt['units']
    return mad


def linearregression(y, axis=None, x=None, error=None,
                     probability=None, nointercept=None, noslope=None):
    """
    Computes the linear regression of y over x or an axis. This function returns
    Values of the slope and intercept, and optionally, Error estimates and
    associated probability distributions for T-value (T-Test) and F-value (for
    analysis of variance f) can be returned. You can choose to return all these
    for either slope or intercept  or both (default behaviour). For theoretical
    details, refer to 'Statistical Methods in Atmospheric Sciences' by
    Daniel S. Wilks, Academic Press, 1995.

    :Example:

        .. doctest:: statistics_linearregression

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> y=f('v')
            >>> result = linearregression(y, axis='x') # default linearregression over x axis

    :param y: Slab to compute linear regression of
    :type y: cdms2.tvariable.TransientVariable or numpy.array

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param x: Slab over which linear regression of y will be computed
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param nointercept: Setting to 0 or None means intercept calculations
        are returned. To turn OFF the intercept computations set
        nointercept = 1.
    :type nointercept: int or None

    :param noslope: Setting to None or 0 means slope calculations are
        returned. To turn OFF the slope computations set noslope to 1.
    :type noslope: int or None

    :param error: If set to 0 or None, no associated errors are
        returned.
        If set to 1, the unadjusted standard error is returned.
        If set to 2, standard error returned. This standard error is
            adjusted using the centered autocorrelation of the residual.
        If set to 3, standard error returned. The standard error here is
            adjusted using the centered autocorrelation of the raw data (y).
    :type error: int or None

    :param probability: If set to 0 or None, no associated probabilities are
            returned. Set this to 1 to compute probabilities.

        .. note::

            Probabilities are returned only if erroroptions are set to one
            of 1, 2, or 3. If it is set to None or 0, then setting
            probabilityoptions has no meaning.
    :type probability: int or None

    :Return Values:

        The returned values depend on the combination of options you select. If
        both slope and intercept are required, a tuple is returned for both Value
        and optionally Error (or optionally associated Probabilities), but single
        values (not tuples) are returned if only one set (slope OR intercept) is
        required. See examples below for more details.

        * When erroroption = 1 (from description above for erroroptions you know
            that means unadjusted standard error) and probabilityoption=1, then
            the following are returned:
            pt1 : The p-value for regression coefficient t-value. (With no
                  adjustment for standard error or critical t-value.)
            None: There is only one p-value to be returned (pf1) but None is
                  returned to keep the length of the returned values consistent.
            pf1 : The p-value for regression coefficient F-value (one-tailed).
            pf2 : The p-value for regression coefficient F-value (two-tailed).

        * When erroroption = 2 or 3 (implying error adjustment using the residual
            or the raw data and probabilityoption = 1, then the following are
            returned:
            pt1 : The p-value for regression coefficient t-value.(With effective
                  sample size adjustment for standard error of slope.
            pt2 : The p-value for regression coefficient t-value.(With effective
                  sample size adjustment for standard error of slope and critical
                  t-value.)
            pf1 : The p-value for regression coefficient F-value (one-tailed).
            pf2 : The p-value for regression coefficient F-value (two-tailed).

        * The values pt1 and pt2 are used to test the null hypothesis that b = 0
            (i.e., y is independent of x).
        * The values pf1 and pf2 are used to test the null hypothesis that the
            regression is linear (goodness of linear fit). For non-replicated
            values of y, the degrees of freedom are 1 and n-2.


        :Examples:

            .. code-block:: python

                # Let us first examine the default behaviour of the linearregression
                # function.
                >>> Values = statistics.linearregression(y)

                # The returned "Values" is actually a tuple consisting of the slope and
                # intercept. They can also be accessed as follows:
                >>> slope, intercept = statistics.linearregression(y)

                # If error estimates are also required, then:
                >>> Values, Errors = linearregression(y, error=1)

                # where "Values" and "Errors" are tuples containing answer for
                # slope AND intercept. You can break them as follows.
                # slope, intercept = Value and slope_error, intercept_error = Errors. i.e.
                >>> (slope, intercept), (slo_error, int_error) = \
                                    linearregression(y, error=1)

                # WARNING: The following will not work.
                >>> slope, intercept, slo_error, int_error = linearregression(y, error=1)

                # To get the standard error non adjusted result for slope only
                >>> slope, slope_error = linearregression(y, error=1, nointercept=1)

                # In the line below all the returned values are tuples.
                >>> Values,Errors,Pt1,Pt2,Pf1,Pf2 = \
                             linearregression(y, error=2,probability=1)
                # That means in the above statement is returning tuples ordered so:
                # (slope, intercept), (slo_error, int_error), (pt1_slo, pt1_int),\
                             (pt2_slo, pt2_int), (pf1_slo, pf1_int), (pf2_slo, pf2_int)

                # If we want results returned for the intercept only.
                >>> intercept,intercept_error,pt1,pt2,pf1,pf2 = \
                          linearregression(y, error=2,probability=1,noslope=1)

    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.linearregression")
    yisV = cdms2.isVariable(y)
    if yisV:
        yatt = y.attributes
    if axis is not None and x is not None:
        raise StatisticsError(
            'Error you cannot pass an indepedent variable and an axis')
    if x is None and axis is None:
        axis = 0
    if axis is not None:
        if not isinstance(axis, type([])):
            ax = cdms2.orderparse(str(axis))
        if len(ax) > 1:
            raise StatisticsError('Error only one dim allowed')
        if yisV:
            ax = y.getAxisList()
            yid = y.id
        y, x, w, axis, axs = __checker(y, x, None, axis)
        if yisV:
            ax2 = ax.pop(axis[0])
            ax.insert(0, ax2)
            y = cdms2.createVariable(y, axes=ax, id=yid, copy=0)
            x = cdms2.MV2.array(ax2[:])
            x.setAxis(0, ax2)
            y, x, w, axis, axs = __checker(y, x, None, 0)
        else:
            x = numpy.ma.arange(y.shape[0])
            x = numpy.ma.resize(x, y.shape)
    else:
        y, x, w, axis, axs = __checker(y, x, None, 0)
        ax2 = cdms2.asVariable(x)
        try:
            ax2.units = x.units
        except BaseException:
            ax2.units = ''

    if error is None or error == 0:
        val = __linearregression(
            y, x, nointercept=nointercept, noslope=noslope)
    elif probability is None or probability == 0:
        val, err = __linearregression(
            y, x, error=error, nointercept=nointercept, noslope=noslope)
    else:
        val, err, pt1, pt2, pf1, pf2 = __linearregression(
            y, x, error=error, probability=probability, nointercept=nointercept, noslope=noslope)
    if yisV:
        if noslope is None or noslope == 0:
            val[0] = cdms2.createVariable(val[0], axes=axs, id='slope', copy=0)
        if nointercept is None or nointercept == 0:
            val[-1] = cdms2.createVariable(val[-1],
                                           axes=axs, id='intercept', copy=0)
        if 'units' in list(yatt.keys()):
            for v in val:
                v.units = yatt['units'] + ' per ' + ax2.units
    if error is None or error == 0:
        if len(val) > 1:
            return val
        else:
            return val[0]
    elif probability is None or probability == 0:
        if yisV:
            if noslope is None or noslope == 0:
                err[0] = cdms2.createVariable(
                    err[0], axes=axs, id='standard_error', copy=0)
                if error == 1:
                    setattr(
                        err[0],
                        'long_name',
                        'standard error for regression coefficient')
                elif error == 2:
                    setattr(
                        err[0],
                        'long_name',
                        'standard error for regression coefficient ' +
                        'adjusted with residual (using centered autocorrelation)')
                elif error == 3:
                    setattr(
                        err[0],
                        'long_name',
                        'standard error for regression coefficient adjusted with y (using centered autocorrelation)')

            if nointercept is None or nointercept == 0:
                err[-1] = cdms2.createVariable(err[-1],
                                               axes=axs, id='standard_error', copy=0)
                if error == 1:
                    setattr(err[-1], 'long_name',
                            'standard error for regression constant')
                elif error == 2:
                    setattr(
                        err[-1], 'long_name', 'standard error for regression constant adjusted with residual ' +
                        '(using centered autocorrelation)')
                elif error == 3:
                    setattr(err[-1],
                            'long_name',
                            'standard error for regression constant adjusted with y (using centered autocorrelation)')
            if 'units' in list(yatt.keys()):
                for e in err:
                    e.units = yatt['units'] + ' per ' + ax2.units
        if len(val) > 1:
            return val, err
        else:
            return val[0], err[0]
    else:
        if yisV:
            if noslope is None or noslope == 0:
                err[0] = cdms2.createVariable(
                    err[0], axes=axs, id='standard_error', copy=0)
                if error == 1:
                    setattr(
                        err[0],
                        'long_name',
                        'standard error for regression coefficient')
                elif error == 2:
                    setattr(
                        err[0],
                        'long_name',
                        'standard error for regression coefficient adjusted with residual' +
                        '(using centered autocorrelation)')
                elif error == 3:
                    setattr(
                        err[0],
                        'long_name',
                        'standard error for regression coefficient adjusted with y (using centered autocorrelation)')

            if nointercept is None or nointercept == 0:
                err[-1] = cdms2.createVariable(err[-1],
                                               axes=axs, id='standard_error', copy=0)
                if error == 1:
                    setattr(err[-1], 'long_name',
                            'standard error for regression constant')
                elif error == 2:
                    setattr(
                        err[-1], 'long_name', 'standard error for regression constant adjusted with residual ' +
                        '(using centered autocorrelation)')
                elif error == 3:
                    setattr(err[-1],
                            'long_name',
                            'standard error for regression constant adjusted with y (using centered autocorrelation)')
            if 'units' in list(yatt.keys()):
                for e in err:
                    e.units = yatt['units'] + ' per ' + ax2.units
            if noslope is None or noslope == 0:
                if error > 1:
                    pt1[0] = cdms2.createVariable(
                        pt1[0], axes=axs, id='p-value', copy=0)
                    pt1[0].units = '-'
                    pt1[0].long_name = 'p-value for regression coefficient t-value.' +\
                        'Effective sample size adjustment for standard error (seb).'
                    pt2[0] = cdms2.createVariable(
                        pt2[0], axes=axs, id='p-value', copy=0)
                    pt2[0].units = '-'
                    pt2[
                        0].long_name = 'p-value for regression coefficient t-value. ' +\
                        'Effective sample size adjustment for standard error (seb) and critical t-value.'
                else:
                    pt1[0] = cdms2.createVariable(
                        pt1[0], axes=axs, id='p-value', copy=0)
                    pt1[0].units = '-'
                    pt1[0].long_name = 'p-value for regression coefficient t-value. ' +\
                        'No adjustment for standard error or critical t-value.'
                pf1[0] = cdms2.createVariable(
                    pf1[0], axes=axs, id='p-value', copy=0)
                pf1[0].unit = '-'
                pf1[0].long_name = 'p-value for regression coefficient F-value (one-tailed)'
                pf2[0] = cdms2.createVariable(
                    pf2[0], axes=axs, id='p-value', copy=0)
                pf2[0].unit = '-'
                pf2[0].long_name = 'p-value for regression coefficient F-value (two-tailed)'
            if nointercept is None or nointercept == 0:
                if error > 1:
                    pt1[-1] = cdms2.createVariable(pt1[-1],
                                                   axes=axs, id='p-value', copy=0)
                    pt1[-1].units = '-'
                    pt1[-1].long_name = 'p-value for regression coefficient t-value. ' +\
                        'Effective sample size adjustment for standard error (seb).'
                    pt2[-1] = cdms2.createVariable(pt2[-1],
                                                   axes=axs, id='p-value', copy=0)
                    pt2[-1].units = '-'
                    pt2[-1].long_name = 'p-value for regression coefficient t-value. ' +\
                        'Effective sample size adjustment for standard error (seb) and critical t-value.'
                else:
                    pt1[-1] = cdms2.createVariable(pt1[-1],
                                                   axes=axs, id='p-value', copy=0)
                    pt1[-1].units = '-'
                    pt1[-1].long_name = 'p-value for regression coefficient t-value. ' +\
                        'No adjustment for standard error or critical t-value.'
                pf1[-1] = cdms2.createVariable(pf1[-1],
                                               axes=axs, id='p-value', copy=0)
                pf1[-1].unit = '-'
                pf1[-1].long_name = 'p-value for regression coefficient F-value (one-tailed)'
                pf2[-1] = cdms2.createVariable(pf2[-1],
                                               axes=axs, id='p-value', copy=0)
                pf2[-1].unit = '-'
                pf2[-1].long_name = 'p-value for regression coefficient F-value (two-tailed)'
        if len(val) > 1:
            return val, err, pt1, pt2, pf1, pf2
        else:
            return val[0], err[0], pt1[0], pt2[0], pf1[0], pf2[0]


def geometricmean(x, axis=0, max_pct_missing=100.):
    """
    Returns the geometric mean over a specified axis.

    :Example:

        .. doctest:: statistics_geometricmean

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('v')
            >>> result = geometricmean(x, axis='x')

    :param x: Slab to compute geometric mean for
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str

    :param max_pct_missing: Maximum percentage of cell which is allowed to be masked (missing).
    :type max_pct_missing: float
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.geometricmean")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, dum, weights, axis, ax = __checker(x, None, None, axis)
    gmean = __geometricmean(x)
    gmean = _treat_missing(gmean, x, max_pct_missing=max_pct_missing)
    if ax is not None:
        gmean = cdms2.createVariable(
            gmean, axes=ax, id='geometric_mean', copy=0)
        if 'units' in list(xatt.keys()):
            gmean.units = xatt['units']
    return gmean


def _percentiles(out, percent):
    if cdms2.isVariable(out):
        out = MV2.sort(out, axis=0).asma()
        ns = MV2.count(out, axis=0).asma()
    else:
        out = numpy.ma.sort(out, axis=0)
        ns = numpy.ma.count(out, axis=0)

    output = None
    for p in percent:
        i = numpy.floor((p / 100.) * (ns - 1))
        try:
            i = i.astype(numpy.int)
        except BaseException:
            i = int(i)
        ii = i + 1
        tmp = (100. * i) / (ns - 1)
        try:
            tmp = tmp.filled(1.E20)
        except BaseException:
            pass
        Ai = numpy.where(numpy.equal(ns, 1), 0., tmp)
        tmp = (100. * ii) / (ns - 1)
        try:
            tmp = tmp.filled(1.E20)
        except BaseException:
            pass
        Aii = numpy.where(numpy.equal(ns, 1), 100., tmp)
        ii = numpy.where(numpy.equal(ii, ns), ns - 1, ii)
        if numpy.rank(ii) > 0:
            ii = ii.astype(numpy.int)
# tmp = (p-Ai)/(Aii-Ai)*array_indexing.extract(out,ii) + \
# (Aii-p)/(Aii-Ai)*array_indexing.extract(out,i)

        tmp = (p - Ai) / (Aii - Ai) * arrayindexing.get(out, ii) + \
            (Aii - p) / (Aii - Ai) * arrayindexing.get(out, i)
        try:
            sh = list(tmp.shape)
            sh.insert(0, 1)
            tmp = numpy.ma.reshape(tmp, sh)
        except BaseException:
            pass
        if output is None:
            output = numpy.ma.array(tmp)
            output = output.astype(out.dtype.char)
        else:
            output = numpy.ma.concatenate(
                (output, tmp.astype(numpy.float32)), 0)
    return output


def percentiles(x, percentiles=[50.], axis=0):
    """
    Returns values at the defined percentiles for an array.

    :Example:

        .. doctest:: statistics_percentiles

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('v')
            >>> result = percentiles(x, percentiles=[25.,50.,75.,100.]) # percentiles across default axis

    :param x: Slab to compute statistics for
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param percentiles: A python list of values
            Default = [50.] (the 50th percentile i.e the median value)
    :type percentiles: list

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
    :type axis: int or str
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.percentiles")
    if cdms2.isVariable(x):
        xatt = x.attributes
    x, dum, weights, axis, ax = __checker(x, None, None, axis)
    p = _percentiles(x, percentiles)
    if ax is not None and ax != []:
        pax = cdms2.createAxis(percentiles)
        pax.id = 'percentiles'
        pax.units = '%'
        ax.insert(0, pax)
        p = MV2.array(p)
        p = cdms2.createVariable(p, axes=ax, id='percentiles', copy=0)
        if 'units' in list(xatt.keys()):
            p.units = xatt['units']
    return p


def median(x, axis=0):
    """
    Returns the median value of an array.

    :Example:

        .. doctest:: statistics_median

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('v')
            >>> result = median(x, axis='x')

    :param x: Slab to compute statistics for
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
            default value = 0. You can pass the name of the dimension or index
            (integer value 0...n) over which you want to compute the statistic.
    :type axis: str or int
    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.median")
    tmp = percentiles(x, percentiles=[50.], axis=axis)
    tmp.id = 'median'
    return tmp


def rank(x, axis=0):
    """
    Returns the rank of each element along the specified axis
    where 0 is lowest value, 100 is maximum value. Handles missing values correctly

    :Example:

        .. doctest:: statistics_rank

            >>> import vcs, cdms2, os
            >>> try:
            ...    os.listdir(vcs.sample_data)
            >>> except:
            ...    vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sample_data + 'clt.nc')
            >>> x=f('v')
            >>> result = median(x, axis='x')

    :param x: Slab to compute statistics for
    :type x: cdms2.tvariable.TransientVariable or numpy.array

    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
            default value = 0. You can pass the name of the dimension or index
            (integer value 0...n) over which you want to compute the statistic.
    :type axis: str or int
    """

    # preprocessing
    #cdat_info.pingPCMDIdb("cdat", "genutil.statistics.rank")
    if cdms2.isVariable(x):
        xatt = x.attributes
        axs = x.getAxisList()
        o = x.getOrder(ids=1)
    x, dum, weights, axis, ax = __checker(x, None, None, axis)

    # First figures out indices to sort
    a0 = numpy.ma.array(numpy.ma.argsort(x.filled(1.E20), axis=0), dtype='i')
    n = a0.shape[0]

    # initialize output array
    b = numpy.ma.zeros(a0.shape, dtype='f')

    # Get the indices
    # Make sure b and a0 are of the right type
    b = array_indexing.rank(b, a0)
    m = x.mask
    if m is not numpy.ma.nomask:
        b = MV2.masked_where(m, b)
    else:
        b = MV2.array(b)
    n = MV2.count(b, axis=0)
    n.setAxisList(b.getAxisList()[1:])
    b, n = grower(b, n)
    b = 100. * b / (n - 1)

    # Now reorders everything
    if ax is not None:
        # First set the unchanged axes
        sh = []
        for i in range(len(ax)):
            sh.append(len(ax[i]))
        # Now figures the other axes to add
        for i in range(len(axis)):
            sh.insert(i, len(axs[axis[i]]))
        b = MV2.reshape(b, sh)
        for i in range(len(ax)):
            b.setAxis(i + len(axis), ax[i])
        for i in range(len(axis)):
            b.setAxis(i, axs[axis[i]])
        b = b(order=o)
        for a in list(xatt.keys()):
            if a[0] != '_':
                setattr(b, a, xatt[a])
        b.units = '%'
    elif len(axis) == 1:
        sh = list(range(b.rank()))
        sh[0] = axis[0]
        sh[axis[0]] = 0
        b = numpy.ma.transpose(b, sh)
    return b

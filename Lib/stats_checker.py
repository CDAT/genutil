import cdms2
import numpy
from .grower import grower
from .averager import __check_weightoptions


class StatisticsError(Exception):
    pass


def __makeweights(x, w, axes):
    """
    This function takes an array and weights options from Krishna\'s averager
    and return an numpy.ma of the coresponding weights
    """

    # Now if weights is a list, uses Krishna's stuff to get the weights....
    tmpaxes = axes
    if isinstance(axes, type(1)):
        tmpaxes = str(axes)
    # First make sure x and w have same dims if w is MV2
    if cdms2.isVariable(w) and cdms2.isVariable(x) and x.shape != w.shape:
        x, w = grower(x, w)
    w = __check_weightoptions(x, tmpaxes, w)
    if not numpy.ma.isarray(w):
        # Ok Krishna returned a list of 1D arrays.... Let's put it together
        axs = x.getAxisList()
        axes = cdms2.order2index(axs, axes)[: len(cdms2.orderparse(axes))]
        endax = []
        for i in range(len(axes)):
            if w[i] == "unweighted":
                w[i] = numpy.ma.ones(len(axs[axes[i]]), dtype=x.dtype.char)
            if i == 0:
                wo = w[i]
                endax.append(axs[axes[i]])
            else:
                wo = wo[..., None] * w[i]
                endax.append(axs[axes[i]])
        w = cdms2.MV2.array(wo)
        w.setAxisList(endax)
    # else:
    # w.setAxisList(x.getAxisList())
    return w


def __checker(x, y, w, axes, smally=0):
    # Are the input Variables ?
    xismv = cdms2.isVariable(x)
    yismv = cdms2.isVariable(y)
    if y is None:
        yismv = 1
    wismv = cdms2.isVariable(w)
    if w is None:
        wismv = 1
    ax = None
    if not numpy.ma.isarray(x):
        x = numpy.ma.array(x, copy=0)
    if not numpy.ma.isarray(y) and y is not None:
        y = numpy.ma.array(y, copy=0)
    if not numpy.ma.isarray(w) and w is not None and not isinstance(w, type("")):
        if not isinstance(w[0], type("")):
            w = numpy.ma.array(w, copy=0)
        else:
            if not xismv:
                raise StatisticsError(
                    "Error if weights are a list then x must be an MV2 !!!"
                )
            w = __makeweights(x, w, axes)
            wismv = 1
    elif w is not None:
        if not xismv:
            raise StatisticsError(
                "Error if weights are a list then x must be an MV2 !!!"
            )
        w = __makeweights(x, w, axes)
        wismv = 1

    if xismv * yismv * wismv != 1:
        # We didn't pass all MV2s shapes have to match (unless None)
        if smally == 0:
            if x.shape != numpy.ma.shape(y) and y is not None:
                raise StatisticsError(
                    "Error x and y shape do not match !" + str(x.shape) + "," + str(numpy.ma.shape(y))
                )
        else:
            shy = list(y.shape)
            shy2 = y.shape
            shx = list(x.shape)
            if isinstance(axes, str):
                myaxes = []
                for i in axes:
                    myaxes.append(eval(i))
            elif isinstance(axes, int):
                myaxes = [
                    axes,
                ]
            else:
                myaxes = list(axes)
            for anaxis in myaxes[::-1]:
                shy.insert(0, shx[anaxis])
            y = numpy.ma.resize(y, shy)
            sh = list(range(len(x.shape)))
            if axes != 0:
                for i in range(len(myaxes)):
                    sh[myaxes[i]] = i
                    sh[i] = myaxes[i]
                y = numpy.ma.transpose(y, sh)
            if x.shape != numpy.ma.shape(y) and y is not None:
                err_msg = "Error x and y shape do not match (y shouldbe 1D less than x) !"
                raise StatisticsError(
                    err_msg + str(x.shape) + "," + str(shy2) + " Remember y must be 1D less than x"
                )
        if x.shape != numpy.ma.shape(w) and w is not None:
            msg1 = "Error x and weights shape do not match !"
            msg2 = " ATTENTION if you are trying to pass a list of 1D arrays for each dim, then x must be an MV2 !!!"
            raise StatisticsError(msg1 + str(x.shape) + "," + str(numpy.ma.shape(w)) + msg2)
        if not isinstance(axes, type([])):
            axes = cdms2.orderparse(str(axes))
        for i in axes:
            if len(x.shape) < i:
                err_msg = "Error you have " + str(len(x.shape)) + " dimensions and try to work on dim:" + str(i)
                raise StatisticsError(err_msg)
    else:
        if y is not None:
            x, y = grower(x, y)
            if x.shape != y.shape:
                raise StatisticsError(
                    "Error x and y have different shapes" + str(x.shape) + ", " + str(y.shape)
                )
        ax = x.getAxisList()
        xorder = x.getOrder(ids=1)
        # Now grows w
        if w is not None:
            worder = w.getOrder(ids=1)
            for o in worder:
                if o not in xorder:
                    raise StatisticsError(
                        "Error weights have a dimension that is neither in x or y:" + o
                    )
            x, w = grower(x, w)
            if x.shape != w.shape:
                raise StatisticsError(
                    "Error x and weights have different shapes" + str(x.shape) + ", " + str(w.shape)
                )
        # Last thing convert the axes input to numbers
        if isinstance(axes, type(1)):
            axes = str(axes)
        if not isinstance(axes, type([])):
            axesparse = cdms2.orderparse(axes)
            naxes = len(axesparse)
            for i in range(naxes):
                o = axesparse[i]
                if isinstance(o, type("")):
                    for j in range(len(xorder)):
                        if xorder[j] == o:
                            axesparse[i] = j
                    # Well it must be a name for x y t....
                    if isinstance(axesparse[i], type("")):
                        for j in range(len(x.shape)):
                            if o[1:-1] == x.getAxis(j).id:
                                axesparse[i] = j
                    # Everything failed the axis id must be not existing in the
                    # slab...
                    if isinstance(axesparse[i], type("")):
                        raise StatisticsError(
                            "Error axis id :" + o + " not found in first slab: " + x.getOrder(ids=1)
                        )
            axes = axesparse
    # Now we have array those shape match, and a nice list of axes let's keep
    # going
    naxes = len(axes)
    n0 = 1
    xsh = x.shape
    xorder = list(range(len(x.shape)))
    forder = []
    for i in range(naxes):
        a = axes[i]
        forder.append(a)
        try:
            n0 = n0 * xsh[a]
        except IndexError:
            raise Exception("Axis {} is out of bounds for dimension {}".format(a, len(xsh)))
    fsh = [n0]
    ax2 = []
    for i in range(len(x.shape)):
        if i not in forder:
            forder.append(i)
            fsh.append(xsh[i])
            if ax is not None:
                ax2.append(ax[i])
    if ax is not None:
        ax = ax2
    x = numpy.ma.transpose(x, forder)
    x = numpy.ma.resize(x, fsh)
    if y is not None:
        y = numpy.ma.transpose(y, forder)
        y = numpy.ma.resize(y, fsh)
    if w is not None:
        w = numpy.ma.transpose(w, forder)
        w = numpy.ma.resize(w, fsh)
    # Now mask everything correctly (union of masks)
    if y is not None:
        m = y.mask
        if m is not numpy.ma.nomask:
            x = numpy.ma.masked_where(m, x)
        m = x.mask
        if m is not numpy.ma.nomask:
            y = numpy.ma.masked_where(m, y)
    if w is not None:
        m = x.mask
        if m is not numpy.ma.nomask:
            w = numpy.ma.masked_where(m, w)

    # IF y has to be 1D less than x, then it is shrunk back
    if smally == 1:
        y = y[0]
    return x, y, w, axes, ax

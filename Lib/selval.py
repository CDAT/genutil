# Adapted for numpy/ma/cdms2 by convertcdms.py
from cdms2.selectors import SelectorComponent
import MV2
import numpy
import cdtime
import cdat_info


class PickComponent(SelectorComponent):
    """
    Let the user pick non contiguous values along an axis
    keyword "match" is reserved for handling of non-existing values
    match=1 : (default): raise an exception if one of the select-values does not exist
    match=0 : replace inexistince selcet-values with missing
    match=-1: skip inexisting select-values
    """

    def __init__(self, *args, **kargs):
        """Initialise some values such as tolerances for equality

        :param kargs: match is the most relevant keyword argument. There is a
            match=1 : (default): raise an exception if one of the select-values does not exist
            match=0 : replace inexistence select-values with missing
            match=-1: skip inexistent select-values
        """
        self.args = args
        self.kargs = kargs
        self.match = kargs.get('match', 1)
        if self.match not in [0, 1, -1]:
            raise Exception(
                'Error match must be 1 (strict matching), 0 (missing value) or -1 (skip inexistant values)')

    def __str__(self):
        s = 'Specific non contiguous values Selector\n'
        if len(self.args) > 0:
            s += 'Arguments: ('
            for a in len(self.args):
                s += str(a) + ', '
            s += ')\n'
        if self.kargs != {}:
            s += 'Keywords:\n'
            for k in list(self.kargs.keys()):
                s += '\t' + str(k) + ':' + str(self.kargs[k]) + '\n'
        return s

    def specify(self, slab, axes, specification, confined_by, aux):
        """ First part: confine the slab within a Domain wide enough to do the exact in post"""
        import copy
        from numpy.ma import minimum, maximum
        # myconfined is for later, we can't confine a dimension twice with an
        # argument plus a keyword or 2 keywords
        myconfined = [None] * len(axes)
        self.aux = copy.copy(specification)
        # First look at the arguments (i.e not keywords) and confine the dimensions
        # in the order of the arguments
        for i in range(len(self.args)):
            if confined_by[i] is None:  # Check it hasn't been confined by somebody else
                myconfined[i] = 1  # dim confined by argument list
                # for cdms I want to confine this dimension
                confined_by[i] = self
                # How do we want to confine this dim ?
                self.aux[i] = specs = list(self.args[i])
                if not (isinstance(specs, list) or isinstance(specs, tuple)):
                    raise Exception(
                        "Error in Selector, you must specify a list or a tuple, you passed:" +
                        str(specs))
                elif isinstance(specs[0], type(cdtime.comptime(1999)))\
                        or isinstance(specs[0], type(cdtime.reltime(0, 'days since 1999')))\
                        or isinstance(specs[0], type('')):
                    list2 = []
                    for l in specs:
                        if not isinstance(l, type('')):
                            list2.append(l.torel('days since 1900').value)
                        else:
                            list2.append(
                                cdtime.s2r(
                                    l, 'days since 1900').value)
                    min = minimum(list2)
                    max = maximum(list2)
                    specification[i] = cdtime.reltime(
                        min, 'days since 1900'), cdtime.reltime(
                        max, 'days since 1900')
                else:  # But if it's not...
                    specification[i] = minimum(specs), maximum(
                        specs)  # sets the specifications
            else:
                return 1
        for kw in list(self.kargs.keys()):
            axis = None
            for i in range(len(axes)):
                if axes[i].id == kw:
                    axis = i
            if axis is None:
                if kw == 'time':
                    for i in range(len(axes)):
                        if axes[i].isTime():
                            axis = i
                elif kw == 'level':
                    for i in range(len(axes)):
                        if axes[i].isLevel():
                            axis = i
                elif kw == 'longitude':
                    for i in range(len(axes)):
                        if axes[i].isLongitude():
                            axis = i
                elif kw == 'latitude':
                    for i in range(len(axes)):
                        if axes[i].isLatitude():
                            axis = i
                # keyword not a recognised keyword or dimension name
                elif kw not in ['match']:
                    raise Exception(
                        'Error, keyword: ' + kw + ' not recognized')
            # At this point, if axis is None:
            # we are dealing with a keyword for the selector
            # so we'll skip it
            if axis is not None:
                if confined_by[axis] is None:
                    confined_by[axis] = self
                    myconfined[axis] = 1
                    self.aux[axis] = specs = list(self.kargs[kw])
                    if isinstance(specs[0], type(cdtime.comptime(1999))) or isinstance(
                            specs[0], type(cdtime.reltime(0, 'days since 1999'))) or isinstance(specs[0], type('')):
                        list2 = []
                        for l in specs:
                            if not isinstance(l, type('')):
                                list2.append(l.torel('days since 1900').value)
                            else:
                                list2.append(
                                    cdtime.s2r(
                                        l, 'days since 1900').value)
                        min = minimum(list2)
                        max = maximum(list2)
                        specification[axis] = cdtime.reltime(
                            min, 'days since 1900'), cdtime.reltime(
                            max, 'days since 1900')
                    else:  # But if it's not...
                        specification[axis] = minimum(specs), maximum(specs)

                else:
                    if myconfined[axis] == 1:
                        raise 'Error you are attempting to set the axis: ' + \
                            str(axes[axis].id) + ' more than once'
                    else:
                        return 1
        return 0

    def post(self, fetched, slab, axes, specifications,
             confined_by, aux, axismap):
        """ Post processing retouches the bounds and later will deal with the mask"""
        import cdms2 as cdms
        fetched = cdms.createVariable(fetched, copy=1)
        faxes = fetched.getAxisList()
        a = None
        for i in range(len(faxes)):
            if confined_by[i] is self:
                newaxvals = []
                bounds = []
                a = None
                sh = list(fetched.shape)
                sh[i] = 1
                for l in self.aux[i]:
                    try:
                        tmp = fetched(**{faxes[i].id: (l, l)})
                        ax = tmp.getAxis(i)
                        # print ax
                        newaxvals.append(ax[0])
                        if ax.getBounds() is not None:
                            bounds.append(ax.getBounds()[0])
                        else:
                            bounds = None
                    except BaseException:
                        if self.match == 1:
                            raise Exception(
                                'Error axis value :' +
                                str(l) +
                                ' was requested but is not present in slab\n(more missing might exists)')
                        elif self.match == 0:
                            tmp = MV2.ones(sh, typecode=MV2.float)
                            tmp = MV2.masked_equal(tmp, 1)
                            if isinstance(l, type(cdtime.comptime(1999))) or isinstance(
                                    l, type(cdtime.reltime(0, 'days since 1999'))) or isinstance(l, type('')):
                                if not isinstance(l, type('')):
                                    newaxvals.append(
                                        l.torel(faxes[i].units).value)
                                else:
                                    newaxvals.append(cdtime.s2r(
                                        l, faxes[i].units).value)
                            else:
                                newaxvals.append(l)
                            if bounds is not None:
                                bounds.append([ax[-1] - 1., ax[-1] + 1])
                        else:
                            tmp = None
                    if tmp is not None:
                        if a is None:
                            a = tmp
                        else:
                            a = MV2.concatenate((a, tmp), i)
                if bounds is not None:
                    newax = cdms.createAxis(
                        numpy.array(newaxvals), bounds=numpy.array(bounds), id=ax.id)
                else:
                    newax = cdms.createAxis(numpy.array(newaxvals), id=ax.id)
                for att in list(faxes[i].attributes.keys()):
                    setattr(newax, att, faxes[i].attributes.get(att))
                for j in range(len(fetched.shape)):
                    if j == i:
                        a.setAxis(i, newax)
                    else:
                        a.setAxis(j, faxes[j])
                fetched = a.astype(fetched.dtype.char)
                faxes = fetched.getAxisList()

        return fetched


def picker(*args, **kargs):
    """
    Let the user pick non contiguous values along an axis
    Usage:
    picker(dim2=list1,dim2=list2)
    keyword 'match' is reserved for handling of inexistent values
    match=1 : (default): raise an exception if one of the select-values does not exist
    match=0 : replace inexistince select-values with missing
    match=-1: skip inexisting select-values

    :Example:

        .. doctest:: selval_picker

            # The following examples show numerous ways of specifying a picker
            >>> f=cdms.open('/pcmdi/obs/mo/ta/rnl_ncep/ta.rnl_ncep.ctl')
            >>> s=f('ta',picker(level=[1000,850,700])) # f first levels are 1000.00, 925.00, 850.00, 700.00
            >>> s=f('ta',picker(level=[1000,700,850]) # different order from first example
            >>> s=f('ta',picker(level=[1000,700,800]) # raise an exception since 800 doesn't exist
            >>> s=f('ta',picker(level=[1000,700,800],match=0) # replace 800 level with missing values
            >>> s=f('ta',picker(level=[1000,700,800],match=-1) # skip 800 level
            >>> s=f('ta',genutil.picker(time=['1987-7','1988-3',cdtime.comptime(1989,3)],level=[1000,700,850]))

    """
    #cdat_info.pingPCMDIdb("cdat", "genutil.picker")
    import cdms2 as cdms
    a = cdms.selectors.Selector(PickComponent(*args, **kargs))
    return a

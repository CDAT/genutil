# Adapted for numpy/ma/cdms2 by convertcdms.py
# Adapted for numpy/ma/cdms2 by convertcdms.py
import cdms2 as cdms


def grower(x, y, singleton=0):
    """
    This function takes 2 transient variables and grows them to
    match their axes.
    The final order will be the order of the first variables
    followed by all/any dimension(s) from the second variable not present in the first variable.

    :Example:

        .. doctest:: genutil_grower

            >>> import cdms2, vcs
            >>> vcs.download_sample_data_files()
            >>> f=cdms2.open(vcs.sampledata + '/clt.nc')
            >>> x=f('clt')
            >>> y=f('v')
            >>> x, y = grower(x, y, singleton=1)

    :param x: First transient variable to grow
    :type x: cdms.tvariable.TransientVariable

    :param y: Second transient variable to grow
    :type y: cdms.tvariable.Transientvariable

    :param singleton: Integer flag to indicate whether to raise an error if either dimension is not a singleton.
        Default = 0
        If singletonoption is set to 1 then an error is raised if one of the dims is not a singleton dimension.
    :type singleton: int

    :returns: x and y, adjusted such that their axes match.
    :rtype: cdms.tvariable.TransientVariable
    """
    # Parse the x axes
    xorder = x.getOrder(ids=1)
    xaxes = x.getAxisList()
    # Parse the y axes
    yorder = y.getOrder(ids=1)
    yaxes = y.getAxisList()

    # Now determine the shape of the final array (matching x and y dims,x
    # first)
    forder = []
    prev = 0
    txt = ''
    for o in xorder:
        if o == '(':
            prev = 1
        elif prev == 1:
            if o != ')':
                txt = txt + o
            else:
                forder.append('(%s)' % txt)
                prev = 0
                txt = ''
        else:
            forder.append(o)

    prev = 0
    txt = ''
    xorder = forder[:]
    nyorder = []
    for o in yorder:
        if o == '(':
            prev = 1
        elif prev == 1:
            if o != ')':
                txt = txt + o
            else:
                nyorder.append('(%s)' % txt)
                if '(%s)' % txt not in forder:
                    forder.append('(%s)' % txt)
                prev = 0
                txt = ''
        else:
            nyorder.append(o)
            if o not in forder:
                forder.append(o)
    yorder = nyorder
    # Now grow x
    # print forder,xorder,yorder,nyorder
    for o in forder:
        if o not in xorder:
            for i in range(len(yorder)):
                if yorder[i] == o:
                    newaxes = x.getAxisList()
                    ax = yaxes[i]
                    if len(ax) > 1 and singleton == 1:
                        raise 'Error, dimension:' + ax.id + 'is not a singleton dimension,(len is:' + \
                              str(len(ax)) + \
                              ') you specified to grow only singleton dims, exiting'
                    xsh = list(x.shape)
                    xsh.insert(0, len(ax))
                    x = cdms.MV2.resize(x, xsh)
                    newaxes.insert(0, ax)
                    x.setAxisList(newaxes)
    xorder = x.getOrder(ids=1)
    sp = xorder.split('(')
    xorder = []
    for s in sp:
        if s.find(')') == -1:
            for t in s:
                xorder.append(t)
        else:
            sp2 = s.split(')')
            xorder.append(sp2[0])
            for t in sp2[1]:
                xorder.append(t)
    xaxes = x.getAxisList()

    # Now grow y
    # print forder,yorder
    for o in forder:
        if o not in yorder:
            for i in range(len(xorder)):
                if o in ['(%s)' % xorder[i], xorder[i]]:
                    newaxes = y.getAxisList()
                    ax = xaxes[i]
                    if len(ax) > 1 and singleton == 1:
                        raise 'Error, dimension:' + ax.id + 'is not a singleton dimension,(len is:' + \
                              str(len(ax)) +\
                            ') you specified to grow only singleton dims, exiting'
                    ysh = list(y.shape)
                    ysh.insert(0, len(ax))
                    y = cdms.MV2.resize(y, ysh)
                    newaxes.insert(0, ax)
                    y.setAxisList(newaxes)
    # Figure out the string to reorder x and y
    # print x.shape,y.shape
    return x(order=forder), y(order=forder)

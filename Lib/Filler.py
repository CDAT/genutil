import cdat_info
class StringConstructor:
    """
    This class aims at spotting keywords in a string and replacing them

    :Usage:

        .. code-block:: python

            >>> template = "templates are strings containing any number of %(keywords) using %(this_format)"
            >>> Filler=StringConstructor(template)
            # or
            >>> Filler=StringConstructor()
            >>> Filler.template=template
            # In order to construct the string (i.e. replace keywords with some values):
            >>> str=Filler(keywords='keywdstr',this_format='')
            # or
            >>> Filler.keyword='kwstr'
            >>> Filler.this_format=''
            >>> str=Filler()



    template is a string of form: 'my string here with %(keywords) in it'

    You can have has many keywords as you want, and use them as many times as you want.
    Keywords are delimited on the left by %( and on the right by ).
    """
    def __init__(self,template=None):
        """
        Instantiates a StringConstructor object.
        :param template: A string used by StringConstructor for keyword replacement.
                template is a string of form: 'my string here with %(keywords) in it'.
                There can be an unlimited number of keywords, delimited by %( on the left and ) on the right.
        """
        cdat_info.pingPCMDIdb("cdat","genutil.StringConstructor")
        self.template=template
        ## ok we need to generate the keys and set them to empty it seems like a better idea
        keys = self.keys()
        for k in keys:
            setattr(self,k,"")

    def keys(self,template=None):
        if template is None:
            template=self.template
        if template is None:
            return []
##         # First sets all the keyword values passed
##         for k in kw.keys():
##             setattr(self,k,kw[k])
        # Now determine the keywords in the template:
        end=0
        s2=template.split('%(')
        keys=[]
        for k in s2:
            sp=k.split(')')
            i=len(sp[0])
            if len(k)>i:
                if k[i]==')' and (not sp[0] in  keys):
                    keys.append(sp[0])
        return keys

    def construct(self,template=None,**kw):
        """
        Accepts a string with an unlimited number of keywords to replace.
        Keywords to replace must be in the format %(keyword) within the string.
        Keyword values are either passed as keyword to the construct function or preset.

        :Example:

            .. doctest:: Filler_construct

                >>> structure='/pcmdi/amip/mo/%(variable)/%(model)/%(variable)_%(model).xml'
                >>> Filler=StringConstructor()
                >>> Filler.variable='tas'
                >>> myfilename=Filler.construct(structure,model='ugamp-98a')
                >>> print myfilename
                '/pcmdi/amip/mo/tas/ugamp-98a/tas_ugamp-98a.xml'

        :param template: A string used by StringConstructor for keyword replacement.
                template is a string of form: 'my string here with %(keywords) in it'.
                There can be an unlimited number of keywords, delimited by %( on the left and ) on the right.
        :type template: str

        :param kw: Comma-delimited list of keyword to string value mappings, i.e.:
                    keyword1='kwd1 string',keyword2='kwd2 string', ...
        :type kw: list
        """
        if template is None:
            template=self.template
        # Now determine the keywords in the template:
        keys = self.keys()
        # Now replace the keywords with their values
        for k in keys:
               template=template.replace('%('+k+')',kw.get(k,getattr(self,k,'')))
##             cmd='template=string.replace(template,\'%('+k+')\',self.'+k+')'
##             exec(cmd)
        return template

    def reverse(self,name,debug=False):
        """
        The reverse function attempts to take a template and derive its keyword values based on name parameter.

        :Example:

            .. doctest:: Filler_reverse

                >>> Filler=StringConstructor(template="%(a).%(b)")
                >>> Filler.reverse("A.B")
                {a:"A", b:"B"}

        :param name: String to test the template's keyword values.
        :type name: str

        :param debug: Boolean flag to indicate whether or not to print debug output.
        :type debug: bool

        :returns: A dictionary mapping the StringConstructor's template's keywords to the corresponding values,
            according to the format of the name parameter.
        :rtype: dict

        .. warning::

            reverse makes its best effort at deriving keyword values from a string, but it is not guaranteed to work.
        """
        out={}
        template = self.template
        for k in self.keys():
            sp=template.split("%%(%s)" % k)
            n = len(sp)
            i1=name.find(sp[0])+len(sp[0])
            j1=sp[1].find("%(")
            if j1==-1:
                if sp[1]=="":
                    val=name[i1:]
                else:
                    i2=name.find(sp[1])
                    val = name[i1:i2]
                if debug:
                    print k,j1,sp[1],"****",sp
                    print k,name[i1:i2]
                    print k,i1,i2,val
            else:
                i2=name[i1:].find(sp[1][:j1])
                val=name[i1:i1+i2]
                if debug:
                    print k,j1,sp[1][:j1]
                    print k,name[i1:]
                    print k,i1,i2,val
            if debug:
                print '-----------------'
            template=template.replace("%%(%s)"%k,val)
            if debug:
                print template
            out[k]=val
        if debug:
            print out
        if self.construct(self.template,**out)!=name:
            raise "Invalid pattern sent"
        return out
    
    def __call__(self,*args,**kw):
        """default call is construct function"""
        return self.construct(*args,**kw)


Filler=StringConstructor()

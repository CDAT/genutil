"""
genutil -- General utility modules for scientific computing
"""
# Lean mode does not install xmgrace module
from grower import grower  # noqa
try:
    import xmgrace  # noqa
except BaseException:
    pass
import statistics  # noqa
from minmax import minmax  # noqa
from statusbar import statusbar  # noqa
from selval import picker  # noqa
import filters  # noqa
import arrayindexing  # noqa
import ASCII  # noqa
import udunits  # noqa
from Filler import Filler, StringConstructor  # noqa
from averager import averager, AveragerError, area_weights, getAxisWeight, getAxisWeightByName, __check_weightoptions  # noqa
import cdat_info  # noqa
from ASCII import get_parenthesis_content  # noqa
import os  # noqa
import sys  # noqa

# udunits bits
from udunits import udunits, addBaseUnit, addDimensionlessUnit, addScaledUnit  # noqa
from udunits import addOffsettedUnit, addMultipliedUnits, addInvertedUnit, addDividedUnits  # noqa
udunits_init = 0  # noqa

xml_pth = os.path.join(sys.prefix, "share", "udunits", "udunits2.xml")
if os.path.exists(xml_pth):
    os.environ["UDUNITS2_XML_PATH"] = xml_pth
else:
    try:
        xml_pth = os.path.join(cdat_info.externals, "share", "udunits", "udunits2.xml")
        if os.path.exists(xml_pth):
            os.environ["UDUNITS2_XML_PATH"] = xml_pth
    except BaseException:
        pass

cdat_info.pingPCMDIdb("cdat", "genutil")

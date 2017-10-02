"""
genutil -- General utility modules for scientific computing
"""
# Lean mode does not install xmgrace module
from grower import grower
try:
    import xmgrace
except BaseException:
    pass
import statistics
from minmax import minmax
from statusbar import statusbar
from selval import picker
import filters
#import salstat
import arrayindexing
import ASCII
import udunits
from Filler import Filler, StringConstructor
from averager import averager, AveragerError, area_weights, getAxisWeight, getAxisWeightByName, __check_weightoptions
import cdat_info
from ASCII import get_parenthesis_content
import os
import sys

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

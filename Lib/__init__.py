"""
genutil -- General utility modules for scientific computing
"""
# Lean mode does not install xmgrace module
from .grower import grower  # noqa
try:
    import xmgrace  # noqa
except BaseException:
    pass
from .minmax import minmax  # noqa
from .statusbar import statusbar  # noqa
from .selval import picker  # noqa
from . import filters  # noqa
from . import arrayindexing  # noqa
from . import statistics  # noqa
from . import ASCII  # noqa
from . import udunits  # noqa
from .Filler import Filler, StringConstructor  # noqa
from .averager import averager, AveragerError, area_weights, getAxisWeight, getAxisWeightByName, __check_weightoptions  # noqa
from .ASCII import get_parenthesis_content  # noqa
import os  # noqa
import sys  # noqa
import cdat_info
from .stats_checker import StatisticsError  # noqa

# udunits bits
from .udunits import udunits, addBaseUnit, addDimensionlessUnit, addScaledUnit  # noqa
from .udunits import addOffsettedUnit, addMultipliedUnits, addInvertedUnit, addDividedUnits  # noqa
udunits_init = 0  # noqa

xml_pth = os.path.join(sys.prefix, "share", "udunits", "udunits2.xml")
if os.path.exists(xml_pth):
    os.environ["UDUNITS2_XML_PATH"] = xml_pth
else:
    try:
        xml_pth = os.path.join(
            cdat_info.externals,
            "share",
            "udunits",
            "udunits2.xml")
        if os.path.exists(xml_pth):
            os.environ["UDUNITS2_XML_PATH"] = xml_pth
    except BaseException:
        pass

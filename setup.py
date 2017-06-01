from numpy.distutils.core import setup, Extension
import os,sys,string
import numpy
from os.path import expanduser,expandvars

try:
    sys.path.append(os.environ['BUILD_DIR'])
    import cdat_info
    Version=cdat_info.Version
except:
    Version="???"
#for mars install
expand_sh = lambda path: expanduser(expandvars(path))

#f2py_options = ['--build-src --inplace', '-c', '--verbose', '--build-dir mars/src/']#, '--debug-capi']
MARS = Extension('genutil.mars',
        sources=['Src/mars/src/mars.f', 'Src/mars/src/mars.pyf'],
        #include_dirs = [os.path.join(sys.prefix,'include')],
        #library_dirs = [os.path.join(sys.prefix,'lib')],
        #f2py_options = ['--debug-capi']
            )
setup (name = "genutil",
       version=Version,
       author='LLNL',
       description = "General utilities for scientific computing",
       url = "http://uvcdat.llnl.gov/software",
       packages = ['genutil', 'unidata'], #'genutil.pymars'],
       package_dir = {'genutil': 'Lib', 'unidata':"unidata"}, #"genutil.pymars":"Lib/pymars/src"},
       include_dirs = [numpy.lib.utils.get_include()],
       ext_modules = [
    Extension('genutil.array_indexing',
              ['Src/array_indexing.c',]
              ),
    Extension('genutil.udunits_wrap',
        ['Src/udunits_wrap.c',
            ## 'Src/utparse.c',
            ## 'Src/utlib.c',
            ## 'Src/utscan.c',
            ],
        include_dirs = [os.path.join(sys.prefix,'include')],
        library_dirs = [os.path.join(sys.prefix,'lib')],
        libraries=['udunits2','expat']
        ),
    MARS
    ],
      data_files=[('share/genutil', ('share/test_data_files.txt',))]
      )

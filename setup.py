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
FMARS = Extension('genutil.Fmars',
        sources=['Src/mars/src/mars_nolog.f', 'Src/mars/src/mrsgo1.f', 'Src/mars/src/marsgo.f', 'Src/mars/src/Fmars.pyf'],
            #include_dirs = [os.path.join(sys.prefix,'include')],
            #library_dirs = [os.path.join(sys.prefix,'lib')],
            #f2py_options = ['--build-src --inplace', '-c', '--verbose', '--build-dir mars/src/']#, '--debug-capi']
            )
setup (name = "genutil",
       version=Version,
       author='LLNL',
       description = "General utilities for scientific computing",
       url = "http://uvcdat.llnl.gov/software",
       packages = ['genutil','unidata','genutil.pymars'],
       package_dir = {'genutil': 'Lib', 'unidata':"unidata", "genutil.pymars":"Lib/pymars/src"},
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
    FMARS
    ],
      data_files=[('share/genutil', ('share/test_data_files.txt',))]
      )
yyy=    Extension('mars',
            ['mars_nolog.f', 'mrsgo1.f', 'marsgo.f', 'Fmars.pyf'],
            #include_dirs = [os.path.join(sys.prefix,'include')],
            #library_dirs = [os.path.join(sys.prefix,'lib')],
            f2py_options = ['-c', '--verbose', '--build-src --inplace', '--build-dir mars/src/']#, '--debug-capi']
            )
xxx=    Extension('mars',
            ['src/mars_nolog.f', 'src/mrsgo1.f', 'src/marsgo.f', 'src/Fmars.pyf'],
            #include_dirs = [os.path.join(sys.prefix,'include')],
            #library_dirs = [os.path.join(sys.prefix,'lib')],
            f2py_options = ['-c', '--verbose', '--build-src --inplace', '--build-dir mars/src/']#, '--debug-capi']
            )

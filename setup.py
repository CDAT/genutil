from distutils.core import setup, Extension
import os,sys,string
import numpy
try:
    sys.path.append(os.environ['BUILD_DIR'])
    import cdat_info
    Version=cdat_info.Version
except:
    Version="???"
setup (name = "genutil",
       version=Version,
       author='PCMDI',
       description = "General utilities for scientific computing",
       url = "http://www-pcmdi.llnl.gov/software",
       packages = ['genutil'],
       package_dir = {'genutil': 'Lib'},
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
        )

    
    ],
      data_files=[('share/genutil', ('share/test_data_files.txt',))]
      )


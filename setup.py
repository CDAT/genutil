from distutils.core import setup, Extension
import os,sys
import numpy
try:
    import cdat_info
    Version=cdat_info.Version
except:
    Version="???"
print("VERSION:",Version)
setup (name = "genutil",
       version=Version,
       author='LLNL',
       description = "General utilities for scientific computing",
       url = "http://uvcdat.llnl.gov/software",
       packages = ['genutil','unidata'],
       package_dir = {'genutil': 'Lib', 'unidata':"unidata"},
       include_dirs = ['Include', numpy.lib.utils.get_include()],
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

#!/usr/bin/env bash
ls
pwd
export PATH=${HOME}/miniconda/bin:${PATH}
conda create -n py3 -c cdat/label/nightly -c uvcdat/label/nightly -c nesii/label/dev-esmf -c conda-forge -c uvcdat cdms2 nose flake8 "python>3" cdat_info udunits2
conda create -n py2 -c cdat/label/nightly -c uvcdat/label/nightly -c conda-forge -c uvcdat cdms2 cdat_info udunits2 nose flake8
export UVCDAT_ANONYMOUS_LOG=False
source activate py3
python setup.py install
source activate py2
rm -rf build
python setup.py install

#!/usr/bin/env bash
ls
pwd
export PATH=${HOME}/miniconda/bin:${PATH}
conda install -c uvcdat/label/nightly -c conda-forge -c uvcdat cdms2 cdat_info udunits2 nose
pip install dropbox
export UVCDAT_ANONYMOUS_LOG=False
python setup.py install

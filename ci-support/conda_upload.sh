#!/usr/bin/env bash
PKG_NAME=genutil
USER=cdat
export PATH="$HOME/miniconda/bin:$PATH"
echo "Trying to upload conda"
ESMF_CHANNEL="nesii/label/dev-esmf"
if [ `uname` == "Linux" ]; then
    OS=linux-64
    echo "Linux OS"
    conda update -y -q conda
else
    echo "Mac OS"
    OS=osx-64
fi

mkdir ~/conda-bld
source activate root
conda install -q anaconda-client conda-build
conda config --set anaconda_upload no
export CONDA_BLD_PATH=${HOME}/conda-bld
#export VERSION=`date +%Y.%m.%d`
echo "Cloning recipes"
git clone git://github.com/UV-CDAT/conda-recipes
cd conda-recipes
# uvcdat creates issues for build -c uvcdat confises package and channel
rm -rf uvcdat
python ./prep_for_build.py
echo "Building and uploading now"
conda build -c ${ESMF_CHANNEL} -c conda-forge -c cdat/label/nightly -c uvcdat/label/nightly -c uvcdat ${PKG_NAME}
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l nightly $CONDA_BLD_PATH/$OS/$PKG_NAME*.tar.bz2 --force

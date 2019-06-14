#!/usr/bin/env bash

# Build script for BioSimSpace Conda installation.

# First install the Python package.
cd python && BSS_CONDA_INSTALL=True python setup.py install --single-version-externally-managed --record=record.txt

# Install additional pip packages. Note that conda-build turns off downloading
# so pip won't be able to get packages from PyPi. We re-enable downloading by
# setting PIP_NO_INDEX to True, i.e. 0.
PIP_NO_INDEX=0 pip install --index-url https://pypi.org/simple --no-deps --ignore-installed --verbose fileupload pygtail pypdb

# Install and enable notebook extensions.
jupyter-nbextension install fileupload --py --sys-prefix
jupyter-nbextension enable fileupload --py --sys-prefix
jupyter-nbextension install nglview --py --sys-prefix
jupyter-nbextension enable nglview --py --sys-prefix

# Now download and unpack KCOMBU.
wget --user=kcombu --password=wakame -r --accept .tar.gz --level 1 --cut-dirs 3 -nH http://strcomp.protein.osaka-u.ac.jp/kcombu/src/
mkdir kcombu
tar -xzf kcombu*.tar.gz -C kcombu
cd kcombu/src

# Update the C compiler on Linux.
if [ "$(uname)" == "Linux" ]; then
    sed -i "s/gcc/$CC/g" Makefile.fkcombu
fi

# Build FKCOMBU.
make -f Makefile.fkcombu

# Copy the FKCOMBU executable to the bin directory.
cp -a ../fkcombu ${PREFIX}/bin

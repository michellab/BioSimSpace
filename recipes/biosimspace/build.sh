#!/usr/bin/env bash

# Build script for BioSimSpace Conda installation.

# First install the Python package.
cd python && BSS_CONDA_INSTALL=True python setup.py install --single-version-externally-managed --record=record.txt

# Install and enable notebook extensions.
jupyter-nbextension install nglview --py --sys-prefix
jupyter-nbextension enable nglview --py --sys-prefix

# Now download and unpack KCOMBU.
curl -O https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/XG5XPIRcsIZsxzwzbrLKuTfJ5Ow6r4AU0YxbU6_Lnl31BrwTFlvToXOOjpbheK2R/n/hugs/b/notebook/o/kcombu-src-20200414.tar.gz
mkdir kcombu
tar -xzf kcombu*.tar.gz -C kcombu
cd kcombu/src

if [ "$(uname)" == "Linux" ]; then
    # Update the C compiler on Linux.
    sed -i "s#gcc#$CC#g" Makefile.fkcombu
else
    # Ignore implicit function declaration errors on macOS.
    sed -i.bak -e "s#Wall#Wno-implicit-function-declaration#g" Makefile.fkcombu
fi

# Build FKCOMBU.
make -f Makefile.fkcombu

# Copy the FKCOMBU executable to the bin directory.
cp -a ../fkcombu ${PREFIX}/bin

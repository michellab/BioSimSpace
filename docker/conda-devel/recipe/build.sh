#!/usr/bin/env bash

# Build script for BioSimSpace Conda installation.

# First install the Python package.
cd python && BSS_CONDA_INSTALL=True python setup.py install --single-version-externally-managed --record=record.txt

# Now download and build fkcombu.
wget --user=kcombu --password=wakame -r --accept .tar.gz --level 1 --cut-dirs 3 -nH http://strcomp.protein.osaka-u.ac.jp/kcombu/src/
mkdir kcombu
tar -xzf kcombu*.tar.gz -C kcombu
cd kcombu/src
make -f Makefile.fkcombu

# Copy the fkcombu executable to the bin directory.
cp -a ../fkcombu ${PREFIX}/bin

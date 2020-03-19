#!/usr/bin/env bash

# Build script for BioSimSpace Conda installation.

# First install the Python package.
cd python && BSS_CONDA_INSTALL=True python setup.py install --single-version-externally-managed --record=record.txt

# Install and enable notebook extensions.
jupyter-nbextension install nglview --py --sys-prefix
jupyter-nbextension enable nglview --py --sys-prefix

# Now download and unpack KCOMBU.
# By default, we download the latest src from http://strcomp.protein.osaka-u.ac.jp
# However, we periodically upload a version to the Oracle Cloud to use as a
# fallback in case the server is down.
wget --user=kcombu --password=wakame -r --accept .tar.gz --level 1 --cut-dirs 3 -nH http://strcomp.protein.osaka-u.ac.jp/kcombu/src/
if ! ls kcombu*.tar.gz > /dev/null 2>&1; then
    curl -O https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/7p17kkk6oW8MI60-tlyQmWjOm3kzsZMzbYKgSd_gWto/n/hugs/b/notebook/o/kcombu-src-20190208.tar.gz
fi
mkdir kcombu
tar -xzf kcombu*.tar.gz -C kcombu
cd kcombu/src

# Update the C compiler on Linux.
if [ "$(uname)" == "Linux" ]; then
    sed -i "s#gcc#$CC#g" Makefile.fkcombu
fi

# Build FKCOMBU.
make -f Makefile.fkcombu

# Copy the FKCOMBU executable to the bin directory.
cp -a ../fkcombu ${PREFIX}/bin

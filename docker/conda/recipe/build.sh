#!/usr/bin/env bash

# Build script for BioSimSpace Conda installation.

# First install the Python package.
cd python && BSS_CONDA_INSTALL=True python setup.py install --single-version-externally-managed --record=record.txt

# Install and enable notebook extensions.
jupyter-nbextension install nglview --py --sys-prefix
jupyter-nbextension enable nglview --py --sys-prefix

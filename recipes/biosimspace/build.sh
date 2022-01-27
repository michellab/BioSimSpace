#!/usr/bin/env bash

# Build script for BioSimSpace Conda installation.

cd python && BSS_CONDA_INSTALL=True python setup.py install --single-version-externally-managed --record=record.txt

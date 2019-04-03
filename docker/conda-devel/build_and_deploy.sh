#!/usr/bin/env bash

# Get the Anaconda Cloud access token.
ANACONDA_TOKEN=${ANACONDA_TOKEN:-$1}

# Check the OS.
if [ "$(uname)" == "Darwin" ]; then
    OS=osx-64
else
    OS=linux-64
fi

# Set the bin directory.
BIN_DIR=$HOME/sire.app/bin

# Set the Conda build directory on macOS.
CONDA_DIR=./docker/conda-devel/recipe

# Linux runs in a docker container from $HOME.
if [ ! -d $CONDA_DIR ]; then
    CONDA_DIR=$HOME/BioSimSpace/docker/conda-devel/recipe
fi

# Move the to build directory.
cd $CONDA_DIR

# Build the Conda package.
$BIN_DIR/conda-build -c conda-forge -c omnia -c michellab .

# Upload the package to the michellab channel on Anaconda Cloud.
$BIN_DIR/anaconda -t $ANACONDA_TOKEN upload --user michellab $HOME/sire.app/conda-bld/$OS/biosimspace-* --label dev

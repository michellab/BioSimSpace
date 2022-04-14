#!/usr/bin/env bash

# Get the path to source directory.
SRC_DIR=${SRC_DIR:-$1}

# Get the path to the Conda installation.
CONDA=${CONDA:-$2}

# Get the Anaconda Cloud access token.
ANACONDA_TOKEN=${ANACONDA_TOKEN:-$3}

# Set the path to the conda-bld directory.
CONDA_BLD=${CONDA}/envs/bss_build/conda-bld

# Check the OS.
if [ "$(uname)" == "Darwin" ]; then
    OS=osx-64
else
    OS=linux-64
fi

# Get the tag associated with the latest commit.
TAG=$(git --git-dir="$SRC_DIR"/.git --work-tree="$SRC_DIR" tag --contains)

# If the tag is not empty, then set the label to main.
if [ ! -z "$TAG" ]; then
    LABEL=main
fi

# Upload the package to the michellab channel on Anaconda Cloud.

# Label release packages with main and dev so that dev is at least as new as main.
if [ "$LABEL" = "main" ]; then
    anaconda \
        --token "$ANACONDA_TOKEN" upload \
        --user michellab \
        --label main \
        --label dev \
        --force \
        "$CONDA_BLD"/"$OS"/biosimspace-*
else
    anaconda \
        --token "$ANACONDA_TOKEN" upload \
        --user michellab \
        --label dev \
        --force \
        "$CONDA_BLD"/"$OS"/biosimspace-*
fi

echo "Package uploaded!"

#!/usr/bin/env bash

# Set the Conda build directory on macOS.
CONDA_DIR=./docker/conda-devel/recipe

# Linux runs in a docker container from $HOME.
if [ ! -d $CONDA_DIR ]; then
    CONDA_DIR=$HOME/BioSimSpace/docker/conda-devel/recipe
fi

# Set the BioSimSpace source directory.
SRC_DIR=.

# Linux runs in a docker container from $HOME.
if [ ! -d $SRC_DIR ]; then
    SRC_DIR=$HOME/BioSimSpace
fi

# Store the name of the recipe and template YAML files.
RECIPE=$CONDA_DIR/meta.yaml
TEMPLATE=$CONDA_DIR/template.yaml

# Overwite the recipe with the template file.
cp $TEMPLATE $RECIPE

# Get the BioSimSpace version.
BSS_VER=$(git --git-dir=$SRC_DIR/.git describe --tags | tr - _)

# Update the BioSimSpace version number.
echo "Updating BioSimSpace version number: '$BSS_VER'"
sed -i.bak "s/VERSION/$BSS_VER/" -- $RECIPE && rm -- $RECIPE.bak

echo "Recipe updated!"

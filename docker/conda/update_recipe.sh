#!/usr/bin/env bash

# Set the source and Conda build directories on macOS.
SRC_DIR=$(pwd)
CONDA_DIR=$SRC_DIR/docker/conda/recipe

# Linux runs in a docker container from $HOME.
if [ ! -d $CONDA_DIR ]; then
    SRC_DIR=$HOME/BioSimSpace
    CONDA_DIR=$HOME/BioSimSpace/docker/conda/recipe
fi

# Store the name of the recipe and template YAML files.
RECIPE=$CONDA_DIR/meta.yaml
TEMPLATE=$CONDA_DIR/template.yaml

# Overwite the recipe with the template file.
cp $TEMPLATE $RECIPE

# Get the BioSimSpace version.
BSS_VER=$(git --git-dir=$SRC_DIR/.git --work-tree=$SRC_DIR describe --tags --abbrev=0)

# Get the build number. (Number of commits since last tag.)
BSS_BUILD=$(git --git-dir=$SRC_DIR/.git --work-tree=$SRC_DIR log --oneline $BSS_VER.. | wc -l)

# Update the BioSimSpace version number.
echo "Updating BioSimSpace version number: '$BSS_VER'"
sed -i.bak -e "s/VERSION/$BSS_VER/" $RECIPE && rm $RECIPE.bak

# Update the build number.
echo "Updating BioSimSpace build number: '$BSS_BUILD'"
sed -i.bak -e "s/BUILD/$BSS_BUILD/" $RECIPE && rm $RECIPE.bak

echo "Recipe updated!"

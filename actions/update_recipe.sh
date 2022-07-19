#!/usr/bin/env bash

# Get the path to the source directories.
BSS_SRC_DIR=${BSS_SRC_DIR:-$1}
SIRE_SRC_DIR=${SIRE_SRC_DIR:-$1}

# Set the Conda build directory path.
CONDA_DIR="$BSS_SRC_DIR"/recipes/biosimspace

# Store the name of the recipe and template YAML files.
RECIPE="$CONDA_DIR"/meta.yaml
TEMPLATE="$CONDA_DIR"/template.yaml

# Overwite the recipe with the template file.
cp "$TEMPLATE" "$RECIPE"

# Get the Sire version. (Latest tag.)
#SIRE_VERSION=$(git --git-dir="$SIRE_SRC_DIR"/.git --work-tree="$SIRE_SRC_DIR" describe --tags --abbrev=0)
SIRE_VERSION=2022.3.0

# Get the next major Sire version number.
NEXT_SIRE_VERSION=$(echo $SIRE_VERSION | cut -d '.' -f1)
NEXT_SIRE_VERSION=$((NEXT_SIRE_VERSION+1))

# Get the BioSimSpace version. (Latest tag.)
BSS_VERSION=$(git --git-dir="$BSS_SRC_DIR"/.git --work-tree="$BSS_SRC_DIR" describe --tags --abbrev=0)

# Get the build number. (Number of commits since last tag.)
BSS_BUILD=$(git --git-dir="$BSS_SRC_DIR"/.git --work-tree="$BSS_SRC_DIR" log --oneline "$BSS_VERSION".. | wc -l | xargs)

# Get the BioSimSpace branch.
BSS_BRANCH=$(git --git-dir="$BSS_SRC_DIR"/.git --work-tree="$BSS_SRC_DIR" rev-parse --abbrev-ref HEAD)

# Update the BioSimSpace version number.
echo "Updating BioSimSpace version number: '$BSS_VERSION'"
sed -i.bak -e "s/BSS_VERSION/$BSS_VERSION/" "$RECIPE" && rm "$RECIPE".bak

# Update the build number.
echo "Updating BioSimSpace build number: '$BSS_BUILD'"
sed -i.bak -e "s/BSS_BUILD/$BSS_BUILD/" "$RECIPE" && rm "$RECIPE".bak

# Update the branch name.
echo "Updating BioSimSpace branch name: '$BSS_BRANCH'"
sed -i.bak -e "s/BSS_BRANCH/$BSS_BRANCH/" "$RECIPE" && rm "$RECIPE".bak

# Update the Sire dependency version.
echo "Updating Sire dependency version: '>=$SIRE_VERSION,<$NEXT_SIRE_VERSION.0a0'"
sed -i.bak -e "s/NEXT_SIRE_VERSION/$NEXT_SIRE_VERSION/" "$RECIPE" && rm "$RECIPE".bak
sed -i.bak -e "s/SIRE_VERSION/$SIRE_VERSION/" "$RECIPE" && rm "$RECIPE".bak

echo "Recipe updated!"

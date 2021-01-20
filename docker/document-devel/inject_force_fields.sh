#!/usr/bin/env bash

# BioSimSpace dynamically generates function names for force fields from the
# Open Force Field initiative. As such, we can't hard-code these into the
# Sphinx docstrings within the BioSimSpace.Parameters __init__ file. Here
# we use a Python script to capture the names of the functions, then inject
# them into the __init__ file using sed.

# Query the force fields that are supported in this build and capture the output.
force_fields=$($HOME/sire.app/bin/python $HOME/BioSimSpace/docker/document-devel/supported_forcefields.py)

# Inject the string into the __init__ of the BioSimSpace.Parameters package.
# This is only run on Linux so the following, non-portable, sed command is fine.
sed -i "s/__FORCE_FIELDS__/$forcefields/g" $HOME/sire.app/lib/python3.7/site-packages/BioSimSpace/Parameters/__init__.py

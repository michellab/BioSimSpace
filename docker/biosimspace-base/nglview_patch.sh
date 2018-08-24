#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "usage: ./nglview_patch /path/to/sire.app"
    exit
fi

# Get the location of sire.app.
SIRE=$1

# Check that the directory exists.
if [ ! -d "$SIRE" ]; then
    echo "$SIRE doesn't exist!"
    exit
fi

# Set the name of the problem file.
FIX_FILE="$SIRE/lib/python3.5/site-packages/nglview/utils/py_utils.py"

# Check that it exists.
if [ ! -f "$FIX_FILE" ]; then
    echo "$FIX_FILE doesn't exist!"
    exit
fi

# Patch the file, accounting for different sed behaviour on macOS and FreeBSD
if [ "$os_name" = "Darwin" ] || [ "$os_name" = "FreeBSD" ]; then
    sed -i "" "s#'rb'#'r'#g" $FIX_FILE
else
    sed -i "s#'rb'#'r'#g" $FIX_FILE
fi

echo "File patched!"

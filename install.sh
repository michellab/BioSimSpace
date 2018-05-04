#!/usr/bin/env bash

# Installer script for BioSimSpace.

# Store the current directory.
CURR_DIR=$(pwd)

# Set the default branch.
BRANCH=master

# Store the OS name.
OS="$(uname)"

# Only Linux and macOS are supported.
if [[ ! ($OS == "Linux" || $OS == "Darwin") ]]; then
    echo "Installer only works on Linux and macOS!"
    exit 1
fi

if [ $# -gt 0 ]; then
    if [[ "$1" == "master" ]]; then
        BRANCH=master
    elif [[ "$1" == "devel" ]]; then
        BRANCH=devel
    else
        echo "Unsupported branch! Defaulting to 'master'"
    fi
fi

if [ -z "$INSTALL_DIR" ]; then
    # Ask the user where they would like to install BioSimSpace. By default
    # we will aim for $HOME/BioSimSpace.app
    echo -n "Where would you like to install BioSimSpace? [$HOME/BioSimSpace.app]: "
    read INSTALL_DIR

    if [ ! ${INSTALL_DIR} ]; then
        INSTALL_DIR=$HOME/BioSimSpace.app
    else
        # Use eval so that we can expand variables such as $HOME
        INSTALL_DIR=`eval echo ${INSTALL_DIR}`
    fi
fi

echo "Installing into directory: ${INSTALL_DIR}"

# Create the installation directory.
mkdir -p ${INSTALL_DIR}

echo "  --> Downloading Sire binary"

# Download latest self-extracting Sire binary.
if [[ $OS == "Linux" ]]; then
    curl -sL http://siremol.org/largefiles/sire_releases/download.php?name=sire_2018_1_1_linux.run -o ${INSTALL_DIR}/sire.run
else
    curl -sL http://siremol.org/largefiles/sire_releases/download.php?name=sire_2018_1_1_osx.run -o ${INSTALL_DIR}/sire.run
    mkdir -p $HOME/.matplotlib
    touch $HOME/.matplotlib/matplotlibrc
    if ! grep -q "backend" $HOME/.matplotlib/matplotlibrc; then
        echo "backend: TkAgg" >> $HOME/.matplotlib/matplotlibrc
    else
        sed -i '' 's/.*backend.*/backend: TkAgg/' $HOME/.matplotlib/matplotlibrc
    fi
fi

# Make the binary executable.
chmod a+x ${INSTALL_DIR}/sire.run

# Store the current DISPLAY variable. This allows us to unset DISPLAY to ensure
# that Sire binary unpacking runs in the current shell.
CURR_DISPLAY=$DISPLAY

echo "  --> Installing Sire application"

# Unpack the binary.
unset DISPLAY && echo "${INSTALL_DIR}/sire.app" | ${INSTALL_DIR}/sire.run > /dev/null 2>&1

# Reset the DISPLAY.
export DISPLAY=$CURR_DISPLAY

# Remove the binary file.
rm ${INSTALL_DIR}/sire.run

echo "  --> Downloading BioSimSpace"

# Download the latest zip archives from the BioSimSpace master branches (source and tests).
curl -sL https://github.com/michellab/BioSimSpace/archive/$BRANCH.zip -o ${INSTALL_DIR}/$BRANCH.zip
curl -sL https://github.com/michellab/BioSimSpaceUnitTests/archive/master.zip -o ${INSTALL_DIR}/tests.zip

# Change to the installation directory.
cd ${INSTALL_DIR}

# Remove old demo directory.
if [ -d demo ]; then
    rm -r demo
fi

# Remove old tests directory.
if [ -d tests ]; then
    rm -r tests
fi

# Unzip the BioSimSpace archives.
unzip -q $BRANCH.zip
unzip -q tests.zip
rm $BRANCH.zip
rm tests.zip

# Change to the python directory.
cd BioSimSpace-$BRANCH/python

echo "  --> Installing BioSimSpace"
${INSTALL_DIR}/sire.app/bin/python setup.py install > /dev/null 2>&1

# Clean up.
cd ${INSTALL_DIR}
mv BioSimSpace-$BRANCH/demo .
rm -r BioSimSpace-$BRANCH
mv BioSimSpaceUnitTests-master/tests .
rm -r BioSimSpaceUnitTests-master

# Switch back to original workspace.
cd $CURR_DIR

# Add aliases to ~/.biosimspacerc
echo "# BioSimSpace aliases." > $HOME/.biosimspacerc
echo "alias bss_python=${INSTALL_DIR}/sire.app/bin/python" >> $HOME/.biosimspacerc
echo "alias bss_ipython=${INSTALL_DIR}/sire.app/bin/ipython" >> $HOME/.biosimspacerc
echo "alias bss_jupyter=${INSTALL_DIR}/sire.app/bin/jupyter" >> $HOME/.biosimspacerc
echo "alias bss_test='cd ${INSTALL_DIR}; sire.app/bin/pytest -v tests; cd -'" >> $HOME/.biosimspacerc

# Store the name of the shell rc file.
SHELL_RC=$HOME/.$(basename "$SHELL")rc

# Source aliases in the SHELL_RC file
if ! grep -q "BioSimSpace" $SHELL_RC; then
    echo -e "\n# Source BioSimSpace aliases." >> $SHELL_RC
    echo "source $HOME/.biosimspacerc" >> $SHELL_RC
fi

echo "Installation complete!"
echo "You'll need to run 'source ~/.biosimspacerc' or restart your shell for aliases to take effect."

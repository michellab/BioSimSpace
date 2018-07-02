#!/usr/bin/env bash

# Installer script for BioSimSpace.

# Check required command-line utilities are installed.

# Unzip is needed for uncompressing archive files from GitHub.
if ! [ -x "$(command -v unzip)" ]; then
    echo "Error: 'unzip' is not installed." >&2
    exit 1
fi

# Curl or wget are needed to download Sire and BioSimSpace.
if [ -x "$(command -v curl)" ]; then
    HAS_CURL=true
elif [ -x "$(command -v wget)" ]; then
    HAS_CURL=false
else
    echo "Error: 'curl' or 'wget' must be installed." >&2
    exit 1
fi

# Download function definition.
download() {
    if [[ "$HAS_CURL" == "true" ]]; then
        if ! curl --silent --insecure --location "$1" --output "$2"; then
            echo "Failed to download: $1" >&2
            exit 1
        fi
    else
        if ! wget --quiet --no-check-certificate "$1" --output-document "$2"; then
            echo "Failed to download: $1" >&2
            exit 1
        fi
    fi
}

# Current Sire version.
SIRE_VER=201811

# Store the current directory.
CURR_DIR=$(pwd)

# Set the default branch.
BRANCH=master

# Store the OS name.
OS="$(uname)"

# Only Linux and macOS are supported.
if [[ ! ($OS == "Linux" || $OS == "Darwin") ]]; then
    echo "Installer only works on Linux and macOS!" >&2
    exit 1
fi

# Simple parser for command-line options.
for i in "$@"; do
    case "$1" in
    -b=*|--branch=*)
        BRANCH="${i#*=}";;
    -c|--clean)
        CLEAN_BUILD=true;;
    -*)
        echo " Unrecognized option \"$1\"" >&2;
        echo " Available options are: 'branch' and 'clean'" >&2
        exit 1;;
    esac
    shift
done

# Check branch is valid.
if [[ ! ($BRANCH == master || $BRANCH == devel) ]]; then
	echo " Unsupported branch: $BRANCH" >&2
	echo " Available options are: 'master' or 'devel'" >&2
	exit 1
fi

if [ -z "$INSTALL_DIR" ]; then
    # Ask the user where they would like to install BioSimSpace. By default
    # we will aim for $HOME/BioSimSpace.app
    echo -n "Where would you like to install BioSimSpace? [$HOME/BioSimSpace.app]: "
    read -r INSTALL_DIR

    if [ ! "$INSTALL_DIR" ]; then
        INSTALL_DIR="$HOME/BioSimSpace.app"
    else
        # Use eval so that we can expand variables such as $HOME
        INSTALL_DIR=$(eval echo "$INSTALL_DIR")
    fi
fi

echo "Installing into directory: $INSTALL_DIR"

# Store the location of the Sire version file.
SIRE_VER_FILE="$INSTALL_DIR/SIRE_VER.txt"

# Set the Sire application directory.
SIRE_DIR="$INSTALL_DIR/sire.app"

# Create the installation directory.
mkdir -p "$INSTALL_DIR"

# Check for existing Sire installation.
if [ -z "$CLEAN_BUILD" ]; then
    if [ -f "$SIRE_VER_FILE" ]; then
        # Get the current version.
        CURR_VER=$(cat "$SIRE_VER_FILE")

        # Is the version up to date?
        if [ "$SIRE_VER" -gt "$CURR_VER" ]; then
            CLEAN_BUILD=true
        fi
    else
        CLEAN_BUILD=true
    fi
fi

# Download and install Sire.
if ! [ -z "$CLEAN_BUILD" ]; then

    echo "  --> Downloading Sire binary"

    # Clean any existing Sire installation.
    if [ -d "$SIRE_DIR" ]; then
        rm -r "$SIRE_DIR"
    fi

    # Download latest self-extracting Sire binary.
    if [[ $OS == "Linux" ]]; then
    download http://siremol.org/largefiles/sire_releases/download.php?name=sire_2018_1_1_linux.run "$INSTALL_DIR/sire.run"
    else
        download http://siremol.org/largefiles/sire_releases/download.php?name=sire_2018_1_1_osx.run "$INSTALL_DIR/sire.run"
        mkdir -p "$HOME/.matplotlib"
        touch "$HOME/.matplotlib/matplotlibrc"
        if ! grep -q "backend" "$HOME/.matplotlib/matplotlibrc"; then
            echo "backend: TkAgg" >> "$HOME/.matplotlib/matplotlibrc"
        else
            sed -i '' 's/.*backend.*/backend: TkAgg/' "$HOME/.matplotlib/matplotlibrc"
        fi
    fi

    # Make the binary executable.
    chmod a+x "$INSTALL_DIR/sire.run"

    # Store the current DISPLAY variable. This allows us to unset DISPLAY to ensure
    # that Sire binary unpacking runs in the current shell.
    CURR_DISPLAY=$DISPLAY

    echo "  --> Installing Sire application"

    # Unpack the binary.
    unset DISPLAY && echo "$SIRE_DIR" | "$INSTALL_DIR/sire.run" > /dev/null 2>&1 && echo $SIRE_VER > "$SIRE_VER_FILE"

    # Reset the DISPLAY.
    export DISPLAY=$CURR_DISPLAY

    # Remove the binary file.
    rm "$INSTALL_DIR"/sire.run
fi

echo "  --> Downloading BioSimSpace"

# Download the latest zip archives from the BioSimSpace master branches (source and tests).
download "https://github.com/michellab/BioSimSpace/archive/$BRANCH.zip" "$INSTALL_DIR/$BRANCH.zip"
download https://github.com/michellab/BioSimSpaceUnitTests/archive/master.zip "$INSTALL_DIR/tests.zip"

# Change to the installation directory.
cd "$INSTALL_DIR" || exit 1

# Remove old demo directory.
if [ -d demo ]; then
    rm -r demo
fi

# Remove old tests directory.
if [ -d tests ]; then
    rm -r tests
fi

# Unzip the BioSimSpace archives.
unzip -q "$BRANCH.zip"
unzip -q tests.zip
rm "$BRANCH.zip"
rm tests.zip

# Change to the python directory.
cd "BioSimSpace-$BRANCH/python" || exit 1

echo "  --> Installing BioSimSpace"
"$SIRE_DIR/bin/python" setup.py install > /dev/null 2>&1

# Clean up.
cd "$INSTALL_DIR" || exit 1
mv "BioSimSpace-$BRANCH/demo" .
rm -r "BioSimSpace-$BRANCH"
mv BioSimSpaceUnitTests-master/tests .
rm -r BioSimSpaceUnitTests-master

# Switch back to original workspace.
cd "$CURR_DIR" || exit 1

# Add aliases to ~/.biosimspacerc
echo "# BioSimSpace aliases." > "$HOME/.biosimspacerc"
{
    echo "alias bss_python='$SIRE_DIR/bin/python'"
    echo "alias bss_ipython='$SIRE_DIR/bin/ipython'"
    echo "alias bss_jupyter='$SIRE_DIR/bin/jupyter'"
    echo "alias bss_test='cd $SIRE_DIR/bin/pytest -v tests; cd -'"
} >> "$HOME/.biosimspacerc"

# Store the name of the shell rc file.
SHELL_RC="$HOME/.$(basename "$SHELL")rc"

# Source aliases in the SHELL_RC file
if ! grep -q "BioSimSpace" "$SHELL_RC"; then
    printf "\n# Source BioSimSpace aliases.\n" >> "$SHELL_RC"
    echo "source $HOME/.biosimspacerc" >> "$SHELL_RC"
fi

echo "Installation complete!"
echo "You'll need to run 'source ~/.biosimspacerc' or restart your shell for aliases to take effect."

echo
echo "For additional support, please install the following packages..."
echo "AMBER:   http://ambermd.org"
echo "GROMACS: http://www.gromacs.org"
echo "NAMD:    http://www.ks.uiuc.edu/Research/namd"

#!/bin/sh

set install_dir=$BIOSIMSPACE_INSTALL_DIR

if [ "x$install_dir" = "x" ]; then

    echo " "
    echo "#######################################################"
    echo "#      Where would you like to install BioSimSpace?"
    echo "#      Note that you can move BioSimSpace to any"
    echo "#      directory you want after installation"
    echo "#######################################################\n"
    echo " "

    echo "Full Path to Install directory ($HOME/biosimspace.app): "

    read install_dir

    if [ "x$install_dir" = "x" ]; then
        echo "Using default install directory: $HOME/biosimspace.app"
        install_dir="$HOME/biosimspace.app"
    fi
else
    echo " "
    echo "#######################################################"
    echo "#      Installing to $install_dir"
    echo "#      If you would like to change this, then unset"
    echo "#      the BIOSIMSPACE_INSTALL_DIR environment variable."
    echo "#######################################################\n"
fi

echo "Going to install BioSimSpace to $install_dir"

if [ -d "$install_dir" ]; then
    echo "There is already a version of BioSimSpace installed at $install_dir"
    echo "Please remove it, or choose another installation directory"
    exit -1
fi

echo " "
echo "Installing BioSimSpace to $install_dir"
mv tmp_sire.app $install_dir

if [ ! -d "$install_dir" ]; then
  echo " "
  echo "************************************************"
  echo "* WARNING - INSTALLATION FAILED"
  echo "* PLEASE CHECK THAT YOU CAN WRITE TO $install_dir"
  echo "* IF YOU CAN, PLEASE CONTACT THE DEVELOPERS AT"
  echo "* http://biosimspace.org"
  echo "************************************************"
  echo " "
  exit -1
fi

$install_dir/bin/python $install_dir/pkgs/sire-*/share/Sire/build/restore_path.py $install_dir

echo " "
echo "##############################################################################"
echo "##"
echo "## CONGRATULATIONS - SUCCESSFUL INSTALLATION"
echo "##"
echo "## BioSimSpace is installed in $install_dir"
echo "## You can run a BioSimSpace python script called script.py by typing"
echo "## $install_dir/bin/python script.py"
echo "##"
echo "## All BioSimSpace binaries are available in "
echo "## $install_dir/bin"
echo "##"
echo "## Everything is contained in this directory, so you can delete BioSimSpace"
echo "## by deleting this directory"
echo "## e.g. rm -rf $install_dir"
echo "##"
echo "## If you have never used BioSimSpace before, please take a look at the "
echo "## BioSimSpace website at http://biosimspace.org"
echo "##"
echo "##############################################################################"
echo " "

"""
@package biosimspace
@author  Lester Hedges
@brief   Functionality for driving molecular dynamics simulations.
"""

from Sire.Base import findExe

import Sire.System

import BioSimSpace.Process as Process
import BioSimSpace.Protocol as Protocol

from os import environ

import path

# A dictionary mapping MD packages to their executable names and GPU support.
#                PACKAGE        EXE           GPU
_md_packages = { "AMBER"   : { "pmemd.cuda" : True, "pmemd" : False, "sander" : False },
                 "GROMACS" : { "gmx" : True },
                 "NAMD"    : { "namd2" : False }
               }

# A dictionary reverse mapping MD packages to their Sire file extensions.
#                    EXTENSION     EXTENSION
_file_extensions = { "PRM7,RST7"    : "AMBER",
                     "PRM7,RST"     : "AMBER",
                     "GroTop,Gro87" : "GROMACS",
                     "PSF,PDB"      : "NAMD" }

def findMDPackage(system, supports_gpu=False):
    """Find a molecular dynamics package on the system and return
       a handle to it as a MDPackage object.

       Positional arguments:

       system       -- The molecular system.

       Keyword arguments:

       supports_gpu -- Whether GPU support is required.
    """

    # Check that the system is valid.
    if system.__class__ is not Sire.System.System:
        raise TypeError("'system' must be of type 'Sire.System._System.System'")

    if type(supports_gpu) is not bool:
        raise TypeError("'supports_gpu' keyword must be of type 'bool'")

    # Get the file format of the molecular system.
    fileformat = system.property("fileformat").toString()

    # Make sure that this format is supported.
    if not fileformat in _file_extensions:
        raise ValueError("Cannot find an MD package that supports format: %s" % fileformat)
    else:
        package = _file_extensions[fileformat]

    # Get a list of the useable executables for this package. The executables
    # should be listed in order of preference, e.g. performance.
    if supports_gpu:
        # Initalise a list of executables with GPU support.
        gpu_executables = []

        # Loop over all executables associated with the package and append
        # those that have GPU support to the list.
        for exe, gpu in _md_packages[package].items():
            if gpu:
                gpu_executables.append(exe)

        # Raise an error if there is no GPU support.
        if len(gpu_executables) == 0:
            raise ValueError("%s has no GPU support!" % package)
        else:
            executables = gpu_executables
    else:
        executables = _md_packages[package].keys()

    # Loop over each executable in turn and see whether it exists on the system.
    for exe in executables:
        if package == "AMBER":
            # Search AMBERHOME, if set.
            if "AMBERHOME" in environ:
                amber_home = environ.get("AMBERHOME")
                exe = "%s/bin/%s" % (amber_home, exe)
                if path.isfile(exe):
                    return MDPackage(package, exe)

            # Search PATH.
            else:
                try:
                    exe = findExe(exe).absoluteFilePath()
                    return MDPackage(package, exe)
                except:
                    pass
        else:
            try:
                exe = findExe(exe).absoluteFilePath()
                return MDPackage(package, exe)
            except:
                pass

    # If we get this far, then no executable has been found.
    raise ValueError("No executable found for package: %s" % package)

class MDPackage():
    """A class for driving molecular dynamics simulations."""

    def __init__(self, name, exe):
        """Constructor.

           Positional arguments:

           name -- The name of the package.
           exe  -- The executable name.
        """

        # Set the package name and executable.
        # Since MDPackage isn't exposed to the user directly, there is no need
        # for type checking (already validated by findMDPackage).
        self._name = name
        self._exe = exe

    def getName(self):
        """Return the name of the package."""
        return self._name

    def getExe(self):
        """Return the name of the executable."""
        return self._exe

######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for configuring and driving molecular dynamics simulations.
"""

from Sire.Base import findExe as _findExe

from BioSimSpace import _amber_home
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from ..Protocol import Custom as _Custom
from ..Protocol._protocol import Protocol as _Protocol
from .._SireWrappers import System as _System

import BioSimSpace.Process as _Process

import os as _os

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["run"]

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

def _find_md_package(system, protocol, use_gpu=True):
    """Find a molecular dynamics package on the system and return
       a handle to it as a MDPackage object.

       Parameters
       ----------

       system : BioSimSpace._SireWrappers.System
           The molecular system.

       protocol : BioSimSpace.Protocol
           The simulation protocol.

       use_gpu : bool
           Whether to use GPU support.

       Returns
       -------

       (package, exe) : tuple
           The name of the MD package and a path to its executable.
    """

    # Check that the system is valid.
    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    # Check that the use_gpu flag is valid.
    if type(use_gpu) is not bool:
        raise TypeError("'use_gpu' keyword must be of type 'bool'")

    # Get the file format of the molecular system.
    fileformat = system.fileFormat()

    # Make sure that this format is supported.
    if not fileformat in _file_extensions:
        raise ValueError("Cannot find an MD package that supports format: %s" % fileformat)
    else:
        package = _file_extensions[fileformat]

    # Loop over each executable in turn and see whether it exists on the system.
    for exe, gpu in _md_packages[package].items():
        if package == "AMBER":
            # Search AMBERHOME, if set.
            if _amber_home is not None:
                _exe = "%s/bin/%s" % (_amber_home, exe)
                if _os.path.isfile(_exe):
                    return (package, _exe)

        # Search system PATH.
        else:
            try:
                exe = _findExe(exe).absoluteFilePath()
                return (package, exe)
            except:
                pass

    # If we get this far, then no executable was found.
    raise _MissingSoftwareError("No executable found for package: '%s'" % package) from None

def run(system, protocol, autostart=True,
        name="md", work_dir=None, seed=None, property_map={}):
    """Constructor.

        Parameters
        ----------

        system : BioSimSpace._SireWrappers.System
            The molecular system.

        protocol : BioSimSpace.Protocol
            The simulation protocol.

        autostart : bool
            Whether to start the process automatically.

        name : str
            The name of the process.

        work_dir : str
            The working directory for the process.

        seed : int
            A random number seed.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        process : BioSimSpace.Process
            A process to run the molecular dynamics protocol.
    """

    # Check that the system is valid.
    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    # Check that the protocol is valid.
    if not isinstance(protocol, _Protocol):
        if type(protocol) is str:
            protocol = _Custom(protocol)
        else:
            raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol' "
                            "or the path to a custom configuration file.")

    # Find a molecular dynamics package and executable.
    package, exe = _find_md_package(system, protocol)

    # Create the process object.

    # AMBER.
    if package == "AMBER":
        process = _Process.Amber(system, protocol, exe, name, work_dir, seed, property_map)

    # GROMACS.
    elif package == "GROMACS":
        process = _Process.Gromacs(system, protocol, exe, name, work_dir, seed, property_map)

    # NAMD.
    elif package == "NAMD":
        process = _Process.Namd(system, protocol, exe, name, work_dir, seed, property_map)

    # Start the process.
    if autostart:
        return process.start()
    else:
        return process

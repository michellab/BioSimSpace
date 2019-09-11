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

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["run"]

import os as _os

from Sire import Base as _SireBase

from BioSimSpace import _amber_home, _gmx_exe
from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace import Process as _Process
from BioSimSpace import Protocol as _Protocol

# A dictionary mapping MD packages to their executable names and GPU support.
#                PACKAGE        EXE               GPU
_md_packages = { "AMBER"   : { "pmemd.cuda.MPI" : True,
                               "pmemd.cuda"     : True,
                               "pmemd.MPI"      : False,
                               "sander"         : False },
                 "GROMACS" : { "gmx"            : True,
                               "gmx_mpi"        : True },
                 "SOMD"    : { "somd"           : True },
                 "NAMD"    : { "namd2"          : False }
               }

# A dictionary reverse mapping MD packages to their supported Sire file extensions.
# Use SOMD as a fall-back where possible. Since we can't guarantee interconversion
# of potentials for CHARMM-PSF format input files, we restrict such simulations to
# only run using NAMD.
#                    EXTENSION        PACKAGES
_file_extensions = { "PRM7,RST7"    : ["AMBER", "GROMACS", "SOMD"],
                     "PRM7,RST"     : ["AMBER", "GROMACS", "SOMD"],
                     "GroTop,Gro87" : ["GROMACS", "AMBER", "SOMD"],
                     "PSF,PDB"      : ["NAMD"]
                   }

# Whether each package supports free energy simulations. This dictionary needs to
# be updated as support for different packages is added.
_free_energy = { "AMBER"   : False,
                 "GROMACS" : True,
                 "SOMD"    : True,
                 "NAMD"    : False
                }

def _find_md_package(system, protocol, gpu_support=False):
    """Find a molecular dynamics package on the system and return
       a handle to it as a MDPackage object.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The molecular system.

       protocol : :class:`Protocol <BioSimSpace.Protocol>`
           The simulation protocol.

       gpu_support : bool
           Whether to use package must have GPU support.

       Returns
       -------

       (package, exe) : (str, str)
           The name of the MD package and a path to its executable.
    """

    # The input has already been validated in the run method, so no need
    # to re-validate here.

    # Get the file format of the molecular system.
    fileformat = system.fileFormat()

    # Make sure that this format is supported.
    if not fileformat in _file_extensions:
        raise ValueError("Cannot find an MD package that supports format: %s" % fileformat)
    else:
        packages = _file_extensions[fileformat]

    # Is this a free energy protocol.
    if type(protocol) is _Protocol.FreeEnergy:
        is_free_energy = True
    else:
        is_free_energy = False

    # Loop over each package that supports the file format.
    for package in packages:
        # If this is free energy protocol, then check that the package has support.
        if not is_free_energy or _free_energy[package]:
            # Check whether this package exists on the system and has the desired
            # GPU support.
            for exe, gpu in _md_packages[package].items():
                # If the user has requested GPU support make sure the package
                # supports it.
                if not gpu_support or gpu:
                    # AMBER
                    if package == "AMBER":
                        # Search AMBERHOME, if set.
                        if _amber_home is not None:
                            _exe = "%s/bin/%s" % (_amber_home, exe)
                            if _os.path.isfile(_exe):
                                return (package, _exe)
                        # Search system PATH.
                        else:
                            try:
                                exe = _SireBase.findExe(exe).absoluteFilePath()
                                return (package, exe)
                            except:
                                pass
                    # GROMACS
                    elif package == "GROMACS":
                        if _gmx_exe is not None:
                            return (package, _gmx_exe)
                    # SOMD
                    elif package == "SOMD":
                        return (package, _SireBase.getBinDir() + "/somd")
                    # Search system PATH.
                    else:
                        try:
                            exe = _SireBase.findExe(exe).absoluteFilePath()
                            return (package, exe)
                        except:
                            pass

    # If we get this far, then no package was found.
    raise _MissingSoftwareError("Couldn't find package to support format: %s" % fileformat)

def run(system, protocol, gpu_support=False, auto_start=True,
        name="md", work_dir=None, seed=None, property_map={}):
    """Auto-configure and run a molecular dynamics process.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The molecular system.

       protocol : :class:`Protocol <BioSimSpace.Protocol>`
           The simulation protocol.

       gpu_support : bool
           Whether to choose a package with GPU support.

       auto_start : bool
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

       process : :class:`Process <BioSimSpace.Process>`
           A process to run the molecular dynamics protocol.
    """

    # Check that the system is valid.
    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    # Check that the protocol is valid.
    if not isinstance(protocol, _Protocol._protocol.Protocol):
        if type(protocol) is str:
            protocol = _Protocol.Custom(protocol)
        else:
            raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol' "
                            "or the path to a custom configuration file.")

    # Validate optional arguments.

    if type(gpu_support) is not bool:
        raise TypeError("'gpu_support' must be of type 'bool'")

    if type(auto_start) is not bool:
        raise TypeError("'auto_start' must be of type 'bool'")

    if type(name) is not str:
        raise TypeError("'name' must be of type 'str'")

    if work_dir is not None:
        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

    if seed is not None:
        if type(seed) is not int:
            raise TypeError("'seed' must be of type 'int'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Find a molecular dynamics package and executable.
    package, exe = _find_md_package(system, protocol, gpu_support)

    # Create the process object.

    # AMBER.
    if package == "AMBER":
        process = _Process.Amber(system, protocol, exe=exe, name=name,
            work_dir=work_dir, seed=seed, property_map=property_map)

    # GROMACS.
    elif package == "GROMACS":
        process = _Process.Gromacs(system, protocol, exe=exe, name=name,
            work_dir=work_dir, seed=seed, property_map=property_map)

    # SOMD.
    elif package == "SOMD":
        if gpu_support:
            platform = "CUDA"
        else:
            platform = "CPU"
        process = _Process.Somd(system, protocol, exe=exe, name=name,
            work_dir=work_dir, seed=seed, property_map=property_map, platform=platform)

    # NAMD.
    elif package == "NAMD":
        process = _Process.Namd(system, protocol, exe=exe, name=name,
            work_dir=work_dir, seed=seed, property_map=property_map)

    # Start the process.
    if auto_start:
        return process.start()
    else:
        return process

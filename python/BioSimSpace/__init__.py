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
Making biomolecular simulation a breeze!

A collection of tools that makes it easy to write robust and interoperable
molecular workflow components.

www.biosimspace.org
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Align",
           "FreeEnergy",
           "Gateway",
           "IO",
           "Metadynamics",
           "MD",
           "Node",
           "Notebook",
           "Parameters",
           "Process",
           "Protocol",
           "Solvent",
           "Trajectory",
           "Types",
           "Units"]

# Make sure we're using the Sire python interpreter.
try:
    import Sire
    del Sire
except ModuleNotFoundError:
    raise ModuleNotFoundError("BioSimSpace currently requires the Sire "
        + "Python interpreter: www.siremol.org")

# Determine whether we're being imported from a Jupyter notebook.
try:
    shell = get_ipython().__class__.__name__
    if shell == 'ZMQInteractiveShell':
        _is_notebook = True   # Jupyter notebook or qtconsole
    elif shell == 'TerminalInteractiveShell':
        _is_notebook = False  # Terminal running IPython
    else:
        _is_notebook = False  # Other type (?)
except NameError:
    _is_notebook = False      # Probably standard Python interpreter

# Determine whether we're being run interactively.
try:
    shell = get_ipython().__class__.__name__
    if shell == 'ZMQInteractiveShell':
        _is_interactive = True   # Jupyter notebook or qtconsole
    elif shell == 'TerminalInteractiveShell':
        _is_interactive = True   # Terminal running IPython
    else:
        _is_interactive = False  # Other type (?)
except NameError:
    _is_interactive = False      # Probably standard Python interpreter

# Default to non-verbose error messages.
_is_verbose = False

def setVerbose(verbose):
    """Set verbosity of error messages.

       Parameters
       ----------

       verbose : bool
           Whether to print verbose error messages.
    """
    if type(verbose) is not bool:
        raise TypeError("'verbose' must be of type 'bool'.")

    global _is_verbose
    _is_verbose = verbose

def _isVerbose():
    """Whether verbose error messages are active.

       Returns
       ------

       is_verbose : bool
           Whether verbose error messages are active.
    """
    global _is_verbose
    return _is_verbose

from os import environ as _environ
from warnings import warn as _warn

# Check to see if AMBERHOME is set.
if "AMBERHOME" in _environ:
    _amber_home = _environ.get("AMBERHOME")
else:
    _amber_home = None

# Check to see if GROMACS is installed.
from Sire import Base as _SireBase
from os import path as _path

# First, let the user tell us where to find GROMACS. This
# assumes that gromacs is installed in $GROMACSHOME/bin/gmx.
_gmx_exe = None
if "GROMACSHOME" in _environ:
    try:
        _gmx_exe = _SireBase.findExe("%s/bin/gmx" % _environ.get("GROMACSHOME")) \
                            .absoluteFilePath()
    except:
        try:
            _gmx_exe = _SireBase.findExe("%s/bin/gmx_mpi" % _environ.get("GROMACSHOME")) \
                                .absoluteFilePath()
        except:
            pass

if _gmx_exe is None:
    # The user has not told us where it is, so need to look in $PATH.
    try:
        _gmx_exe = _SireBase.findExe("gmx").absoluteFilePath()
    except:
        try:
            _gmx_exe = _SireBase.findExe("gmx_mpi").absoluteFilePath()
        except:
            pass

# Set the bundled GROMACS topology file directory.
_gromacs_path = _path.dirname(_SireBase.getBinDir()) + "/share/gromacs/top"
del _environ
del _SireBase

if not _path.isdir(_gromacs_path):
    _gromacs_path = None

    # Try using the GROMACS exe to get the location of the data directory.
    if _gmx_exe is not None:

        import subprocess as _subprocess

        # Generate the shell command.
        _command = "%s -h 2>&1 | grep 'Data prefix' | awk -F ':' '{print $2}'" % _gmx_exe

        # Run the command.
        _proc = _subprocess.run(_command, shell=True, stdout=_subprocess.PIPE)

        del _command

        # Get the data prefix.
        if _proc.returncode == 0:
            _gromacs_path = _proc.stdout.decode("ascii").strip() + "/share/gromacs/top"
            # Check for the topology file directory.
            if not _path.isdir(_gromacs_path):
                _gromacs_path = None

        del _path
        del _proc
        del _subprocess

from . import Align
from . import FreeEnergy
from . import Gateway
from . import IO
from . import Metadynamics
from . import MD
from . import Node
from . import Notebook
from . import Parameters
from . import Process
from . import Protocol
from . import Solvent
from . import Trajectory
from . import Types
from . import Units

from ._version import get_versions
__version__ = get_versions()['version']
del _version
del get_versions

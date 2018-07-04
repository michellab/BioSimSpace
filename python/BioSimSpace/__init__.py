######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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

# Make sure we're using the Sire python interpreter and store the
# locate of the Sire binary directory.
try:
    import Sire
    del(Sire)
except ModuleNotFoundError:
    raise ModuleNotFoundError("BioSimSpace currently requires the Sire Python interpreter: www.siremol.org")

# Determine whether we're being imported from a Jupyter notebook.
def _is_notebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

# Determine whether we're being run interactively.
def _is_interactive():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return True   # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

from os import environ as _environ
from warnings import warn as _warn

# Check to see if AMBERHOME is set.
if "AMBERHOME" in _environ:
    _amber_home = _environ.get("AMBERHOME")
else:
    _amber_home = None
    _warn("Missing 'AMBERHOME' environment variable. Please download and install AMBER from: http://ambermd.org")
del(_environ)

# Check to see if GROMACS is installed.
import Sire.Base as _SireBase
from os import path as _path
try:
    _gmx_exe = _SireBase.findExe("gmx").absoluteFilePath()
except:
    _gmx_exe = None
    _warn("Cannot find 'gmx' executable. Please download and install GROMACS from: http://gromacs.org")

# Set the bundled GROMACS topology file directory.
_gromacs_path = _path.dirname(_SireBase.getBinDir()) + "/share/gromacs/top"
del(_SireBase)

if not _path.isdir(_gromacs_path):
    # Try using the GROMACS exe to get the location of the data directory.
    if _gmx_exe is not None:

        import subprocess as _subprocess

        # Generate the shell command.
        _command = "%s -h 2>&1 | grep 'Data prefix' | awk -F ':' '{print $2}'" % _gmx_exe

        # Run the command.
        _proc = _subprocess.run(_command, shell=True, stdout=_subprocess.PIPE)

        # Get the data prefix.
        if _proc.returncode == 0:
            _gromacs_path = _proc.stdout.decode("ascii").strip()
            _not_found = False
        else:
            _not_found = True

    if _not_found:
        _warnings.warn("Could not locate GROMACS topology file directory!")

    del(_command)
    del(_not_found)
    del(_proc)
    del(_subprocess)

from BioSimSpace.MD import MD
from BioSimSpace.Trajectory import Trajectory

import BioSimSpace.Gateway as Gateway
import BioSimSpace.IO as IO
import BioSimSpace.Notebook as Notebook
import BioSimSpace.Process as Process
import BioSimSpace.Parameters as Parameters
import BioSimSpace.Protocol as Protocol
import BioSimSpace.Solvent as Solvent
import BioSimSpace.Types as Types

# Top-level functions.

def viewMolecules( files, idxs=None ):
    """View the molecules contained in the passed file(s). Optionally supply
       a list of indices of molecules you want to view. This views the molecules
       and also returns a view object that will allow you to change the view,
       e.g. choosing different molecules to view etc.
    """

    if not _is_notebook():
        _warn("You can only view molecules from within a Jupyter notebook.")
        return None

    if isinstance(files, str):
        files = [files]

    print("Reading molecules from '%s'" % files)
    s = readMolecules(files)

    print("Rendering the molecules...")
    v = BioSimSpace.Notebook.View(s)

    if idxs:
        v.molecules(idxs)
    else:
        v.molecules()

    return v

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
del _warn

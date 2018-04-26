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

from BioSimSpace.MD import MD
from BioSimSpace.Trajectory import Trajectory

import BioSimSpace.Gateway as Gateway
import BioSimSpace.IO as IO
import BioSimSpace.Notebook as Notebook
import BioSimSpace.Process as Process
import BioSimSpace.Protocol as Protocol

from warnings import warn as _warn

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

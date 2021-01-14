######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
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
Functionality for handling parameterisation protocols
for force fields from the Open Force Field Initiative.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["OpenForceField"]

# To override any protocols, just implement a custom "run" method in any
# of the classes.

import os as _os
import queue as _queue
import subprocess as _subprocess
import warnings as _warnings

from Sire import IO as _SireIO
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace import _isVerbose
from BioSimSpace import IO as _IO
from BioSimSpace._Exceptions import ParameterisationError as _ParameterisationError
from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace.Parameters._utils import formalCharge as _formalCharge
from BioSimSpace.Types import Charge as _Charge

from . import _protocol

class OpenForceField(_protocol.Protocol):
    """A class for handling protocols for Open Force Field models."""

    def __init__(self, forcefield, property_map={}):
        """Constructor.

           Parameters
           ----------

           forcefield : str
               The name of the force field.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield=forcefield, property_map=property_map)

        # Set the compatibility flags.
        self._tleap = False
        self._pdb2gmx = False

    def run(self, molecule, work_dir=None, queue=None):
        """Run the parameterisation protocol.

           Parameters
           ----------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule to apply the parameterisation protocol to.

           work_dir : str
               The working directory.

           queue : queue.Queue
               The thread queue is which this method has been run.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The parameterised molecule.
        """

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

        if type(work_dir) is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        if type(queue) is not None and type(queue) is not _queue.Queue:
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Create the file prefix.
        prefix = work_dir + "/"

        # Create a copy of the molecule.
        new_mol = molecule.copy()

        if queue is not None:
            queue.put(new_mol)
        return new_mol

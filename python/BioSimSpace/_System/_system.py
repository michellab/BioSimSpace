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
A thin wrapper around Sire.System. This is an internal package and should
not be directly exposed to the user.

Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Mol as _SireMol
import Sire.System as _SireSystem

__all__ = ["System"]

class _MolWithResName(_SireMol.MolWithResID):
    def __init__(self, resname):
        super().__init__(_SireMol.ResName(resname))

class System():
    """A container class for storing molecular systems."""

    def __init__(self, system):
        """Constructor.

           Positional arguments:

           system -- A Sire System object.
        """

        # Check that the system is valid.
        if not isinstance(system, _SireSystem.System):
            raise TypeError("'system' must be of type 'Sire.System._System.System'")

        # Set the system.
        self._system = system

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def nMolecules(self):
        """Return the number of molecules in the system."""
        return self._system.nMolecules()

    def fileFormat(self):
        """Return the file formats associated with the system."""
        return self._system.property("fileformat").value()

    def getMolWithResName(self, resname):
        """Return the molecule containing the given residue.

           Positional arguments:

           resname -- The name of a residue unique to the molecule.
        """
        try:
            return self._system[_MolWithResName(resname)]
        except:
            raise KeyError("System does not contain residue '%s'" % resname)

    def _getSireSystem(self):
        """Return the full Sire System object."""
        return self._system

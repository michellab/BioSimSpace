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
A thin wrapper around Sire.Mol.Residue. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Residue"]

import os.path as _path
import random as _random
import string as _string

import Sire.Maths as _SireMaths
import Sire.Mol as _SireMol

from ._sire_wrapper import SireWrapper as _SireWrapper

class Residue(_SireWrapper):
    """A class for storing a residue."""

    def __init__(self, residue):
        """Constructor.

           Parameters
           ----------

           atom : Sire.Mol.Residue, Atom :class:`Atom <BioSimSpace._SireWrappers.Atom>`
               A Sire or BioSimSpace Residue object.
        """

        # Check that the residue is valid.

        # A Sire Residue object.
        if type(residue) is _SireMol.Residue:
            sire_object = residue

        # Another BioSimSpace Residue object.
        elif type(residue) is Residue:
            sire_object = residue._sire_object

        # Invalid type.
        else:
            raise TypeError("'residue' must be of type 'Sire.Mol.Residue' "
                            "or 'BioSimSpace._SireWrappers.Residue'.")

        # Call the base class constructor.
        super().__init__(sire_object)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Residue: name=%r, molecule=%d, index=%d, nAtoms=%d>" \
            % (self.name(), self.moleculeNumber(), self.index(), self.nAtoms())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Residue: name=%r, molecule=%d, index=%d, nAtoms=%d>" \
            % (self.name(), self.moleculeNumber(), self.index(), self.nAtoms())

    def name(self):
        """Return the name of the residue.

           Returns
           -------

           name : str
               The name of the residue.
        """
        return self._sire_object.name().value()

    def index(self):
        """Return the index of the residue.

           Returns
           -------

           index : int
               The index of the residue.
        """
        return self._sire_object.index().value()

    def moleculeNumber(self):
        """Return the number of the molecule to which this residue belongs.

           Returns
           -------

           number : int
               The number of the molecule to which the residue belongs.
        """
        return self._sire_object.molecule().number().value()

    def nAtoms(self):
        """Return the number of atoms in the residue.

           Returns
           -------

           num_atoms : int
               The number of atoms in the system.
        """
        return self._sire_object.nAtoms()

    def getAtoms(self):
        """Return a list containing all of the atoms in the residue.

           Parameters
           ----------

           Returns
           -------

           atoms : [:class:`Atoms <BioSimSpace._SireWrappers.Atom>`]
               The list of atoms in the residue.
        """
        atoms = []
        for atom in self._sire_object.atoms():
            atoms.append(_Atom(atom))
        return atoms

    def toMolecule(self):
        """Convert a single Residue to a Molecule.

           Returns
           -------

           system : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
        """
        return _Molecule(_SireMol.PartialMolecule(self._sire_object).extract().molecule())

    def search(self, query):
        """Search the residue for atoms and residues.

           Parameters
           ----------

           query : str
               The search query.

           Returns
           -------

           results : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`]
               A list of objects matching the search query.

           Examples
           --------

           Search for all oxygen or hydrogen atoms.

           >>> result = residue.search("element oxygen or element hydrogen")

           Search for atom index 23.

           >>> result = residue.search("atomidx 23")
        """

        if type(query) is not str:
            raise TypeError("'query' must be of type 'str'")

        # Intialise a list to hold the search results.
        results = []

        try:
            # Query the Sire system.
            search_result = self._sire_object.search(query)

        except:
            raise ValueError("'Invalid search query: %r" % query) from None

        return _SearchResult(search_result)

# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecule import Molecule as _Molecule
from ._search_result import SearchResult as _SearchResult

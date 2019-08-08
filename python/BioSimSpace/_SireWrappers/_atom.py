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
A thin wrapper around Sire.Mol.Atom. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Atom"]

import os.path as _path
import random as _random
import string as _string

from Sire import Mol as _SireMol

from ._sire_wrapper import SireWrapper as _SireWrapper

class Atom(_SireWrapper):
    """A class for storing an atom."""

    def __init__(self, atom):
        """Constructor.

           Parameters
           ----------

           atom : Sire.Mol.Atom, :class:`Atom <BioSimSpace._SireWrappers.Atom>`
               A Sire or BioSimSpace Atom object.
        """

        # Check that the atom is valid.

        # A Sire Atom object.
        if type(atom) is _SireMol.Atom:
            sire_object = atom

        # Another BioSimSpace Atom object.
        elif type(atom) is Atom:
            sire_object = atom._sire_object

        # Invalid type.
        else:
            raise TypeError("'atom' must be of type 'Sire.Mol.Atom' "
                            "or 'BioSimSpace._SireWrappers.Atom'.")

        # Call the base class constructor.
        super().__init__(sire_object)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Atom: name=%r, molecule=%d index=%d>" \
            % (self.name(), self.moleculeNumber(), self.index())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Atom: name=%r, molecule=%d index=%d>" \
            % (self.name(), self.moleculeNumber(), self.index())

    def name(self):
        """Return the name of the atom.

           Returns
           -------

           name : str
               The name of the atom.
        """
        return self._sire_object.name().value()

    def index(self):
        """Return the index of the atom.

           Returns
           -------

           index : int
               The index of the atom.
        """
        return self._sire_object.index().value()

    def moleculeNumber(self):
        """Return the number of the molecule to which this atom belongs.

           Returns
           -------

           number : int
               The number of the molecule to which the atom belongs.
        """
        return self._sire_object.molecule().number().value()

    def element(self, property_map={}):
        """Return the element.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           element : str
               The element.
        """

        prop = property_map.get("element", "element")

        # Calculate the element.
        try:
            element = self._sire_object.property(prop).toString()
        except:
            element = ""

        # Return the element.
        return element

    def toMolecule(self):
        """Convert a single Atom to a Molecule.

           Returns
           -------

           system : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
        """
        return _Molecule(_SireMol.PartialMolecule(self._sire_object).extract().molecule())

# Import at bottom of module to avoid circular dependency.
from ._molecule import Molecule as _Molecule

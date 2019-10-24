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
Base class for wrapped Sire objects. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["SireWrapper"]

from Sire import Maths as _SireMaths
from Sire import Mol as _SireMol
from Sire import Vol as _SireVol

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace import Units as _Units

class SireWrapper():
    """A base class for wrapping Sire objects."""

    def __init__(self, object):
        """Constructor.

           Parameters
           ----------

           object : Sire.System.System, Sire.Mol.Molecule, Sire.Mol.Residue, Sire.Mol.Atom
               A Sire object.
        """

        # Store a deep copy of the Sire object.
        self._sire_object = object.__deepcopy__()

        # Intialise flags.
        self._is_multi_atom = False
        self._is_merged = False

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object == other._sire_object
        else:
            return False

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object != other._sire_object
        else:
            return False

    def __hash__(self):
        """Hash operator."""
        return hash(self._sire_object)

    def copy(self):
        """Return a copy of this object. The return type is same as the object
           on which copy is called.

           Returns
           -------

           system : :class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                    :class:`Residue <BioSimSpace._SireWrappers.Residue>`, \
                    :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                    :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                    :class:`System <BioSimSpace._SireWrappers.System>`
               A copy of the object.
        """
        return type(self)(self)

    def charge(self, property_map={}, is_lambda1=False):
        """Return the charge.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           is_lambda1 : bool
              Whether to use the charge at lambda = 1 if the object is peturbable.

           Returns
           -------

           charge : :class:`Charge <BioSimSpace.Types.Charge>`
               The charge.
        """

        # Copy the map.
        _property_map = property_map.copy()

        # This is a merged molecule.
        if self._is_merged:
            if is_lambda1:
                _property_map = { "charge" : "charge1" }
            else:
                _property_map = { "charge" : "charge0" }

        # Calculate the charge.
        try:
            charge = self._sire_object.evaluate().charge(_property_map).value()
        except:
            charge = 0

        # Return the charge.
        return charge * _Units.Charge.electron_charge

    def translate(self, vector, property_map={}):
        """Translate the object.

           Parameters
           ----------

           vector : [:class:`Length <BioSimSpace.Types.Length>`]
               The translation vector in Angstroms.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Convert tuple to a list.
        if type(vector) is tuple:
            vector = list(vector)

        # Validate input.
        if type(vector) is list:
            vec = []
            for x in vector:
                if type(x) is int:
                    vec.append(float(x))
                elif type(x) is float:
                    vec.append(x)
                elif type(x) is _Length:
                    vec.append(x.angstroms().magnitude())
                else:
                    raise TypeError("'vector' must contain 'int', 'float', or "
                                    "'BioSimSpace.Types.Length' types only!")
        else:
            raise TypeError("'vector' must be of type 'list' or 'tuple'")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        try:
            # Make a local copy of the property map.
            _property_map = property_map.copy()

            if "coordinates" not in property_map and self._is_merged:
                _property_map["coordinates"] = "coordinates0"

            # Perform the translation.
            self._sire_object = self._sire_object                                     \
                                    .move()                                           \
                                    .translate(_SireMaths.Vector(vec), _property_map) \
                                    .commit()
        except UserWarning as e:
            msg = "Cannot compute axis-aligned bounding box " + \
                   "since the object has no 'coordinates' property."
            if _isVerbose():
                raise _IncompatibleError(msg) from e
            else:
                raise _IncompatibleError(msg) from None

    def getAxisAlignedBoundingBox(self, property_map={}):
        """Get the axis-aligned bounding box enclosing the object.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           box_min : [:class:`Length <BioSimSpace.Types.Length>`]
               The minimum coordinates of the axis-aligned bounding box in
               each dimension.

           box_max : [:class:`Length <BioSimSpace.Types.Length>`]
               The minimum coordinates of the axis-aligned bounding box in
               each dimension.
        """
        aabox = self._getAABox(property_map)

        box_min = [x*_Units.Length.angstrom for x in aabox.minCoords()]
        box_max = [x*_Units.Length.angstrom for x in aabox.maxCoords()]

        return box_min, box_max

    def _getSireObject(self):
        """Return the underlying Sire object.

           Returns
           -------

           object : Sire.System.System, Sire.Mol.Molecule, Sire.Mol.Residue, Sire.Mol.Atom
        """
        return self._sire_object

    def _getAABox(self, property_map={}):
        """Get the axis-aligned bounding box for the object.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           aabox : Sire.Vol.AABox
               The axis-aligned bounding box for the object.
        """

        # Initialise the coordinates vector.
        coord = []

        prop = property_map.get("coordinates", "coordinates")

        # Handle merged molecules.
        if self._is_merged:
            prop = "coordinates0"

        # Residues need to be converted to molecules to have a
        # coordinates property.
        try:
            c = self._sire_object.property(prop)
        except:
            try:
                c = self.toMolecule()._sire_object.property(prop)
            except Exception as e:
                msg = "Cannot compute axis-aligned bounding box " + \
                      "since the object has no 'coordinates' property."
                if _isVerbose():
                    print(msg)
                    raise _IncompatibleError(msg) from e
                else:
                    raise _IncompatibleError(msg) from None

        # We have a vector of coordinates. (Multiple atoms)
        if self._is_multi_atom:
            coord.extend(c.toVector())
        # Convert to a list.
        else:
            coord = [c]

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)

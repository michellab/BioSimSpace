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
An energy type.
"""

import Sire.Units as _Units

from ._type import Type as _Type

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Energy"]

class Energy(_Type):
    """An energy type."""

    # Dictionary of allowed units.
    _supported_units = { "KILO CALORIES PER MOL" : _Units.kcal_per_mol,
                         "KILO JOULES PER MOL"   : _Units.kJ_per_mol,
                         "KT"                    : 2.479 * _Units.kJ_per_mol }

    # Map unit abbreviations to the full name.
    _abbreviations = { "KCAL/MOL" : "KILO CALORIES PER MOL",
                       "KJ/MOL"   : "KILO JOULES PER MOL",
                       "KT"       : "KT" }

    # Print formatting.
    _print_format = { "KILO CALORIES PER MOL" : "kcal/mol",
                      "KILO JOULES PER MOL"   : "kJ/mol",
                      "KT"                    : "KT" }

    # Documentation strings.
    _doc_strings = { "KILO CALORIES PER MOL" : "An energy in kcal per mol.",
                     "KILO JOULES PER MOL"   : "An energy in kJ per mol.",
                     "KT"                    : "An energy in KT." }

    def __init__(self, *args):
        """Constructor.

           Parameters
           ----------

           magnitude : float
               The magnitude.

           unit : str
               The unit.

           string : str
               A string representation of the energy.
        """

        # Call the base class constructor.
        super().__init__(*args)

    def kcal_per_mol(self):
        """Return the energy in kcal per mol.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The energy in kcal per mol.
        """
        return Energy((self._magnitude * self._supported_units[self._unit]).to(_Units.kcal_per_mol), "KILO CALORIES PER MOL")

    def kj_per_mol(self):
        """Return the energy in kJ per mol.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The energy in kJ per mol.
        """
        return Energy((self._magnitude * self._supported_units[self._unit]).to(_Units.kJ_per_mol), "KILO JOULES PER MOL")

    def kt(self):
        """Return the energy in KT.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The energy in KT.
        """
        return Energy((self._magnitude * self._supported_units[self._unit]).to(2.479 * _Units.kJ_per_mol), "KT")

    def _default_unit(self, mag=None):
        """Internal method to return an object of the same type in the default unit.

           Parameters
           ----------

           mag : float
               The magnitude (optional).

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The energy in the default unit of kcal per mol.
        """
        if mag is None:
            return self.kcal_per_mol()
        else:
            return Energy(mag, "KILO CALORIES PER MOL")

    def _convert_to(self, unit):
        """Return the energy in a different unit.

           Parameters
           ----------

           unit : str
               The unit to convert to.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The energy in the specified unit.
        """
        if unit == "KILO CALORIES PER MOL":
            return self.kcal_per_mol()
        elif unit == "KILO JOULES PER MOL":
            return self.kj_per_mol()
        elif unit == "KT":
            return self.kt()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace all instances of KILO with K.
        unit = unit.replace("KILO", "K")

        # Replace all instances of PER with the / character.
        unit = unit.replace("PER", "/")

        # Replace all instance of MOLE with MOL.
        unit = unit.replace("MOLE", "MOL")

        # Replace all instances of CALORIES with CAL.
        unit = unit.replace("CALORIES", "CAL")

        # Replace all instances of JOULES with J.
        unit = unit.replace("JOULES", "J")

        # Check that the unit is supported.
        if unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

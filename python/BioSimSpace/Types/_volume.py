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
A volume type.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Volume"]

from Sire import Units as _SireUnits

from ._type import Type as _Type

class Volume(_Type):
    """A volume type."""

    # Dictionary of allowed units.
    _supported_units = { "METER3"      : _SireUnits.meter3,
                         "NANOMETER3"  : _SireUnits.nanometer3,
                         "ANGSTROM3"   : _SireUnits.angstrom3,
                         "PICOMETER3"  : _SireUnits.picometer3 }

    # Map unit abbreviations to the full name.
    _abbreviations = { "M3"  : "METER3",
                       "NM3" : "NANOMETER3",
                       "A3"  : "ANGSTROM3",
                       "PM3" : "PICOMETER3" }

    # Dictionary of allowed units.
    _print_format = { "METER3"      : "m^3",
                      "NANOMETER3"  : "nm^3",
                      "ANGSTROM3"   : "A^3",
                      "PICOMETER3"  : "pm^3" }

    # Documentation strings.
    _doc_strings = { "METER3"      : "A volume in cube meters.",
                     "NANOMETER3"  : "A volume in cube nanometers.",
                     "ANGSTROM3"   : "A volume in cube Angstrom.",
                     "PICOMETER3"  : "A volume in cube picometers." }

    # Null type unit for avoiding issue printing configargparse help.
    _null_unit = "ANGSTROM3"

    def __init__(self, *args):
        """Constructor.

           ``*args`` can be a magnitude and unit, or a string representation
           of the volume, e.g. "100 pm^3".

           Parameters
           ----------

           magnitude : float
               The magnitude.

           unit : str
               The unit.

           string : str
               A string representation of the volume.

           Examples
           --------

           Create an object representing a volume of 100 cube nanometers then
           print the volume in cube Angstrom.

           >>> import BioSimSpace as BSS
           >>> volume = BSS.Types.Volume(100, "nm^3")
           >>> print(volume.angstroms3())

           The same as above, except passing a string representation of the
           volume to the constructor.

           >>> import BioSimSpace as BSS
           >>> volume = BSS.Types.Volume("100 nm^3")
           >>> print(volume.angstroms3())

           The string matching is extremeley flexible, so all of the following
           would be valid arguments: "100 nm^3", "100 cube nanometers",
           "100 nanometers cubed".
        """

        # Call the base class constructor.
        super().__init__(*args)

        # Don't support negative volumes.
        if self._magnitude < 0:
            raise ValueError("The volume cannot be negative!")

    def __truediv__(self, other):
        """Division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if type(other) is float:
            mag = self._magnitude / other
            return Volume(mag, self._unit)

        # Division by another Volume.
        elif type(other) is Volume:
            return self.angstroms3().magnitude() / other.angstroms3().magnitude()

        # Division by an Area.
        elif type(other) is _Area:
            mag = self.angstroms3().magnitude() / other.angstroms2().magnitude()
            return _Length(mag, "A")

        # Division by a Length.
        elif type(other) is _Length:
            mag = self.angstroms3().magnitude() / other.angstroms().magnitude()
            return _Area(mag, "A2")

        # Division by a string.
        elif type(other) is str:
            try:
                length = _Length(other)
                return self / length
            except:
                try:
                    area = _Area(other)
                    return self / area
                except:
                    try:
                        volume = Volume(other)
                        return self / volume
                    except:
                        raise ValueError("Could not convert the string to a "
                                         "'BioSimSpace.Types.Length', 'BioSimSpace.Types.Area', "
                                         "or 'BioSimSpace.Types.Volume'")
        else:
            raise TypeError("unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def meters3(self):
        """Return the volume in cubic meters.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
               The volume in cubic meters.
        """
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_SireUnits.meter3), "METER3")

    def nanometers3(self):
        """Return the volume in cubic nanometers.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
               The volume in cubic nanometers.
        """
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_SireUnits.nanometer3), "NANOMETER3")

    def angstroms3(self):
        """Return the volume in cubic Angstrom.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
               The volume in cubic Angstrom.
        """
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_SireUnits.angstrom3), "ANGSTROM3")

    def picometers3(self):
        """Return the volume in cubic picometers.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
               The volume in cubic picometers.
        """
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_SireUnits.picometer3), "PICOMETER3")

    def _default_unit(self, mag=None):
        """Internal method to return an object of the same type in the default unit.

           Parameters
           ----------

           mag : float
               The magnitude (optional).

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
               The volume in the default unit of cubic Angstrom.
        """
        if mag is None:
            return self.angstroms3()
        else:
            return Volume(mag, "ANGSTROM3")

    def _convert_to(self, unit):
        """Return the volume in a different unit.

           Parameters
           ----------

           unit : str
               The unit to convert to.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
               The volume in the specified unit.
        """
        if unit == "METER3":
            return self.meters3()
        elif unit == "NANOMETER3":
            return self.nanometers3()
        elif unit == "ANGSTROM3":
            return self.angstroms3()
        elif unit == "PICOMETER3":
            return self.picometers3()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace any occurence of cubed with 3.
        unit = unit.replace("CUBED", "3")
        unit = unit.replace("CUBE", "3")

        # Strip "^" character.
        unit = unit.replace("^", "")

        # Strip any "S" characters.
        unit = unit.replace("S", "")

        # Fix for ANGSTROM (since it contains an "S").
        if unit[0:3] == "ANG":
            unit = "ANGS" + unit[3:]

        # Make sure that the "3" character appears last. This allows the user
        # to write, e.g. "cube nm" or "nm cubed".
        index = unit.find("3")
        if index != -1:
            unit = unit[0:index] + unit[index+1:] + "3"

        # Check that the unit is supported.
        if unit in self._supported_units:
            return unit
        elif unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

# Import at bottom of module to avoid circular dependency.
from ._area import Area as _Area
from ._length import Length as _Length

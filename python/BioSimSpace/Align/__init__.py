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
Functionality for aligning molecules.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from .._SireWrappers import Molecule as _Molecule

def matchAtoms(molecule0, molecule1, scoring_function="", matches=1):
    """Find mappings from atoms in molecule0 to those in molecule1.

       Positional arguments:

       molecule0        -- The reference molecule.
       molecule1        -- The target molecule.

       Keyword arguments:

       scoring_function -- The scoring function used to match atoms. Available
                           options are: "something", "something else", ...
       matches          -- The number of matches to return. (Sorted in order of score).
    """

    # A list of supported scoring functions.
    scoring_functions = ["something", "something else..."]

    if type(molecule0) is not _Molecule:
        raise TypeError("'molecule0' must be of type 'BioSimSpace.Molecule'")

    if type(molecule1) is not _Molecule:
        raise TypeError("'molecule1' must be of type 'BioSimSpace.Molecule'")

    if type(scoring_function) is not str:
        raise TypeError("'scoring_function' must be of type 'str'")
    else:
        if not scoring_function in scoring_functions:
            raise ValueError("Unsupported scoring function '%s'. Options are: %s"
                % (scoring_function, scoring_functions))

    if type(matches) is not int:
        raise TypeError("'matches' must be of type 'int'")
    else:
        if matches < 0:
            raise ValueError("'matches' must be positive!")

    return None

def rmsdAlign(molecule0, molecule1, mapping):
    """Align atoms in molecule0 to those in molecule1.

       Positional arguments:

       molecule0 -- The reference molecule.
       molecule1 -- The target molecule.
       mapping   -- A dictionary mapping atoms in molecule0 to those in molecule1.
    """

    if type(molecule0) is not _Molecule:
        raise TypeError("'molecule0' must be of type 'BioSimSpace.Molecule'")

    if type(molecule1) is not _Molecule:
        raise TypeError("'molecule1' must be of type 'BioSimSpace.Molecule'")

    if type(mapping) is not dict:
        raise TypeError("'mapping' must be of type 'dict'.")

    return None

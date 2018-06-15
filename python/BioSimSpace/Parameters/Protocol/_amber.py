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
Functionality for handling parameterisation protocols
for AMBER force field models.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

# To override any protocols, just implement a custom "run" method in any
# of the classes. Note that all input/out files must be prefixed with
# "input" and "output" respectively. See the base class "run" methods
# in _protocol.py for examples.

from . import _protocol
from ..._SireWrappers import Molecule as _Molecule

__all__ = ["FF03", "FF99", "FF99SB", "FF14SB", "GAFF", "GAFF2"]

class FF03(_protocol.Protocol):
    """A class for handling protocols for the FF03 force field model."""

    def __init__(self):
        """Constructor."""

        # Call the base class constructor.
        super().__init__(forcefield="ff03")

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF99(_protocol.Protocol):
    """A class for handling protocols for the FF99 force field model."""

    def __init__(self):
        """Constructor."""

        # Call the base class constructor.
        super().__init__(forcefield="ff99")

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF99SB(_protocol.Protocol):
    """A class for handling protocols for the FF99SB force field model."""

    def __init__(self):
        """Constructor."""

        # Call the base class constructor.
        super().__init__(forcefield="ff99SB")

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF14SB(_protocol.Protocol):
    """A class for handling protocols for the FF14SB force field model."""

    def __init__(self):
        """Constructor."""

        # Call the base class constructor.
        super().__init__(forcefield="ff14SB")

        # Set the compatibility flags.
        self._tleap = True

class GAFF(_protocol.Protocol):
    """A class for handling protocols for the GAFF force field model."""

    def __init__(self):
        """Constructor."""

        # Call the base class constructor.
        super().__init__(forcefield="gaff")

        self._version = 1

    def run(self):
        """Run the protocotol.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
        """

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

class GAFF2(_protocol.Protocol):
    """A class for handling protocols for the GAFF2 force field model."""

    # Copy the GAFF run method.
    run = GAFF.run

    def __init__(self):
        """Constructor."""

        # Call the base class constructor.
        super().__init__(forcefield="gaff2")

        self._version = 2

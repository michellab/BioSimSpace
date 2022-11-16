######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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
Functionality for minimisation protocols.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Minimisation"]

import warnings as _warnings

from ._protocol import Protocol as _Protocol
from ._position_restraint import _PositionRestrain


class Minimisation(_Protocol, _PositionRestrain):
    """A class for storing minimisation protocols."""

    def __init__(
        self,
        steps=10000,
        restraint=None,
        force_constant=None,
    ):
        """Constructor.

        Parameters
        ----------

        steps : int
            The maximum number of steps to perform.

        restraint : str, [int]
            The type of restraint to perform. This should be one of the
            following options:
                "backbone"
                     Protein backbone atoms. The matching is done by a name
                     template, so is unreliable on conversion between
                     molecular file formats.
                "heavy"
                     All non-hydrogen atoms that aren't part of water
                     molecules or free ions.
                "all"
                     All atoms that aren't part of water molecules or free
                     ions.
            Alternatively, the user can pass a list of atom indices for
            more fine-grained control. If None, then no restraints are used.

        force_constant : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
            The force constant for the restraint potential. If a 'float' is
            passed, then default units of 'kcal_per_mol / angstrom**2' will
            be used.
        """

        # Call the base class constructor.
        _Protocol.__init__(self)
        _PositionRestrain.__init__(self, restraint, force_constant)

        # Set the number of steps.
        self.setSteps(steps)

    def _get_parm(self):
        """Return a string representation of the parameters."""
        return f"steps={self._steps}, " + _PositionRestrain._get_parm(self)

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol.Minimisation: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "BioSimSpace.Protocol.Custom"
        else:
            return f"BioSimSpace.Protocol.Minimisation({self._get_parm()})"

    def getSteps(self):
        """Return the maximum number of steps.

        Returns
        -------

        steps : int
            The maximum number of minimisation steps.
        """
        return self._steps

    def setSteps(self, steps):
        """Set the maximum number of steps.

        Parameters
        ----------

        steps : int
            The maximum number of minimisation steps.
        """
        if not type(steps) is int:
            raise TypeError("'steps' must be of type 'int'")

        if steps <= 0:
            _warnings.warn(
                "Number of steps must be greater than zero. Using default (10000)."
            )
            self._steps = 10000

        else:
            self._steps = steps

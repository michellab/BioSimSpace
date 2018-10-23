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
Functionality for running solvation free energy calculations.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from . import _free_energy
from .._SireWrappers import System as _System

class Solvation(_free_energy.FreeEnergy):
    """A class for configuring and running solvation free energy simulations."""

    def __init__(self, system, protocol=None, work_dir=None):
        """Constructor.

           Positional arguments
           --------------------

           system : BioSimSpace._SireWrappers.System
               The molecular system.


           Keyword arguments
           -----------------

           protocol : BioSimSpace.Protocol.FreeEnergy
               The simulation protocol.

           work_dir : str
               The working directory for the simulation.
        """

        # Call the base class constructor.
        super().__init__(protocol, work_dir)

        # Validate the input.

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of the solvated system.
            self._system_solvated = _System(system._getSireSystem().__deepcopy__())

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError("The system must contain a single perturbable molecule!")

            # The system must be solvated.
            if system.nWaterMolecules() == 0:
                raise ValueError("The system must be solvated!")

            # Create the vacuum system.
            self._system_vacuum = _System(system._getSireSystem().__deepcopy__())
            self._system_vacuum.removeWaterMolecules()

        # Initialise the process runner with all of the simulations required
        # for each leg.
        self._initialise_runner(self._system_solvated, self._system_vacuum)

    def analyse(self):
        """Analyse the solvation free energy data.

           Returns
           -------

           pmf_free : [ ( float, BioSimSpace.Types.Energy, BioSimSpace.Types.Energy ) ]
               The potential of mean force (PMF) for the free leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf_vacuum : [ ( float, BioSimSpace.Types.Energy, BioSimSpace.Types.Energy ) ]
               The potential of mean force (PMF) for the vacuum leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : ( BioSimSpace.Types.Energy, BioSimSpace.Types.Energy )
               The solvation free energy difference and its associated error.
        """

        # This method is just a wrapper to provide simulation specific doc
        # strings. We just call the base class method, which is aware of
        # the simulation type.
        return super().analyse()

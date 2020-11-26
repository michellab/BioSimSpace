######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2020
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

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Solvation"]

from BioSimSpace._SireWrappers import System as _System

from . import _free_energy

class Solvation(_free_energy.FreeEnergy):
    """A class for configuring and running solvation free energy simulations."""

    def __init__(self, system, protocol=None, vacuum_leg=True,
            work_dir=None, engine=None, multistep=False, steps=None):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`
               The simulation protocol.

           vacuum_leg : bool
               Whether to simulation the vacuum leg of the simulation. Set to False
               if you only wish to run the free leg.

           work_dir : str
               The working directory for the simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "GROMACS", or "SOMD". If this argument is omitted then
               BioSimSpace will choose an appropriate engine for you.

           multistep : bool
               Whether or not to perform a multistep perturbation which splits perturbing
               the terms (LJ/charge/bond) instead of all at the same time.  

           steps: dict
               Optional dictionary of specified terms to perturb as keys and number of lambda windows
               per term as values. The full dictionary should look like:
               {
               term1 : n_1,
               term2 : n_2,
               term3 : n_3,
               term4 : n_4,
               term5 : n_5,
               term6 : n_6                             
               }
               if not defined while multistep==True, will default to setting up all terms with 5 lambda
               windows each.

        """

        # Call the base class constructor.
        super().__init__(protocol, work_dir, engine)

        # Validate the input.

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of the solvated system. (Used for the first leg.)
            self._system0 = system.copy()

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError("The system must contain a single perturbable molecule!")

            # The system must be solvated.
            if system.nWaterMolecules() == 0:
                raise ValueError("The system must be solvated!")

            # Create the vacuum system. (Used for the second leg.)
            self._system1 = system.copy()
            self._system1.removeWaterMolecules()

        if type(vacuum_leg) is not bool:
            raise TypeError("'vacuum_leg' must be of type 'bool.")
        else:
            if not vacuum_leg:
                self._is_dual = False

        if type(multistep) is not bool:
          raise TypeError("'multistep' must be of type 'bool'.")

        if steps:
          if type(steps) is not dict:
            raise TypeError("'steps' must be of type 'dict'.")

        # Initialise the process runner with all of the simulations required
        # for each leg.
        if not multistep:
          self._initialise_runner(self._system0, self._system1)

        # multistep approach.
        else:
          if not steps:

            # do stuff here. How should this call change?
            self._initialise_runner(self._system0, self._system1)
            # assume all terms should be perturbed with 5 lambda windows each.





          else:
            # perturb user-specified terms only.
            raise NotImplementedError("Custom multistep approach not yet implemented.")



    def analyse(self):
        """Analyse the solvation free energy data.

           Returns
           -------

           pmf_free : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the free leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf_vacuum : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the vacuum leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The solvation free energy difference and its associated error.
        """

        # This method is just a wrapper to provide simulation specific doc
        # strings. We just call the base class method, which is aware of
        # the simulation type.
        if self._engine == "SOMD":
            return super()._analyse_somd()
        elif self._engine == "GROMACS":
            return super()._analyse_gromacs()

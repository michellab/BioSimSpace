######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
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
__email__ = "lester.hedges@gmail.com"

__all__ = ["Solvation"]

from BioSimSpace._SireWrappers import System as _System

from . import _free_energy

class Solvation(_free_energy.FreeEnergy):
    """A class for configuring and running solvation free energy simulations."""

    def __init__(self, system, protocol=None, vacuum_leg=True,
            work_dir=None, engine=None, setup_only=False,
            ignore_warnings=False, show_errors=True):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`, \
                     [:class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`,
                      :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`]
               The simulation protocol. If one obect is passed, then the same
               protocol will be used for both legs of the simulation. Passing
               two objects enables a different protocol for each leg, e.g. a
               different lambda schedule, or run time.

           vacuum_leg : bool
               Whether to simulation the vacuum leg of the simulation. Set to False
               if you only wish to run the free leg.

           work_dir : str
               The working directory for the simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "GROMACS", or "SOMD". If this argument is omitted then
               BioSimSpace will choose an appropriate engine for you.

           setup_only: bool
               Whether to only support simulation setup. If True, then no
               simulation processes objects will be created, only the directory
               hierarchy and input files to run a simulation externally. This
               can be useful when you don't intend to use BioSimSpace to run
               the simulation. Note that a 'work_dir' must also be specified.

           ignore_warnings : bool
               Whether to ignore warnings when generating the binary run file.
               This option is specific to GROMACS and will be ignored when a
               different molecular dynamics engine is chosen.

           show_errors : bool
               Whether to show warning/error messages when generating the binary
               run file. This option is specific to GROMACS and will be ignored
               when a different molecular dynamics engine is chosen.
        """

        # Call the base class constructor.
        super().__init__(protocol, work_dir, engine, setup_only=setup_only,
            ignore_warnings=ignore_warnings, show_errors=show_errors)

        # Validate the input.

        if type(vacuum_leg) is not bool:
            raise TypeError("'vacuum_leg' must be of type 'bool.")
        else:
            if not vacuum_leg:
                self._is_dual = False

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of the solvated system. (Used for the first leg.)
            self._system0 = system.copy()

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError("The system must contain a single perturbable molecule! "
                                 "Use the 'BioSimSpace.Align' package to map and merge molecules.")

            # The system must be solvated.
            if system.nWaterMolecules() == 0:
                raise ValueError("The system must be solvated! Use the 'BioSimSpace.Solvent' "
                                 "package to solvate your system.")

            # Create the vacuum system. (Used for the second leg.)
            if vacuum_leg:
                self._system1 = system.copy()
                self._system1.removeWaterMolecules()
            else:
                self._system1 = None

        # Initialise the process runner with all of the simulations required
        # for each leg.
        self._initialise_runner(self._system0, self._system1)

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

           overlap_free : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the free leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.

           overlap_vacuum : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the vacuum leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.
        """

        # This method is just a wrapper to provide simulation specific doc
        # strings. We just call the base class method, which is aware of
        # the simulation type.
        if self._engine == "SOMD":
            return super()._analyse_somd(self._work_dir, self._dir0, self._dir1, self._is_dual)
        elif self._engine == "GROMACS":
            return super()._analyse_gromacs(self._work_dir, self._dir0, self._dir1, self._is_dual)

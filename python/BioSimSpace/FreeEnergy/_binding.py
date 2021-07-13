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
Functionality for running binding free energy calculations.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Binding"]

import warnings as _warnings

from BioSimSpace._SireWrappers import System as _System
from BioSimSpace import Solvent as _Solvent
from BioSimSpace import Types as _Types
from BioSimSpace import Units as _Units

from . import _free_energy

class Binding(_free_energy.FreeEnergy):
    """A class for configuring and running binding free energy simulations."""

    def __init__(self, system0, system1=None, protocol=None, work_dir=None,
            engine=None, setup_only=False, property_map={},
            ignore_warnings=False, show_errors=True):
        """Constructor.

           Parameters
           ----------

           system0 : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the bound leg.

           system1 : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the free leg. If None, then the free
               leg is omitted.

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`, \
                     [:class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`,
                      :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`]
               The simulation protocol. If one obect is passed, then the same
               protocol will be used for both legs of the simulation. Passing
               two objects enables a different protocol for each leg, e.g. a
               different lambda schedule, or run time.

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

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

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

        if type(system0) is not _System:
            raise TypeError("'system0' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of the bound system. (Used for the first leg.)
            self._system0 = system0.copy()

            # The system must have a single perturbable molecule.
            if system0.nPerturbableMolecules() != 1:
                raise ValueError("The bound system must contain a single perturbable molecule! "
                                 "Use the 'BioSimSpace.Align' package to map and merge molecules.")

            # The system must be solvated.
            if system0.nWaterMolecules() == 0:
                raise ValueError("The bound system must be solvated! Use the 'BioSimSpace.Solvent' "
                                 "package to solvate your system.")

            # There must be at least one additional molecule in the system.
            if system0.nMolecules() == (system0.nWaterMolecules() + system0.nPerturbableMolecules()):
                raise ValueError("The bound system does not contain a protein molecule!")

        if system1 is not None:
            if type(system1) is not _System:
                raise TypeError("'system1' must be of type 'BioSimSpace._SireWrappers.System'")
            else:
                # Store a copy of the bound system. (Used for the first leg.)
                self._system1 = system1.copy()

                # The system must have a single perturbable molecule.
                if system1.nPerturbableMolecules() != 1:
                    raise ValueError("The free system must contain a single perturbable molecule! "
                                     "Use the 'BioSimSpace.Align' package to map and merge molecules.")

                # The system must be solvated.
                if system1.nWaterMolecules() == 0:
                    raise ValueError("The free system must be solvated! Use the 'BioSimSpace.Solvent' "
                                     "package to solvate your system.")
        else:
            self._is_dual = False
            self._system1 = system1

        # Initialise the process runner with all of the simulations required
        # for each leg.
        self._initialise_runner(self._system0, self._system1)

    def analyse(self):
        """Analyse the binding free energy data.

           Returns
           -------

           pmf_free : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the free leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf_bound : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the bound leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The binding free energy difference and its associated error.

           overlap_free : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the free leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.

           overlap_bound : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the bound leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.
        """

        # This method is just a wrapper to provide simulation specific doc
        # strings. We just call the base class method, which is aware of
        # the simulation type.
        if self._engine == "SOMD":
            return super()._analyse_somd(self._work_dir, self._dir0, self._dir1, self._is_dual)
        elif self._engine == "GROMACS":
            return super()._analyse_gromacs(self._work_dir, self._dir0, self._dir1, self._is_dual)

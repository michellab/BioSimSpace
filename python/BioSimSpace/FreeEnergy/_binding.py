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
Functionality for running binding free energy calculations.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Binding"]

import warnings as _warnings

from BioSimSpace._SireWrappers import System as _System
from BioSimSpace import Solvent as _Solvent
from BioSimSpace import Types as _Types
from BioSimSpace import Units as _Units

from . import _free_energy

class Binding(_free_energy.FreeEnergy):
    """A class for configuring and running binding free energy simulations."""

    def __init__(self, system, protocol=None, box=None, work_dir=None, engine=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`
               The simulation protocol.

           box : [:class:`Length <BioSimSpace.Types.Length>`]
               A list containing the box size in each dimension. This box will be
               used for the "free" leg of the simulation, which typically can be
               run with a significantly smaller box than the "bound" leg.

           work_dir : str
               The working directory for the simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "GROMACS", or "SOMD". If this argument is omitted then
               BioSimSpace will choose an appropriate engine for you.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(protocol, work_dir, engine)

        # Validate the input.

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of the bound system. (Used for the first leg.)
            self._system0 = system.copy()

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError("The system must contain a single perturbable molecule!")

            # The system must be solvated.
            if system.nWaterMolecules() == 0:
                raise ValueError("The system must be solvated!")

            # There must be at least one additional molecule in the system.
            if system.nMolecules() == (system.nWaterMolecules() + system.nPerturbableMolecules()):
                raise ValueError("The system does not contain a protein molecule!")

            # Extract the perturbable molecule.
            molecule = system.getPerturbableMolecules()[0]

            # Use the box size of the original system.
            if box is None:
                try:
                    prop = property_map.get("space", "space")
                    box = system._sire_object.property(prop).dimensions()
                    box = [_Units.Length.angstrom * x for x in box]
                except:
                    raise ValueError("The solvated protein-ligand system has no box information!")

            # Solvate using the user specified box.
            else:
                if len(box) != 3:
                    raise ValueError("The 'box' must have x, y, and z size information.")
                else:
                    if not all(isinstance(x, _Types.Length) for x in box):
                        raise ValueError("The box dimensions must be of type 'BioSimSpace.Types.Length'")

            # Try to get the water model used to solvate the system.
            try:
                water_model = system._sire_object.property("water_model").toString()
            # If the system wasn't solvated by BioSimSpace, e.g. read from file, then try
            # to guess the water model from the topology.
            except:
                num_point = self._system0.getWaterMolecules()[0].nAtoms()

                if num_point == 3:
                    # TODO: Assume TIP3P. Not sure how to detect SPC/E.
                    water_model = "tip3p"
                elif num_point == 4:
                    water_model = "tip4p"
                elif num_point == 5:
                    water_model = "tip5p"
                else:
                    raise RuntimeError("Unsupported %d-point water model!" % num_point)

                # Warn the user that we've guessed the water topology.
                _warnings.warn("Guessed water topology: %r" % water_model)

            # Solvate the perturbable molecule using the same water model as
            # the original system. (This is used for the second leg.)
            self._system1 = _Solvent.solvate(water_model, molecule=molecule, box=box)

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

           pmf_vacuum : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the vacuum leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The binding free energy difference and its associated error.
        """

        # This method is just a wrapper to provide simulation specific doc
        # strings. We just call the base class method, which is aware of
        # the simulation type.
        if self._engine == "SOMD":
            return super()._analyse_somd()
        elif self._engine == "GROMACS":
            return super()._analyse_gromacs()

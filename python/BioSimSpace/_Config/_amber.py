######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
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

"""Functionality for generating configuration files for AMBER."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Amber"]

import math as _math

from sire.legacy import Units as _SireUnits

from .. import Protocol as _Protocol

from ._config import Config as _Config


class Amber(_Config):
    """A class for generating configuration files for AMBER."""

    def __init__(
        self,
        system,
        protocol,
        is_pmemd=False,
        extra_options={},
        extra_lines=[],
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the process.

        is_pmemd : bool
            Whether the configuration is for a simulation using PMEMD.

        extra_options : dict
            A dictionary containing extra options. Overrides the defaults generated
            by the protocol.

        extra_lines : [str]
            A list of extra lines to put at the end of the configuration file.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(
            system,
            protocol,
            extra_options=extra_options,
            extra_lines=extra_lines,
            property_map=property_map,
        )

        if not isinstance(is_pmemd, bool):
            raise TypeError("'is_pmemd' must be of type 'bool'.")

        self._is_pmemd = is_pmemd

    def createConfig(self):
        """
        Create the list of configuration strings.

        Returns
        -------

        config : [str]
            The list of AMBER format configuration strings.
        """

        # Initialise the protocol lines.
        protocol_lines = []

        # Define some miscellaneous defaults.
        protocol_dict = {
            # Interval between reporting energies.
            "ntpr": 200,
            # Interval between saving restart files.
            "ntwr": self.restartInterval(),
            # Trajectory sampling frequency.
            "ntwx": self.restartInterval(),
            # Output coordinates as NetCDF.
            "ntxo": 2,
            # Whether to restart.
            "irest": int(self.isRestart()),
        }

        # Input.
        if self.isRestart():
            # Read coordinates and velocities.
            protocol_dict["ntx"] = 5
        else:
            # Only read coordinates from file.
            protocol_dict["ntx"] = 1

        # Minimisation.
        if isinstance(self._protocol, _Protocol.Minimisation):
            # Work out the number of steepest descent cycles.
            # This is 1000 or 10% of the number of steps, whichever is larger.
            if self.steps() <= 1000:
                num_steep = self.steps()
            else:
                num_steep = _math.ceil(self.steps() / 10)
                if num_steep < 1000:
                    num_steep = 1000

            # Minimisation simulation.
            protocol_dict["imin"] = 1
            # Set the minimisation method to XMIN.
            protocol_dict["ntmin"] = 2
            # Set the number of integration steps.
            protocol_dict["maxcyc"] = self.steps()
            # Set the number of steepest descent steps.
            protocol_dict["ncyc"] = num_steep
        else:
            # Define the timestep
            timestep = self._protocol.getTimeStep().picoseconds().value()
            # Set the integration time step.
            protocol_dict["dt"] = f"{timestep:.3f}"
            # Number of integration steps.
            protocol_dict["nstlim"] = self.steps()

        # Constraints.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            # Enable SHAKE.
            protocol_dict["ntc"] = 2
            # Don't calculate forces for constrained bonds.
            protocol_dict["ntf"] = 2

        # PBC.
        if not self.hasBox() or not self.hasWater():
            # No periodic box.
            protocol_dict["ntb"] = 0
            # Non-bonded cut-off.
            protocol_dict["cut"] = "999."
            if self._is_pmemd:
                # Use vacuum generalised Born model.
                self.addToConfig("  igb=6,")
        else:
            # Non-bonded cut-off.
            protocol_dict["cut"] = "8.0"
            # Wrap the coordinates.
            protocol_dict["iwrap"] = 1

        # Restraints.
        if isinstance(self._protocol, _Protocol.Equilibration):
            # Get the restraint.
            restraint = self._protocol.getRestraint()

            if restraint is not None:
                # Get the indices of the atoms that are restrained.
                if type(restraint) is str:
                    atom_idxs = self.system.getRestraintAtoms(restraint)
                else:
                    atom_idxs = restraint

                # Don't add restraints if there are no atoms to restrain.
                if len(atom_idxs) > 0:
                    # Generate the restraint mask based on atom indices.
                    restraint_mask = self._amber_mask_from_indices(
                        [i + 1 for i in atom_idxs]
                    )

                    # The restraintmask cannot be more than 256 characters.
                    if len(restraint_mask) > 256:

                        # AMBER has a limit on the length of the restraintmask
                        # so it's easy to overflow if we are matching by index
                        # on a large protein. As such, handle "backbone" and
                        # "heavy" restraints using a non-interoperable name mask.
                        if type(restraint) is str:
                            if restraint == "backbone":
                                restraint_mask = "@CA,C,O,N"
                            elif restraint == "heavy":
                                restraint_mask = "!:WAT & !@H"
                            elif restraint == "all":
                                restraint_mask = "!:WAT"

                        # We can't do anything about a custom restraint, since we don't
                        # know anything about the atoms.
                        else:
                            raise ValueError(
                                "AMBER atom 'restraintmask' exceeds 256 character limit!"
                            )

                    protocol_dict["ntr"] = 1
                    force_constant = self._protocol.getForceConstant()._sire_unit
                    force_constant = force_constant.to(
                        _SireUnits.kcal_per_mol / _SireUnits.angstrom2
                    )
                    protocol_dict["restraint_wt"] = force_constant
                    protocol_dict["restraintmask"] = f'"{restraint_mask}"'

        # Pressure control.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self.hasBox() and self.hasWater():
                    # Isotropic pressure scaling.
                    protocol_dict["ntp"] = 1
                    # Pressure in bar.
                    protocol_dict[
                        "pres0"
                    ] = f"{self._protocol.getPressure().bar().value():.5f}"
                    if isinstance(self._protocol, _Protocol.Equilibration):
                        # Berendsen barostat.
                        protocol_dict["barostat"] = 1
                    else:
                        # Monte Carlo barostat.
                        protocol_dict["barostat"] = 2
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation"
                    )

        # Temperature control.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            # Langevin dynamics.
            protocol_dict["ntt"] = 3
            # Collision frequency (ps).
            protocol_dict["gamma_ln"] = 2

            if isinstance(self._protocol, _Protocol.Equilibration):
                temp0 = self._protocol.getStartTemperature().kelvin().value()
                temp1 = self._protocol.getEndTemperature().kelvin().value()
                if not self._protocol.isConstantTemp():
                    # Initial temperature.
                    protocol_dict["tempi"] = f"{temp0:.2f}"
                    # Final temperature.
                    protocol_dict["temp0"] = f"{temp1:.2f}"
                    protocol_dict["nmropt"] = 1
                    protocol_lines += [
                        f"&wt TYPE='TEMP0', istep1=0, istep2={self.steps()}, value1={temp0:.2f}, value2={temp1:.2f} /"
                    ]
                else:
                    if not self.isRestart():
                        # Initial temperature.
                        protocol_dict["tempi"] = f"{temp0:.2f}"
                    # Constant temperature.
                    protocol_dict["temp0"] = f"{temp0:.2f}"
            else:
                temp = self._protocol.getTemperature().kelvin().value()
                if not self.isRestart():
                    # Initial temperature.
                    protocol_dict["tempi"] = f"{temp:.2f}"
                # Final temperature.
                protocol_dict["temp0"] = f"{temp:.2f}"

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **self._extra_options}
        dict_lines = [self._protocol.__class__.__name__, "&cntrl"]
        dict_lines += [
            f"   {k}={v}," for k, v in total_dict.items() if v is not None
        ] + ["/"]
        total_lines = protocol_lines + self._extra_lines
        if total_lines:
            total_lines += ["&wt TYPE='END' /"]
        total_lines = dict_lines + total_lines

        return total_lines

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

import itertools as _it
import math as _math
import warnings as _warnings

from .. import Protocol as _Protocol
from ..Align._merge import _squash
from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Units.Time import nanosecond as _nanosecond

from sire.legacy import Units as _SireUnits

from .. import Protocol as _Protocol
from ..Protocol._position_restraint_mixin import _PositionRestraintMixin
from ..Protocol._hmr_mixin import _HmrMixin

from ._config import Config as _Config


class Amber(_Config):
    """A class for generating configuration files for AMBER."""

    def __init__(self, system, protocol, property_map={}):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the process.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(system, protocol, property_map=property_map)

    def createConfig(
        self, version=None, is_pmemd=False, extra_options={}, extra_lines=[]
    ):
        """
        Create the list of configuration strings.

        version : float
            The AMBER version.

        is_pmemd : bool
            Whether the configuration is for a simulation using PMEMD.

        extra_options : dict
            A dictionary containing extra options. Overrides the defaults generated
            by the protocol.

        extra_lines : [str]
            A list of extra lines to put at the end of the configuration file.

        Returns
        -------

        config : [str]
            The list of AMBER format configuration strings.
        """

        # Validate input.

        if version and not isinstance(version, float):
            raise TypeError("'version' must be of type 'float'.")

        if not isinstance(is_pmemd, bool):
            raise TypeError("'is_pmemd' must be of type 'bool'.")

        if not isinstance(extra_options, dict):
            raise TypeError("'extra_options' must be of type 'dict'.")
        else:
            keys = extra_options.keys()
            if not all(isinstance(k, str) for k in keys):
                raise TypeError("Keys of 'extra_options' must be of type 'str'.")

        if not isinstance(extra_lines, list):
            raise TypeError("'extra_lines' must be of type 'list'.")
        else:
            if not all(isinstance(line, str) for line in extra_lines):
                raise TypeError("Lines in 'extra_lines' must be of type 'str'.")

        # Initialise the protocol lines.
        protocol_lines = []

        # Define some miscellaneous defaults.
        protocol_dict = {
            # Interval between reporting energies.
            "ntpr": self.reportInterval(),
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
            if isinstance(self._protocol, _Protocol._FreeEnergyMixin):
                protocol_dict["ntmin"] = 2 # if using sc potentials (ifsc=1) as with free energies, must be steepest descent
            else:
                protocol_dict["ntmin"] = 1
            # Set the number of integration steps.
            protocol_dict["maxcyc"] = self.steps()
            # Set the number of steepest descent steps.
            protocol_dict["ncyc"] = num_steep
            # Report energies every 100 steps.
            protocol_dict["ntpr"] = 100
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

        # Periodic boundary conditions.
        if not self.hasBox() or not self.hasWater():
            # No periodic box.
            protocol_dict["ntb"] = 0
            # Non-bonded cut-off.
            protocol_dict["cut"] = "999."
            if is_pmemd:
                # Use vacuum generalised Born model.
                self.addToConfig("  igb=6,")
        else:
            # Non-bonded cut-off.
            protocol_dict["cut"] = "10.0"
            # Wrap the coordinates.
            protocol_dict["iwrap"] = 1

        # Postion restraints.
        if isinstance(self._protocol, _PositionRestraintMixin):
            # Get the restraint.
            restraint = self._protocol.getRestraint()

            if restraint is not None:
                # Get the indices of the atoms that are restrained.
                if type(restraint) is str:
                    atom_idxs = self._system.getRestraintAtoms(restraint)
                else:
                    atom_idxs = restraint

                # Don't add restraints if there are no atoms to restrain.
                if len(atom_idxs) > 0:
                    # Generate the restraint mask based on atom indices.
                    restraint_mask = self._create_restraint_mask(atom_idxs)

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
        if not isinstance(
            self._protocol,
            (_Protocol.Metadynamics, _Protocol.Steering, _Protocol.Minimisation),
        ):
            # Langevin dynamics.
            protocol_dict["ntt"] = 3
            # Collision frequency (1 / ps).
            protocol_dict["gamma_ln"] = "{:.5f}".format(
                1 / self._protocol.getThermostatTimeConstant().picoseconds().value()
            )

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


        # Free energies.
        if isinstance(self._protocol, _Protocol._FreeEnergyMixin):
            # Free energy mode.
            protocol_dict["icfe"] = 1
            # Use softcore potentials.
            protocol_dict["ifsc"] = 1
            # Remove SHAKE constraints.
            protocol_dict["ntf"] = 1
            protocol = [str(x) for x in self._protocol.getLambdaValues()]
            # Number of lambda values.
            protocol_dict["mbar_states"] = len(protocol)
            # Lambda values.
            protocol_dict["mbar_lambda"] = ", ".join(protocol)
            # Current lambda value.
            protocol_dict["clambda"] = self._protocol.getLambda()
            if isinstance(self._protocol, _Protocol.Production):
                # Calculate MBAR energies.
                protocol_dict["ifmbar"] = 1
                # Output dVdl
                protocol_dict["logdvdl"] = 1
            # Atom masks based on if HMR.
            hmr = self._protocol.getHmr()
            protocol_dict = {**protocol_dict, **
                             self._generate_amber_fep_masks(hmr)}
            # Do not wrap the coordinates for equilibration steps 
            # as this can create issues for the free leg of FEP runs.
            if isinstance(self._protocol, _Protocol.Equilibration) or isinstance(self._protocol, _Protocol.Minimisation):
                protocol_dict["iwrap"] = 0
                
        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        dict_lines = [self._protocol.__class__.__name__, "&cntrl"]
        dict_lines += [
            f"   {k}={v}," for k, v in total_dict.items() if v is not None
        ] + ["/"]
        total_lines = protocol_lines + extra_lines
        if total_lines:
            total_lines += ["&wt TYPE='END' /"]
        total_lines = dict_lines + total_lines

        return total_lines

    def _create_restraint_mask(self, atom_idxs):
        """
        Internal helper function to create an AMBER restraint mask from a
        list of atom indices.

        Parameters
        ----------

        atom_idxs : [int]
            A list of atom indices.

        Returns
        -------

        restraint_mask : str
            The AMBER restraint mask.
        """

        if not isinstance(atom_idxs, (list, tuple)):
            raise TypeError("'atom_idxs' must be a list of 'int' types.")

        if not all(type(x) is int for x in atom_idxs):
            raise TypeError("'atom_idxs' must be a list of 'int' types.")

        # AMBER has a restriction on the number of characters in the restraint
        # mask (not documented) so we can't just use comma-separated atom
        # indices. Instead we loop through the indices and use hyphens to
        # separate contiguous blocks of indices, e.g. 1-23,34-47,...

        # Create a set to sort and ensure no duplicates, then convert back to a list.
        # This should already by done, but do so again in case the user is accessing
        # the method directly.
        atom_idxs = list(set(atom_idxs))
        atom_idxs.sort()

        # Handle single atom restraints differently.
        if len(atom_idxs) == 1:
            restraint_mask = f"@{atom_idxs[0]+1}"

        else:
            # Start the mask with the first atom index. (AMBER is 1 indexed.)
            restraint_mask = f"@{atom_idxs[0]+1}"

            # Store the current index.
            prev_idx = atom_idxs[0]

            # Store the lead index for this block.
            lead_idx = prev_idx

            # Loop over all other indices.
            for idx in atom_idxs[1:]:
                # There is a gap in the indices.
                if idx - prev_idx > 1:
                    if prev_idx != lead_idx:
                        restraint_mask += f"{prev_idx+1},{idx+1}"
                    else:
                        restraint_mask += f",{idx+1}"
                    lead_idx = idx
                else:
                    # This is the first index beyond the lead.
                    if idx - lead_idx == 1:
                        restraint_mask += "-"
                # Set the value of the previous index.
                prev_idx = idx

            # Add the final atom to the mask.
            if idx - atom_idxs[-2] == 1:
                restraint_mask += f"{idx+1}"
            else:
                if idx != lead_idx:
                    restraint_mask += f",{idx+1}"

        return restraint_mask


    def _amber_mask_from_indices(self, atom_idxs):
        """Internal helper function to create an AMBER restraint mask from a
           list of atom indices.

           Parameters
           ----------

           atom_idxs : [int]
               A list of atom indices.

           Returns
           -------

           restraint_mask : str
               The AMBER restraint mask.
        """
        # AMBER has a restriction on the number of characters in the restraint
        # mask (not documented) so we can't just use comma-separated atom
        # indices. Instead we loop through the indices and use hyphens to
        # separate contiguous blocks of indices, e.g. 1-23,34-47,...

        if atom_idxs:
            atom_idxs = sorted(list(set(atom_idxs)))
            if not all(isinstance(x, int) for x in atom_idxs):
                raise TypeError("'atom_idxs' must be a list of 'int' types.")
            groups = []
            initial_idx = atom_idxs[0]
            for prev_idx, curr_idx in _it.zip_longest(atom_idxs, atom_idxs[1:]):
                if curr_idx != prev_idx + 1 or curr_idx is None:
                    if initial_idx == prev_idx:
                        groups += [str(initial_idx)]
                    else:
                        groups += [f"{initial_idx}-{prev_idx}"]
                    initial_idx = curr_idx
            mask = "@" + ",".join(groups)
        else:
            mask = ""

        return mask

    def _generate_amber_fep_masks(self, hmr):
        """Internal helper function which generates timasks and scmasks based on the system.

           Parameters
           ----------

           hmr : bool
               Whether HMR has been applied or not.
               Generates a different mask based on this.

           Returns
           -------

           option_dict : dict
               A dictionary of AMBER-compatible options.
        """

        # Squash the system into an AMBER-friendly format.
        squashed_system = _squash(self._system)

        # When HMR is used, there can be no noshakemask
        HMR_on = hmr

        # Extract all perturbable molecules from both the squashed and the merged systems.
        perturbable_mol_mask = []
        mols_hybr = []
        nondummy_indices0, nondummy_indices1 = [], []
        for mol in self._system.getMolecules():
            if mol._is_perturbable:
                perturbable_mol_mask += [0, 1]
                mols_hybr += [mol]
                nondummy_indices0 += [[atom.index() for atom in mol.getAtoms()
                                       if "du" not in atom._sire_object.property("ambertype0")]]
                nondummy_indices1 += [[atom.index() for atom in mol.getAtoms()
                                       if "du" not in atom._sire_object.property("ambertype1")]]
            else:
                perturbable_mol_mask += [None]
        mols0 = [squashed_system.getMolecule(i) for i, mask in enumerate(
            perturbable_mol_mask) if mask == 0]
        mols1 = [squashed_system.getMolecule(i) for i, mask in enumerate(
            perturbable_mol_mask) if mask == 1]

        # Find the perturbed atom indices withing the squashed system.
        mols0_indices = [squashed_system.getIndex(
            atom) + 1 for mol in mols0 for atom in mol.getAtoms()]
        mols1_indices = [squashed_system.getIndex(
            atom) + 1 for mol in mols1 for atom in mol.getAtoms()]

        # Find the dummy indices within the squashed system.
        offsets = [0] + list(_it.accumulate(mol.nAtoms()
                             for mol in squashed_system.getMolecules()))
        offsets0 = [offsets[i]
                    for i, mask in enumerate(perturbable_mol_mask) if mask == 0]
        offsets1 = [offsets[i]
                    for i, mask in enumerate(perturbable_mol_mask) if mask == 1]
        dummy0_indices = [offset + idx_map.index(atom.index()) + 1
                          for mol, offset, idx_map in zip(mols_hybr, offsets0, nondummy_indices0)
                          for atom in mol.getAtoms()
                          if "du" in atom._sire_object.property("ambertype1")]
        dummy1_indices = [offset + idx_map.index(atom.index()) + 1
                          for mol, offset, idx_map in zip(mols_hybr, offsets1, nondummy_indices1)
                          for atom in mol.getAtoms()
                          if "du" in atom._sire_object.property("ambertype0")]

        # If it is HMR
        if HMR_on == True:
            no_shake_mask = ""
        else:
            no_shake_mask = self._amber_mask_from_indices(
                mols0_indices + mols1_indices)

        # Create an option dict with amber masks generated from the above indices.
        option_dict = {
            "timask1": f"\"{self._amber_mask_from_indices(mols0_indices)}\"",
            "timask2": f"\"{self._amber_mask_from_indices(mols1_indices)}\"",
            "scmask1": f"\"{self._amber_mask_from_indices(dummy0_indices)}\"",
            "scmask2": f"\"{self._amber_mask_from_indices(dummy1_indices)}\"",
            "noshakemask": f"\"{no_shake_mask}\"",
        }

        return option_dict

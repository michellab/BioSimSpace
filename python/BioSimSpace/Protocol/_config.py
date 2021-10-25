import itertools as _it
import math as _math
import warnings as _warnings

from BioSimSpace.Align._merge import _squash

from BioSimSpace import Protocol as _Protocol


class ConfigFactory:
    # TODO: Integrate this class better into the other Protocols.
    """A class for generating a config based on a template protocol."""

    def __init__(self, system, protocol):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol <BioSimSpace.Protocol>`
        """
        self.system = system
        self.protocol = protocol

    @property
    def _has_box(self):
        """Return whether the current system has a box."""
        if "space" in self.system._sire_object.propertyKeys():
            has_box = True
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False
        return has_box

    @property
    def _has_water(self):
        """Return whether the current system has any water molecules."""
        return self.system.nWaterMolecules() > 0

    @property
    def _report_interval(self):
        """Return the report interval based on the protocol value."""
        if isinstance(self.protocol, _Protocol.Minimisation):
            report_interval = 100
        else:
            report_interval = self.protocol.getReportInterval()
            if report_interval > self._steps:
                report_interval = self._steps
        return report_interval

    @property
    def _restart(self):
        """Return whether this is a restart simulation."""
        if isinstance(self.protocol, (_Protocol.Equilibration, _Protocol.Minimisation)):
            return False
        elif isinstance(self.protocol, _Protocol.FreeEnergy):
            return True
        else:
            return self.protocol.isRestart()

    @property
    def _restart_interval(self):
        """Return the restart interval based on the protocol value."""
        if isinstance(self.protocol, _Protocol.Minimisation):
            restart_interval = None
        else:
            restart_interval = self.protocol.getRestartInterval()
            if restart_interval > self._steps:
                restart_interval = self._steps
        return restart_interval

    @property
    def _steps(self):
        # Return the number of steps based on the protocol value.
        if isinstance(self.protocol, _Protocol.Minimisation):
            steps = self.protocol.getSteps()
        else:
            steps = _math.ceil(self.protocol.getRunTime() / self.protocol.getTimeStep())
        return steps

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

    def _generate_amber_fep_masks(self):
        """Internal helper function which generates timasks and scmasks based on the system.

           Returns
           -------

           option_dict : dict
               A dictionary of AMBER-compatible options.
        """

        # Squash the system into an AMBER-friendly format.
        squashed_system = _squash(self.system)

        # Extract all perturbable molecules from both the squashed and the merged systems.
        perturbable_mol_mask = []
        mols_hybr = []
        nondummy_indices0, nondummy_indices1 = [], []
        for mol in self.system.getMolecules():
            if mol._is_perturbable:
                perturbable_mol_mask += [0, 1]
                mols_hybr += [mol]
                nondummy_indices0 += [[atom.index() for atom in mol.getAtoms()
                                       if "du" not in atom._sire_object.property("ambertype0")]]
                nondummy_indices1 += [[atom.index() for atom in mol.getAtoms()
                                       if "du" not in atom._sire_object.property("ambertype1")]]
            else:
                perturbable_mol_mask += [None]
        mols0 = [squashed_system.getMolecule(i) for i, mask in enumerate(perturbable_mol_mask) if mask == 0]
        mols1 = [squashed_system.getMolecule(i) for i, mask in enumerate(perturbable_mol_mask) if mask == 1]

        # Find the perturbed atom indices withing the squashed system.
        mols0_indices = [squashed_system.getIndex(atom) + 1 for mol in mols0 for atom in mol.getAtoms()]
        mols1_indices = [squashed_system.getIndex(atom) + 1 for mol in mols1 for atom in mol.getAtoms()]

        # Find the dummy indices within the squashed system.
        offsets = [0] + list(_it.accumulate(mol.nAtoms() for mol in squashed_system.getMolecules()))
        offsets0 = [offsets[i] for i, mask in enumerate(perturbable_mol_mask) if mask == 0]
        offsets1 = [offsets[i] for i, mask in enumerate(perturbable_mol_mask) if mask == 1]
        dummy0_indices = [offset + idx_map.index(atom.index()) + 1
                          for mol, offset, idx_map in zip(mols_hybr, offsets0, nondummy_indices0)
                          for atom in mol.getAtoms()
                          if "du" in atom._sire_object.property("ambertype1")]
        dummy1_indices = [offset + idx_map.index(atom.index()) + 1
                          for mol, offset, idx_map in zip(mols_hybr, offsets1, nondummy_indices1)
                          for atom in mol.getAtoms()
                          if "du" in atom._sire_object.property("ambertype0")]

        # Create an option dict with amber masks generated from the above indices.
        option_dict = {
            "timask1": f"\"{self._amber_mask_from_indices(mols0_indices)}\"",
            "timask2": f"\"{self._amber_mask_from_indices(mols1_indices)}\"",
            "scmask1": f"\"{self._amber_mask_from_indices(dummy0_indices)}\"",
            "scmask2": f"\"{self._amber_mask_from_indices(dummy1_indices)}\"",
        }

        return option_dict

    def generateAmberConfig(self, extra_options=None, extra_lines=None):
        """Outputs the current protocol in a format compatible with AMBER.

           Parameters
           ----------

           extra_options : dict
               A dictionary containing extra options. Overrides the ones generated from the protocol.

           extra_lines : list
               A list of extra lines to be put at the end of the script.

           Returns
           -------

           config : list
               The generated config list in an AMBER format.
        """

        extra_options = extra_options if extra_options is not None else {}
        extra_lines = extra_lines if extra_lines is not None else []

        protocol_lines = []

        # Define some miscellaneous defaults.
        protocol_dict = {
            "ntpr": self._report_interval,          # Interval between reporting energies.
            "ntwr": self._restart_interval,         # Interval between saving restart files.
            "ntwx": self._restart_interval,         # Trajectory sampling frequency.
            "ntxo": 1,                              # Output coordinates in ASCII.
            "irest": int(self._restart),            # Whether to restart.
        }

        # Input.
        if self._restart:
            protocol_dict["ntx"] = 5                # Read coordinates and velocities.
        else:
            protocol_dict["ntx"] = 1                # Only read coordinates from file.

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            # Work out the number of steepest descent cycles.
            # This is 1000 or 10% of the number of steps, whichever is larger.
            if self._steps <= 1000:
                num_steep = self._steps
            else:
                num_steep = _math.ceil(self._steps / 10)
                if num_steep < 1000:
                    num_steep = 1000

            protocol_dict["imin"] = 1               # Minimisation simulation.
            protocol_dict["maxcyc"] = self._steps   # Set the number of steps.
            protocol_dict["ncyc"] = num_steep       # Set the number of steepest descent steps.
        else:
            protocol_dict["dt"] = f"{self.protocol.getTimeStep().picoseconds().magnitude():.3f}"  # Time step.
            protocol_dict["nstlim"] = self._steps   # Number of integration steps.

        # Constraints.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntc"] = 2                # Enable SHAKE.
            protocol_dict["ntf"] = 2                # Don't calculate forces for constrained bonds.

        # PBC.
        if not self._has_box or not self._has_water:
            protocol_dict["ntb"] = 0                # No periodic box.
            protocol_dict["cut"] = "999."           # Non-bonded cut-off.
        else:
            protocol_dict["cut"] = "8.0"            # Non-bonded cut-off.

        # Restraints.
        if isinstance(self.protocol, _Protocol.Equilibration):
            # Restrain the backbone.
            restraint = self.protocol.getRestraint()

            if restraint is not None:
                # Get the indices of the atoms that are restrained.
                if type(restraint) is str:
                    atom_idxs = self.system.getRestraintAtoms(restraint)
                else:
                    atom_idxs = restraint

                # Don't add restraints if there are no atoms to restrain.
                if len(atom_idxs) > 0:
                    # Generate the restraint mask based on atom indices.
                    restraint_mask = self._amber_mask_from_indices(atom_idxs)

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
                            raise ValueError("AMBER atom 'restraintmask' exceeds 256 character limit!")

                    protocol_dict["ntr"] = 1
                    protocol_dict["restraint_wt"] = 10
                    protocol_dict["restraintmask"] = restraint_mask

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    protocol_dict["ntp"] = 1        # Isotropic pressure scaling.
                    protocol_dict["pres0"] = f"{self.protocol.getPressure().bar().magnitude():.5f}"  # Pressure in bar.
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntt"] = 3                # Langevin dynamics.
            protocol_dict["gamma_ln"] = 2           # Collision frequency (ps).
            if isinstance(self.protocol, _Protocol.Equilibration):
                temp0 = self.protocol.getStartTemperature().kelvin().magnitude()
                temp1 = self.protocol.getEndTemperature().kelvin().magnitude()
                if not self.protocol.isConstantTemp():
                    protocol_dict["tempi"] = f"{temp0:.2f}"  # Initial temperature.
                    protocol_dict["temp0"] = f"{temp1:.2f}"  # Final temperature.
                    protocol_dict["nmropt"] = 1
                    protocol_lines += [
                        f"&wt TYPE='TEMP0', istep1=0, istep2={self._steps}, value1={temp0:.2f}, value2={temp1:.2f} /"
                    ]
                else:
                    protocol_dict["temp0"] = f"{temp0:.2f}"  # Constant temperature.
            else:
                temp = self.protocol.getTemperature().kelvin().magnitude()
                if not self._restart:
                    protocol_dict["tempi"] = f"{temp:.2f}"  # Initial temperature.
                protocol_dict["temp0"] = f"{temp:.2f}"      # Final temperature.

        # Free energies.
        if isinstance(self.protocol, _Protocol.FreeEnergy):
            protocol_dict["icfe"] = 1                                               # Free energy mode.
            protocol_dict["ifsc"] = 1                                               # Use softcore potentials.
            protocol_dict["ntf"] = 1                                                # Remove SHAKE constraints.
            protocol_dict["ifmbar"] = 1                                             # Calculate MBAR energies.
            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            protocol_dict["mbar_states"] = len(protocol)                            # Number of lambda values.
            protocol_dict["mbar_lambda"] = ", ".join(protocol)                      # Lambda values.
            protocol_dict["clambda"] = self.protocol.getLambda()                    # Current lambda value.
            protocol_dict = {**protocol_dict, **self._generate_amber_fep_masks()}   # Atom masks.

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        dict_lines = [self.protocol.__class__.__name__, "&cntrl"]
        dict_lines += [f"   {k}={v}," for k, v in total_dict.items() if v is not None] + ["/"]
        total_lines = protocol_lines + extra_lines
        if total_lines:
            total_lines += ["&wt TYPE='END' /"]
        total_lines = dict_lines + total_lines

        return total_lines
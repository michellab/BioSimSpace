import itertools as _it
import math as _math
import warnings as _warnings

from sire import units as _SireUnits

from ..Align._merge import _squash
from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Units.Time import nanosecond as _nanosecond

from .. import Protocol as _Protocol

class ConfigFactory:
    # TODO: Integrate this class better into the other Protocols.
    """A class for generating a config based on a template protocol."""

    def __init__(self, system, protocol, squashed_system=None):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol <BioSimSpace.Protocol>`

           squashed_system : :class:`System <BioSimSpace._SireWrappers.System>`
               (Optional) an AMBER-compatible squashed system.
        """
        self.system = system
        self.protocol = protocol
        self.squashed_system = squashed_system

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
        try:
            return self.protocol.isRestart()
        except:
            pass
        return False

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
            # AMBER masks are 1-indexed, while BioSimSpace indices are 0-indexed.
            atom_idxs = [x + 1 for x in sorted(list(set(atom_idxs)))]
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

    def _generate_amber_fep_masks(self, timestep):
        """Internal helper function which generates timasks and scmasks based on the system.

           Parameters
           ----------

           timestep : [float]
               The timestep in ps for the FEP perturbation. Generates a different mask based on this.

           Returns
           -------

           option_dict : dict
               A dictionary of AMBER-compatible options.
        """
        # Squash the system into an AMBER-friendly format.
        if self.squashed_system is None:
            self.squashed_system, _ = _squash(self.system)

        # Get the perturbed molecules and the corresponding squashed molecules.
        pertmols = self.system.getPerturbableMolecules()
        pertmol_offset = len(self.system) - len(pertmols)
        squashed_pertmols = self.squashed_system[pertmol_offset:]
        atom_offsets = [0] + list(_it.accumulate(mol.nAtoms() for mol in self.squashed_system.getMolecules()))

        # Find the perturbed atom indices withing the squashed system.
        mols0_indices, mols1_indices = [], []
        dummy0_indices, dummy1_indices = [], []
        for i, pertmol in enumerate(pertmols):
            mol0 = squashed_pertmols[2 * i]
            mol1 = squashed_pertmols[2 * i + 1]
            atom_offset0 = atom_offsets[pertmol_offset + 2 * i]
            atom_offset1 = atom_offset0 + mol0.nAtoms()
            mols0_indices += list(range(atom_offset0, atom_offset1))
            mols1_indices += list(range(atom_offset1, atom_offset1 + mol1.nAtoms()))
            nondummy_indices0 = [atom.index() for atom in pertmol.getAtoms()
                                 if "du" not in atom._sire_object.property("ambertype0")]
            nondummy_indices1 = [atom.index() for atom in pertmol.getAtoms()
                                 if "du" not in atom._sire_object.property("ambertype1")]
            dummy0_indices += [atom_offset0 + nondummy_indices0.index(atom.index())
                               for atom in pertmol.getAtoms()
                               if "du" in atom._sire_object.property("ambertype1")]
            dummy1_indices += [atom_offset1 + nondummy_indices1.index(atom.index())
                               for atom in pertmol.getAtoms()
                               if "du" in atom._sire_object.property("ambertype0")]

        # Define whether HMR is used based on the timestep.
        # When HMR is used, there can be no SHAKE.
        if timestep >= 0.004:
            no_shake_mask = ""
        else:
            no_shake_mask = self._amber_mask_from_indices(mols0_indices + mols1_indices)

        # Create an option dict with amber masks generated from the above indices.
        option_dict = {
            "timask1": f"\"{self._amber_mask_from_indices(mols0_indices)}\"",
            "timask2": f"\"{self._amber_mask_from_indices(mols1_indices)}\"",
            "scmask1": f"\"{self._amber_mask_from_indices(dummy0_indices)}\"",
            "scmask2": f"\"{self._amber_mask_from_indices(dummy1_indices)}\"",
            "noshakemask": f"\"{no_shake_mask}\"",
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
            "ntpr": 200,                            # Interval between reporting energies.
            "ntwr": self._restart_interval,         # Interval between saving restart files.
            "ntwx": self._restart_interval,         # Trajectory sampling frequency.
            "ntxo": 2,                              # Output coordinates as NetCDF.
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
            protocol_dict["ntmin"] = 2              # Set the minimisation method to XMIN
            protocol_dict["maxcyc"] = self._steps   # Set the number of steps.
            protocol_dict["ncyc"] = num_steep       # Set the number of steepest descent steps.
            # FIX need to remove and fix this, only for initial testing
            timestep = 0.004
        else:
            # Define the timestep
            timestep = self.protocol.getTimeStep().picoseconds().value() # Get the time step in ps
            protocol_dict["dt"] = f"{timestep:.3f}" # Time step.
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
            protocol_dict["iwrap"] = 1              # Wrap the coordinates.

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
                    restraint_mask = self._amber_mask_from_indices([i + 1 for i in atom_idxs])

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
                    force_constant = self.protocol.getForceConstant()._sire_unit
                    force_constant = force_constant.to(_SireUnits.kcal_per_mol / _SireUnits.angstrom2)
                    protocol_dict["restraint_wt"] = force_constant
                    protocol_dict["restraintmask"] = f"\"{restraint_mask}\""

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    protocol_dict["ntp"] = 1        # Isotropic pressure scaling.
                    protocol_dict["pres0"] = f"{self.protocol.getPressure().bar().value():.5f}"  # Pressure in bar.
                    if isinstance(self.protocol, _Protocol.Equilibration):
                        protocol_dict["barostat"] = 1         # Berendsen barostat.
                    else:
                        protocol_dict["barostat"] = 2         # Monte Carlo barostat.
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntt"] = 3                # Langevin dynamics.
            protocol_dict["gamma_ln"] = 2           # Collision frequency (ps).
            if isinstance(self.protocol, _Protocol.Equilibration):
                temp0 = self.protocol.getStartTemperature().kelvin().value()
                temp1 = self.protocol.getEndTemperature().kelvin().value()
                if not self.protocol.isConstantTemp():
                    protocol_dict["tempi"] = f"{temp0:.2f}"  # Initial temperature.
                    protocol_dict["temp0"] = f"{temp1:.2f}"  # Final temperature.
                    protocol_dict["nmropt"] = 1
                    protocol_lines += [
                        f"&wt TYPE='TEMP0', istep1=0, istep2={self._steps}, value1={temp0:.2f}, value2={temp1:.2f} /"
                    ]
                else:
                    if not self._restart:
                        protocol_dict["tempi"] = f"{temp0:.2f}"  # Initial temperature.
                    protocol_dict["temp0"] = f"{temp0:.2f}"  # Constant temperature.
            else:
                temp = self.protocol.getTemperature().kelvin().value()
                if not self._restart:
                    protocol_dict["tempi"] = f"{temp:.2f}"   # Initial temperature.
                protocol_dict["temp0"] = f"{temp:.2f}"       # Final temperature.

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            protocol_dict["icfe"] = 1                                               # Free energy mode.
            protocol_dict["ifsc"] = 1                                               # Use softcore potentials.
            protocol_dict["ntf"] = 1                                                # Remove SHAKE constraints.
            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            protocol_dict["mbar_states"] = len(protocol)                            # Number of lambda values.
            protocol_dict["mbar_lambda"] = ", ".join(protocol)                      # Lambda values.
            protocol_dict["clambda"] = self.protocol.getLambda()                    # Current lambda value.
            if isinstance(self.protocol, _Protocol.Production):
                protocol_dict["ifmbar"] = 1                                         # Calculate MBAR energies.
                protocol_dict["logdvdl"] = 1                                        # Output dVdl
            protocol_dict = {**protocol_dict, **self._generate_amber_fep_masks(timestep)}   # Atom masks.

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        dict_lines = [self.protocol.__class__.__name__, "&cntrl"]
        dict_lines += [f"   {k}={v}," for k, v in total_dict.items() if v is not None] + ["/"]
        total_lines = protocol_lines + extra_lines
        if total_lines:
            total_lines += ["&wt TYPE='END' /"]
        total_lines = dict_lines + total_lines

        return total_lines

    def generateGromacsConfig(self, extra_options=None, extra_lines=None):
        """Outputs the current protocol in a format compatible with GROMACS.

        Parameters
        ----------

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        Returns
        -------

        config : list
            The generated config list in a GROMACS format.
        """

        extra_options = extra_options if extra_options is not None else {}
        extra_lines = extra_lines if extra_lines is not None else []

        # Define some miscellaneous defaults.
        protocol_dict = {
            "nstlog": self._report_interval,                                        # Interval between writing to the log file.
            "nstenergy": self._restart_interval,                                    # Interval between writing to the energy file.
            "nstxout-compressed": self._restart_interval,                                      # Interval between writing to the trajectory file.
        }

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["integrator"] = "steep"                                   # Minimisation simulation.
        else:
            timestep = self.protocol.getTimeStep().picoseconds().value()        # Define the timestep in picoseconds
            protocol_dict["dt"] = f"{timestep:.3f}"                                 # Integration time step.
        protocol_dict["nsteps"] = self._steps                                       # Number of integration steps.

        # Constraints.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["constraints"] = "h-bonds"                            # Rigid bonded hydrogens.
            protocol_dict["constraint-algorithm"] = "LINCS"                         # Linear constraint solver.

        # PBC.
        protocol_dict["pbc"] = "xyz"                                                # Simulate a fully periodic box.
        protocol_dict["cutoff-scheme"] = "Verlet"                                   # Use Verlet pair lists.
        if self._has_box and self._has_water:
            protocol_dict["ns-type"] = "grid"                                       # Use a grid to search for neighbours.
            protocol_dict["nstlist"] = "20"                                         # Rebuild neighbour list every 20 steps. Recommended in the manual for parallel simulations and/or non-bonded force calculation on the GPU.
            protocol_dict["rlist"] = "0.8"                                          # Set short-range cutoff.
            protocol_dict["rvdw"] = "0.8"                                           # Set van der Waals cutoff.
            protocol_dict["rcoulomb"] = "0.8"                                       # Set Coulomb cutoff.
            protocol_dict["coulombtype"] = "PME"                                    # Fast smooth Particle-Mesh Ewald.
            protocol_dict["DispCorr"] = "EnerPres"                                  # Dispersion corrections for energy and pressure.
        else:
            # Perform vacuum simulations by implementing pseudo-PBC conditions,
            # i.e. run calculation in a near-infinite box (333.3 nm).
            # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588
            protocol_dict["nstlist"] = "1"                                          # Single neighbour list (all particles interact).
            protocol_dict["rlist"] = "333.3"                                        # "Infinite" short-range cutoff.
            protocol_dict["rvdw"] = "333.3"                                         # "Infinite" van der Waals cutoff.
            protocol_dict["rcoulomb"] = "333.3"                                     # "Infinite" Coulomb cutoff.
            protocol_dict["coulombtype"] = "Cut-off"                                # Plain cut-off.
        protocol_dict["vdwtype"] = "Cut-off"                                        # Twin-range van der Waals cut-off.

        # Restraints.
        if isinstance(self.protocol, _Protocol.Equilibration):
            protocol_dict["refcoord-scaling"] = "com"                               # The actual restraints need to be defined elsewhere.

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    protocol_dict["pcoupl"] = "c-rescale"                       # Barostat type.
                    protocol_dict["tau-p"] = 1                                      # 1ps time constant for pressure coupling.
                    protocol_dict["ref-p"] = f"{self.protocol.getPressure().bar().value():.5f}"  # Pressure in bar.
                    protocol_dict["compressibility"] = "4.5e-5"                     # Compressibility of water.
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["integrator"] = "md"                                      # leap-frog dynamics.
            protocol_dict["tcoupl"] = "v-rescale"
            protocol_dict["tc-grps"] = "system"                                     # A single temperature group for the entire system.
            protocol_dict["tau-t"] = 2                                              # Collision frequency (ps).

            if isinstance(self.protocol, _Protocol.Equilibration):
                if self.protocol.isConstantTemp():
                    temp = "%.2f" % self.protocol.getStartTemperature().kelvin().value()
                    protocol_dict["ref-t"] = temp
                    protocol_dict["gen-vel"] = "yes"
                    protocol_dict["gen-temp"] = temp
                else:
                    #still need a reference temperature for each group, even when heating/cooling
                    protocol_dict["ref-t"] = "%.2f" % self.protocol.getEndTemperature().kelvin().value()
                    # Work out the final time of the simulation.
                    timestep = self.protocol.getTimeStep().picoseconds().value()
                    end_time = _math.floor(timestep * self._steps)

                    protocol_dict["annealing"] = "single"                               # Single sequence of annealing points.
                    protocol_dict["annealing-npoints"] = 2                              # Two annealing points for "system" temperature group.

                    # Linearly change temperature between start and end times.
                    protocol_dict["annealing-time"] = "0 %d" % end_time
                    protocol_dict["annealing-temp"] = "%.2f %.2f" % (
                        self.protocol.getStartTemperature().kelvin().value(),
                        self.protocol.getEndTemperature().kelvin().value(),
                    )
            else:
                protocol_dict["ref-t"] = "%.2f" % self.protocol.getTemperature().kelvin().value()

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            protocol_dict["free-energy"] = "yes"                                    # Free energy mode.
            nDecoupledMolecules = self.system.nDecoupledMolecules()
            if nDecoupledMolecules == 1:
                [mol, ] = self.system.getDecoupledMolecules()
                decouple_dict = mol._sire_object.property("decouple")
                protocol_dict["couple-moltype"] = mol._sire_object.name().value()
                def tranform(charge, LJ):
                    if charge and LJ:
                        return 'vdw-q'
                    elif charge and not LJ:
                        return 'q'
                    elif not charge and LJ:
                        return 'vdw'
                    else:
                        return 'none'
                protocol_dict["couple-lambda0"] = tranform(
                    decouple_dict["charge"][0],
                    decouple_dict["LJ"][0]
                )
                protocol_dict["couple-lambda1"] = tranform(
                    decouple_dict["charge"][1],
                    decouple_dict["LJ"][1]
                )
                # Add the soft-core parameters for the ABFE
                protocol_dict['sc-alpha'] = 0.5
                protocol_dict['sc-power'] = 1
                protocol_dict['sc-sigma'] = 0.3
                if decouple_dict['intramol'].value():
                    # The intramol is being coupled to the lambda change and thus being annihilated.
                    protocol_dict["couple-intramol"] = 'yes'
                else:
                    protocol_dict["couple-intramol"] = 'no'
            elif nDecoupledMolecules > 1:
                raise ValueError('Gromacs cannot handle more than one decoupled molecule.')
            protocol_dict["calc-lambda-neighbors"] = -1                             # Calculate MBAR energies.
            LambdaValues = self.protocol.getLambdaValues(type='dataframe')
            for name in ['fep', 'bonded', 'coul', 'vdw', 'restraint', 'mass',
                           'temperature']:
                if name in LambdaValues:
                    protocol_dict['{:<20}'.format("{}-lambdas".format(name))] = \
                        ' '.join(list(map('{:.5f}'.format, LambdaValues[name].to_list())))
            protocol_dict["init-lambda-state"] = self.protocol.getLambdaIndex()     # Current lambda value.
            protocol_dict["nstcalcenergy"] = 200                                    # Calculate energies every 200 steps.
            protocol_dict["nstdhdl"] = 200                                          # Write gradients every 200 steps.

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [f"{k} = {v}" for k, v in total_dict.items() if v is not None] + extra_lines

        return total_lines

    def generateSomdConfig(self, extra_options=None, extra_lines=None):
        """Outputs the current protocol in a format compatible with SOMD.

        Parameters
        ----------

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        Returns
        -------

        config : list
            The generated config list in a SOMD format.
        """

        extra_options = extra_options if extra_options is not None else {}
        extra_lines = extra_lines if extra_lines is not None else []

        # Define some miscellaneous defaults.
        protocol_dict = {"save coordinates": True}                                  # Save molecular coordinates.

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["minimise"] = True                                        # Minimisation simulation.
            protocol_dict["minimise maximum iterations"] = self._steps              # Maximum number of steps.
            protocol_dict["minimise tolerance"] = 1                                 # Convergence tolerance.
            protocol_dict["ncycles"] = 1                                            # Perform a single SOMD cycle.
            protocol_dict["nmoves"] = 1                                             # Perform a single MD move.
        else:
            # Get the report and restart intervals.
            report_interval = self._report_interval
            restart_interval = self._restart_interval

            # The restart and report intervals must be a multiple of the energy frequency,
            # which is 200 steps.
            if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                report_interval = int(200 * _math.ceil(report_interval / 200))
                restart_interval = int(200 * _math.ceil(restart_interval / 200))

            # The number of moves per cycle.
            nmoves = report_interval

            # The number of cycles, so that nmoves * ncycles is equal to self._steps.
            ncycles = max(1, self._steps // nmoves)

            # How many cycles need to pass before we write a trajectory frame.
            cycles_per_frame = max(1, restart_interval // nmoves)

            # How many time steps need to pass before we write a trajectory frame.
            buffer_freq = int(nmoves * ((restart_interval / nmoves) % 1))

            protocol_dict["ncycles"] = ncycles                                  # The number of SOMD cycles.
            protocol_dict["nmoves"] = nmoves                                    # The number of moves per cycle.
            protocol_dict["ncycles_per_snap"] = cycles_per_frame                # Cycles per trajectory write.
            protocol_dict["buffered coordinates frequency"] = buffer_freq       # Buffering frequency.
            timestep = self.protocol.getTimeStep().femtoseconds().value()
            protocol_dict["timestep"] = "%.2f femtosecond" % timestep           # Integration time step.

            # Use the Langevin Middle integrator if it is a 4 fs timestep
            if timestep >= 4.00:
                protocol_dict["integrator_type"] = "langevinmiddle"             # Langevin middle integrator
            else:
                pass

        # PBC.
        if self._has_water:
            protocol_dict["reaction field dielectric"] = "78.3"                 # Solvated box.
        if not self._has_box or not self._has_water:
            protocol_dict["cutoff type"] = "cutoffnonperiodic"                  # No periodic box.
        else:
            protocol_dict["cutoff type"] = "cutoffperiodic"                     # Periodic box.
        protocol_dict["cutoff distance"] = "8 angstrom"                        # Non-bonded cut-off.

        # Restraints.
        if isinstance(self.protocol, _Protocol.Equilibration) and self.protocol.getRestraint() is not None:
            raise _IncompatibleError("We currently don't support restraints with SOMD.")

        # Pressure control.
        protocol_dict["barostat"] = False                                       # Disable barostat (constant volume).
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    protocol_dict["barostat"] = True                            # Enable barostat.
                    pressure = self.protocol.getPressure().atm().value()
                    protocol_dict["pressure"] = "%.5f atm" % pressure           # Presure in atmosphere.
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if isinstance(self.protocol, _Protocol.Equilibration) and not self.protocol.isConstantTemp():
                raise _IncompatibleError("SOMD only supports constant temperature equilibration.")

            protocol_dict["thermostat"] = "True"                                # Turn on the thermostat.
            if not isinstance(self.protocol, _Protocol.Equilibration):
                protocol_dict["temperature"] = "%.2f kelvin" % self.protocol.getTemperature().kelvin().value()
            else:
                protocol_dict["temperature"] = "%.2f kelvin" % self.protocol.getStartTemperature().kelvin().value()

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            if not isinstance(self.protocol, _Protocol.Minimisation):
                protocol_dict["constraint"] = "hbonds-notperturbed"             # Handle hydrogen perturbations.
                protocol_dict["energy frequency"] = 200                         # Write gradients every 200 steps.

            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            protocol_dict["lambda array"] = ", ".join(protocol)
            protocol_dict["lambda_val"] = self.protocol.getLambda()             # Current lambda value.
            res_num = self.system.search("perturbable").residues()[0]._sire_object.number().value()
            protocol_dict["perturbed residue number"] = res_num                 # Perturbed residue number.

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [f"{k} = {v}" for k, v in total_dict.items() if v is not None] + extra_lines

        return total_lines

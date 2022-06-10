import itertools as _it
import math as _math
import warnings as _warnings

from ..Align._merge import _squash
from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Units.Time import nanosecond as _nanosecond

from .. import Protocol as _Protocol


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
            _warnings.warn(
                "No simulation box found. Assuming gas phase simulation.")
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
            steps = _math.ceil(self.protocol.getRunTime() /
                               self.protocol.getTimeStep())
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
        squashed_system = _squash(self.system)

        # When HMR is used, there can be no noshakemask
        HMR_on = hmr

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
            # Interval between reporting energies.
            "ntpr": self._report_interval,
            # Interval between saving restart files.
            "ntwr": self._restart_interval,
            # Trajectory sampling frequency.
            "ntwx": self._restart_interval,
            # Output coordinates as NetCDF.
            "ntxo": 2,
            "irest": int(self._restart),            # Whether to restart.
        }

        # Input.
        if self._restart:
            # Read coordinates and velocities.
            protocol_dict["ntx"] = 5
        else:
            # Only read coordinates from file.
            protocol_dict["ntx"] = 1

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
            # Set the minimisation method to XMIN
            protocol_dict["ntmin"] = 2
            protocol_dict["maxcyc"] = self._steps   # Set the number of steps.
            # Set the number of steepest descent steps.
            protocol_dict["ncyc"] = num_steep
        else:
            # Define the timestep
            timestep = self.protocol.getTimeStep().picoseconds().value()  # Get the time step in ps
            protocol_dict["dt"] = f"{timestep:.3f}"  # Time step.
            # Number of integration steps.
            protocol_dict["nstlim"] = self._steps

        # Constraints.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntc"] = 2                # Enable SHAKE.
            # Don't calculate forces for constrained bonds.
            protocol_dict["ntf"] = 2

        # PBC.
        if not self._has_box or not self._has_water:
            protocol_dict["ntb"] = 0                # No periodic box.
            protocol_dict["cut"] = "999."           # Non-bonded cut-off.
        else:
            protocol_dict["cut"] = "8.0"            # Non-bonded cut-off.
            protocol_dict["iwrap"] = 1              # Wrap the coordinates.

        # Restraints.
        if isinstance(self.protocol, _Protocol.Equilibration) or isinstance(self.protocol, _Protocol.Minimisation):
            
            # Get the restraint.
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
                    restraint_mask = self._amber_mask_from_indices(
                        [i + 1 for i in atom_idxs])

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
                                "AMBER atom 'restraintmask' exceeds 256 character limit!")

                    protocol_dict["ntr"] = 1
                    protocol_dict["restraint_wt"] = f"{self.protocol.getForceConstant().value():.1f}"
                    protocol_dict["restraintmask"] = f"\"{restraint_mask}\""

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    # Isotropic pressure scaling.
                    protocol_dict["ntp"] = 1
                    # Pressure in bar.
                    protocol_dict["pres0"] = f"{self.protocol.getPressure().bar().value():.5f}"
                    if isinstance(self.protocol, _Protocol.Equilibration):
                        # Monte Carlo barostat.
                        protocol_dict["barostat"] = 2
                    else:
                        # Monte Carlo barostat.
                        protocol_dict["barostat"] = 2
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntt"] = 3                # Langevin dynamics.
            protocol_dict["gamma_ln"] = 2           # Collision frequency (ps).
            if isinstance(self.protocol, _Protocol.Equilibration):
                temp0 = self.protocol.getStartTemperature().kelvin().value()
                temp1 = self.protocol.getEndTemperature().kelvin().value()
                if not self.protocol.isConstantTemp():
                    # Initial temperature.
                    protocol_dict["tempi"] = f"{temp0:.2f}"
                    # Final temperature.
                    protocol_dict["temp0"] = f"{temp1:.2f}"
                    protocol_dict["nmropt"] = 1
                    protocol_lines += [
                        f"&wt TYPE='TEMP0', istep1=0, istep2={self._steps}, value1={temp0:.2f}, value2={temp1:.2f} /"
                    ]
                else:
                    if not self._restart:
                        # Initial temperature.
                        protocol_dict["tempi"] = f"{temp0:.2f}"
                    # Constant temperature.
                    protocol_dict["temp0"] = f"{temp0:.2f}"
            else:
                temp = self.protocol.getTemperature().kelvin().value()
                if not self._restart:
                    # Initial temperature.
                    protocol_dict["tempi"] = f"{temp:.2f}"
                # Final temperature.
                protocol_dict["temp0"] = f"{temp:.2f}"

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            # Free energy mode.
            protocol_dict["icfe"] = 1
            # Use softcore potentials.
            protocol_dict["ifsc"] = 1
            # Remove SHAKE constraints.
            protocol_dict["ntf"] = 1
            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            # Number of lambda values.
            protocol_dict["mbar_states"] = len(protocol)
            # Lambda values.
            protocol_dict["mbar_lambda"] = ", ".join(protocol)
            # Current lambda value.
            protocol_dict["clambda"] = self.protocol.getLambda()
            if isinstance(self.protocol, _Protocol.Production):
                # Calculate MBAR energies.
                protocol_dict["ifmbar"] = 1
                # Output dVdl
                protocol_dict["logdvdl"] = 1
            # Atom masks based on if HMR.
            hmr = self.protocol.getHmr()
            protocol_dict = {**protocol_dict, **
                             self._generate_amber_fep_masks(hmr)}
            # Do not wrap the coordinates for equilibration steps 
            # as this can create issues for the free leg of FEP runs.
            if isinstance(self.protocol, _Protocol.Equilibration) or isinstance(self.protocol, _Protocol.Minimisation):
                protocol_dict["iwrap"] = 0

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        dict_lines = [self.protocol.__class__.__name__, "&cntrl"]
        dict_lines += [f"   {k}={v}," for k,
                       v in total_dict.items() if v is not None] + ["/"]
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
            # Interval between writing to the log file.
            "nstlog": self._report_interval,
            # Interval between writing to the energy file.
            "nstenergy": self._restart_interval,
            # Interval between writing to the trajectory file.
            "nstxout": self._restart_interval,
        }

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            # Minimisation simulation.
            protocol_dict["integrator"] = "steep"
        else:
            # Define the timestep in picoseconds
            timestep = self.protocol.getTimeStep().picoseconds().value()
            # Integration time step.
            protocol_dict["dt"] = f"{timestep:.3f}"
        # Number of integration steps.
        protocol_dict["nsteps"] = self._steps

        # Constraints.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if timestep >= 0.004:
                # If HMR, all constraints
                protocol_dict["constraints"] = "all-bonds"
            else:
                # Rigid bonded hydrogens.
                protocol_dict["constraints"] = "h-bonds"
            # Linear constraint solver.
            protocol_dict["constraint-algorithm"] = "LINCS"

        # PBC.
        # Simulate a fully periodic box.
        protocol_dict["pbc"] = "xyz"
        # Use Verlet pair lists.
        protocol_dict["cutoff-scheme"] = "Verlet"
        if self._has_box and self._has_water:
            # Use a grid to search for neighbours.
            protocol_dict["ns-type"] = "grid"
            # Rebuild neighbour list every 20 steps. Recommended in the manual for parallel simulations and/or non-bonded force calculation on the GPU.
            protocol_dict["nstlist"] = "20"
            # Set short-range cutoff.
            protocol_dict["rlist"] = "0.8"
            # Set van der Waals cutoff.
            protocol_dict["rvdw"] = "0.8"
            # Set Coulomb cutoff.
            protocol_dict["rcoulomb"] = "0.8"
            # Fast smooth Particle-Mesh Ewald.
            protocol_dict["coulombtype"] = "PME"
            # Dispersion corrections for energy and pressure.
            protocol_dict["DispCorr"] = "EnerPres"
        else:
            # Perform vacuum simulations by implementing pseudo-PBC conditions,
            # i.e. run calculation in a near-infinite box (333.3 nm).
            # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588
            # Single neighbour list (all particles interact).
            protocol_dict["nstlist"] = "1"
            # "Infinite" short-range cutoff.
            protocol_dict["rlist"] = "333.3"
            # "Infinite" van der Waals cutoff.
            protocol_dict["rvdw"] = "333.3"
            # "Infinite" Coulomb cutoff.
            protocol_dict["rcoulomb"] = "333.3"
            # Plain cut-off.
            protocol_dict["coulombtype"] = "Cut-off"
        # Twin-range van der Waals cut-off.
        protocol_dict["vdwtype"] = "Cut-off"

        # Restraints.
        # The actual restraints are added in _gromacs.py
        # TODO check for min if need refcoord
        if isinstance(self.protocol, _Protocol.Equilibration) or isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["refcoord-scaling"] = "all"

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    if isinstance(self.protocol, _Protocol.Equilibration):
                        # Barostat type.
                        protocol_dict["pcoupl"] = "c-rescale"
                    else:
                        # Barostat type.
                        protocol_dict["pcoupl"] = "parrinello-rahman"
                    # 1ps time constant for pressure coupling.
                    protocol_dict["tau-p"] = 1
                    # Pressure in bar.
                    protocol_dict["ref-p"] = f"{self.protocol.getPressure().bar().value():.5f}"
                    # Compressibility of water.
                    protocol_dict["compressibility"] = "4.5e-5"
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            # Langevin dynamics.
            protocol_dict["integrator"] = "sd"
            # A single temperature group for the entire system.
            protocol_dict["tc-grps"] = "system"
            # Collision frequency (ps).
            protocol_dict["tau-t"] = 2

            if not isinstance(self.protocol, _Protocol.Equilibration):
                protocol_dict["ref-t"] = "%.2f" % self.protocol.getTemperature().kelvin().value()
            elif self.protocol.isConstantTemp():
                protocol_dict["ref-t"] = "%.2f" % self.protocol.getStartTemperature().kelvin().value()

            # Heating/cooling protocol.
            elif not self.protocol.isConstantTemp():
                # still need a reference temperature for each group, even when heating/cooling
                protocol_dict["ref-t"] = "%.2f" % self.protocol.getEndTemperature().kelvin().value()
                # Work out the final time of the simulation.
                timestep = self.protocol.getTimeStep().picoseconds().value()
                end_time = _math.floor(timestep * self._steps)

                # Single sequence of annealing points.
                protocol_dict["annealing"] = "single"
                # Two annealing points for "system" temperature group.
                protocol_dict["annealing-npoints"] = 2

                # Linearly change temperature between start and end times.
                protocol_dict["annealing-time"] = "0 %d" % end_time
                protocol_dict["annealing-temp"] = "%.2f %.2f" % (
                    self.protocol.getStartTemperature().kelvin().value(),
                    self.protocol.getEndTemperature().kelvin().value(),
                )

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            # Free energy mode.
            protocol_dict["free-energy"] = "yes"
            # Calculate MBAR energies.
            protocol_dict["calc-lambda-neighbors"] = -1
            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            # Lambda values.
            protocol_dict["fep-lambdas"] = " ".join(protocol)
            lam = self.protocol.getLambda()
            idx = self.protocol.getLambdaValues().index(lam)
            # Current lambda value.
            protocol_dict["init-lambda-state"] = idx
            # Calculate energies every report interval steps.
            protocol_dict["nstcalcenergy"] = self._report_interval
            # Write gradients every report interval steps.
            protocol_dict["nstdhdl"] = self._report_interval

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [f"{k} = {v}" for k,
                       v in total_dict.items() if v is not None] + extra_lines

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
        # Save molecular coordinates.
        protocol_dict = {"save coordinates": True}

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            # Minimisation simulation.
            protocol_dict["minimise"] = True
            # Maximum number of steps.
            protocol_dict["minimise maximum iterations"] = self._steps
            # Convergence tolerance.
            protocol_dict["minimise tolerance"] = 1
            # Perform a single SOMD cycle.
            protocol_dict["ncycles"] = 1
            # Perform a single MD move.
            protocol_dict["nmoves"] = 1
        else:
            # Get the report and restart intervals.
            report_interval = self._report_interval
            restart_interval = self._restart_interval

            # The restart and report intervals must be a multiple of the energy frequency,
            # which is 200 steps.
            if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                report_interval = int(200 * _math.ceil(report_interval / 200))
                restart_interval = int(
                    200 * _math.ceil(restart_interval / 200))

            # The number of moves per cycle - want about 1 cycle per 1 ns.
            nmoves = int(
                max(1, ((self._steps) // ((self.protocol.getRunTime())/(1*_nanosecond)))))

            # The number of cycles, so that nmoves * ncycles is equal to self._steps.
            ncycles = int(max(1, self._steps // nmoves))

            # How many cycles need to pass before we write a trajectory frame.
            cycles_per_frame = max(1, restart_interval // nmoves)

            # How many time steps need to pass before we write a trajectory frame.
            buffer_freq = int(nmoves * ((restart_interval / nmoves) % 1))

            # The number of SOMD cycles.
            protocol_dict["ncycles"] = ncycles
            # The number of moves per cycle.
            protocol_dict["nmoves"] = nmoves
            # Cycles per trajectory write.
            protocol_dict["ncycles_per_snap"] = cycles_per_frame
            # Buffering frequency.
            protocol_dict["buffered coordinates frequency"] = buffer_freq
            timestep = self.protocol.getTimeStep().femtoseconds().value()
            # Integration time step.
            protocol_dict["timestep"] = "%.2f femtosecond" % timestep

            # Use the Langevin Middle integrator if it is a 4 fs timestep
            # Apply HMR if timestep is greater than 4 fs.
            if timestep >= 4.00:
                # Langevin middle integrator
                protocol_dict["integrator_type"] = "langevinmiddle"
                # Repartitioning factor done based on hmr protocol as extra option
            else:
                pass

        # PBC.
        if self._has_water:
            # Solvated box.
            protocol_dict["reaction field dielectric"] = "78.3"
        if not self._has_box or not self._has_water:
            # No periodic box.
            protocol_dict["cutoff type"] = "cutoffnonperiodic"
        else:
            # Periodic box.
            protocol_dict["cutoff type"] = "cutoffperiodic"
        # Non-bonded cut-off.
        protocol_dict["cutoff distance"] = "8 angstrom"

        # Restraints.
        if isinstance(self.protocol, _Protocol.Equilibration) or isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getRestraint() is not None:
                raise _IncompatibleError(
                    "We currently don't support restraints with SOMD.")

        # Pressure control.
        # Disable barostat (constant volume).
        protocol_dict["barostat"] = False
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    # Enable barostat.
                    protocol_dict["barostat"] = True
                    pressure = self.protocol.getPressure().atm().value()
                    # Presure in atmosphere.
                    protocol_dict["pressure"] = "%.5f atm" % pressure
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation")

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if isinstance(self.protocol, _Protocol.Equilibration) and not self.protocol.isConstantTemp():
                raise _IncompatibleError(
                    "SOMD only supports constant temperature equilibration.")

            # Turn on the thermostat.
            protocol_dict["thermostat"] = "True"

            if protocol_dict["integrator_type"] == "langevinmiddle":
                # Turn off the thermostat for langevin middle integrator.
                protocol_dict["thermostat"] = "False"
            else:
                pass

            if not isinstance(self.protocol, _Protocol.Equilibration):
                protocol_dict["temperature"] = "%.2f kelvin" % self.protocol.getTemperature(
                ).kelvin().value()
            else:
                protocol_dict["temperature"] = "%.2f kelvin" % self.protocol.getStartTemperature(
                ).kelvin().value()

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            if not isinstance(self.protocol, _Protocol.Minimisation):
                # Handle hydrogen perturbations.
                protocol_dict["constraint"] = "hbonds-notperturbed"
                # Write gradients every report interval steps.
                protocol_dict["energy frequency"] = self._report_interval

            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            protocol_dict["lambda array"] = ", ".join(protocol)
            # Current lambda value.
            protocol_dict["lambda_val"] = self.protocol.getLambda()
            # protocol_dict["minimise"] = True                                   # minimise at each window
            # Find the ligand, which will have the name LIG if created using BSS.
            lig_res_num = self.system.search(f"resname LIG")[
                0]._sire_object.number().value()
            protocol_dict["perturbed residue number"] = lig_res_num

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [f"{k} = {v}" for k,
                       v in total_dict.items() if v is not None] + extra_lines

        return total_lines

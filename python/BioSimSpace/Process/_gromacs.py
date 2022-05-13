######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for running simulations with GROMACS.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Gromacs"]

from .._Utils import _try_import

import math as _math
import os as _os
_pygtail = _try_import("pygtail")
import shutil as _shutil
import shlex as _shlex
import subprocess as _subprocess
import timeit as _timeit
import warnings as _warnings

from Sire import Base as _SireBase
from Sire import IO as _SireIO
from Sire import Maths as _SireMaths
from Sire import Units as _SireUnits
from Sire import Vol as _SireVol

from .. import _gmx_exe, _gmx_version
from .. import _isVerbose
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from ..Types._type import Type as _Type

from .. import IO as _IO
from .. import Protocol as _Protocol
from .. import Trajectory as _Trajectory
from .. import Types as _Types
from .. import Units as _Units
from .. import _Utils

from . import _process

from ._plumed import Plumed as _Plumed

class Gromacs(_process.Process):
    """A class for running simulations using GROMACS."""

    def __init__(self, system, protocol, exe=None, name="gromacs",
            work_dir=None, seed=None, property_map={}, ignore_warnings=False, show_errors=True):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol <BioSimSpace.Protocol>`
               The protocol for the GROMACS process.

           exe : str
               The full path to the GROMACS executable.

           name : str
               The name of the process.

           work_dir :
               The working directory for the process.

           seed : int
               A random number seed.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           ignore_warnings : bool
               Whether to ignore warnings when generating the binary run file
               with 'gmx grompp'. By default, these warnings are elevated to
               errors and will halt the program.

           show_errors : bool
               Whether to show warning/error messages when generating the binary
               run file.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed, property_map)

        # Set the package name.
        self._package_name = "GROMACS"

        # This process can generate trajectory data.
        self._has_trajectory = True

        if _gmx_exe is not None:
            self._exe = _gmx_exe
        else:
            if exe is not None:
                # Make sure executable exists.
                if _os.path.isfile(exe):
                    self._exe = exe
                else:
                    raise IOError("GROMACS executable doesn't exist: '%s'" % exe)
            else:
                raise _MissingSoftwareError("'BioSimSpace.Process.Gromacs' is not supported. "
                                            "Please install GROMACS (http://www.gromacs.org).")

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool.")
        self._show_errors = show_errors

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = _process._MultiDict()

        # Store the name of the GROMACS log file.
        self._log_file = "%s/%s.log" % (self._work_dir, name)

        # The names of the input files.
        self._gro_file = "%s/%s.gro" % (self._work_dir, name)
        self._top_file = "%s/%s.top" % (self._work_dir, name)

        # The name of the trajectory file.
        self._traj_file = "%s/%s.trr" % (self._work_dir, name)

        # Set the path for the GROMACS configuration file.
        self._config_file = "%s/%s.mdp" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._gro_file, self._top_file]

        # Initialise the PLUMED interface object.
        self._plumed = None

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # Create a copy of the system.
        system = self._system.copy()

        if isinstance(self._protocol, _Protocol.FreeEnergy):
            # Check that the system contains a perturbable molecule.
            if self._system.nPerturbableMolecules() == 0:
                raise ValueError("'BioSimSpace.Protocol.FreeEnergy' requires a "
                                 "perturbable molecule!")

            # Check that the perturbation type is supported..
            if self._protocol.getPerturbationType() != "full":
                msg = ("'BioSimSpace.Process.Gromacs' currently only supports the 'full' "
                       "perturbation type. Please use 'BioSimSpace.Process.Somd' "
                       "for multistep perturbation types.")
                raise NotImplementedError(msg)

        else:
            # Check for perturbable molecules and convert to the chosen end state.
            system = self._checkPerturbable(system)

        # Convert the water model topology so that it matches the GROMACS naming convention.
        system._set_water_topology("GROMACS")

        # GRO87 file.
        gro = _SireIO.Gro87(system._sire_object, self._property_map)
        gro.writeToFile(self._gro_file)

        # TOP file.
        top = _SireIO.GroTop(system._sire_object, self._property_map)
        top.writeToFile(self._top_file)

        # Create the binary input file name.
        self._tpr_file = "%s/%s.tpr" % (self._work_dir, self._name)
        self._input_files.append(self._tpr_file)

        # Generate the GROMACS configuration file.
        # Skip if the user has passed a custom config.
        if isinstance(self._protocol, _Protocol.Custom):
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Generate the dictionary of command-line arguments.
        self._generate_args()

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate GROMACS configuration file strings."""

        # Clear the existing configuration list.
        self._config = []

        # Check whether the system contains periodic box information.
        # For now, well not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if "space" in self._system._sire_object.propertyKeys():
            has_box = True
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # The list of configuration strings.
        # We don't repeatedly call addToConfig since this will run grommp
        # to re-compile the binary run input file each time.
        config = []

        # While the configuration parameters below share a lot of overlap,
        # we choose the keep them separate so that the user can modify options
        # for a given protocol in a single place.

        # Add configuration variables for a minimisation simulation.
        if isinstance(self._protocol, _Protocol.Minimisation):
            config.append("integrator = steep")             # Use steepest descent.
            config.append("nsteps = %d"
                % self._protocol.getSteps())                # Set the number of steps.
            config.append("nstxout = %d"
                % self._protocol.getSteps())                # Only write the final coordinates.
            if has_box and self._has_water:
                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("ns-type = grid")             # Use a grid to search for neighbours.
                config.append("rlist = 1.2")                # Set short-range cutoff.
                config.append("rvdw = 1.2")                 # Set van der Waals cutoff.
                config.append("rcoulomb = 1.2")             # Set Coulomb cutoff.
                config.append("coulombtype = PME")          # Fast smooth Particle-Mesh Ewald.
                config.append("DispCorr = EnerPres")        # Dispersion corrections for energy and pressure.
            else:
                # Perform vacuum simulations by implementing pseudo-PBC conditions,
                # i.e. run calculation in a near-infinite box (333.3 nm).
                # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588

                # Create a copy of the system.
                system = self._system.copy()

                # Create a 999.9 nm periodic box and apply to the system.
                space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
                system._sire_object.setProperty(self._property_map.get("space", "space"), space)

                # Re-write the GRO file.
                gro = _SireIO.Gro87(system._sire_object, self._property_map)
                gro.writeToFile(self._gro_file)

                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("nstlist = 1")                # Single neighbour list (all particles interact).
                config.append("rlist = 333.3")              # "Infinite" short-range cutoff.
                config.append("rvdw = 333.3")               # "Infinite" van der Waals cutoff.
                config.append("rcoulomb = 333.3")           # "Infinite" Coulomb cutoff.
                config.append("coulombtype = Cut-off")      # Plain cut-off.
            config.append("vdwtype = Cut-off")              # Twin-range van der Waals cut-off.

        # Add configuration variables for an equilibration simulation.
        elif isinstance(self._protocol, _Protocol.Equilibration):

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            config.append("integrator = sd")                    # Leap-frog stochastic dynamics.
            config.append("ld-seed = %d" % seed)                # Random number seed.
            config.append("dt = %.3f" % timestep)               # Integration time step.
            config.append("nsteps = %d" % steps)                # Number of integration steps.
            config.append("nstlog = %d" % report_interval)      # Interval between writing to the log file.
            config.append("nstenergy = %d" % report_interval)   # Interval between writing to the energy file.
            config.append("nstxout = %d" % restart_interval)    # Interval between writing to the trajectory file.
            if has_box and self._has_water:
                config.append("pbc = xyz")                      # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")         # Use Verlet pair lists.
                config.append("ns-type = grid")                 # Use a grid to search for neighbours.
                config.append("rlist = 1.2")                    # Set short-range cutoff.
                config.append("rvdw = 1.2")                     # Set van der Waals cutoff.
                config.append("rcoulomb = 1.2")                 # Set Coulomb cutoff.
                config.append("coulombtype = PME")              # Fast smooth Particle-Mesh Ewald.
                config.append("DispCorr = EnerPres")            # Dispersion corrections for energy and pressure.
            else:
                # Perform vacuum simulations by implementing pseudo-PBC conditions,
                # i.e. run calculation in a near-infinite box (333.3 nm).
                # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588

                # Create a copy of the system.
                system = self._system.copy()

                # Create a 999.9 nm periodic box and apply to the system.
                space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
                system._sire_object.setProperty(self._property_map.get("space", "space"), space)

                # Re-write the GRO file.
                gro = _SireIO.Gro87(system._sire_object, self._property_map)
                gro.writeToFile(self._gro_file)

                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("nstlist = 1")                # Single neighbour list (all particles interact).
                config.append("rlist = 333.3")              # "Infinite" short-range cutoff.
                config.append("rvdw = 333.3")               # "Infinite" van der Waals cutoff.
                config.append("rcoulomb = 333.3")           # "Infinite" Coulomb cutoff.
                config.append("coulombtype = Cut-off")      # Plain cut-off.
            config.append("vdwtype = Cut-off")              # Twin-range van der Waals cut-off.
            config.append("constraints = h-bonds")          # Rigid water molecules.
            config.append("constraint-algorithm = LINCS")   # Linear constraint solver.

            # Temperature control.
            # No need for "berendsen" with integrator "sd".
            config.append("tc-grps = system")               # A single temperature group for the entire system.
            config.append("tau-t = 2.0")                    # 2ps time constant for temperature coupling.
                                                            # Set the reference temperature.
            config.append("ref-t = %.2f" % self._protocol.getEndTemperature().kelvin().value())

            # Heating/cooling protocol.
            if not self._protocol.isConstantTemp():
                # Work out the final time of the simulation.
                end_time = _math.floor(timestep*steps)

                config.append("annealing = single")         # Single sequence of annealing points.
                config.append("annealing-npoints = 2")      # Two annealing points for "system" temperature group.

                # Linearly change temperature between start and end times.
                config.append("annealing-time = 0 %d" % end_time)
                config.append("annealing-temp = %.2f %.2f"
                    % (self._protocol.getStartTemperature().kelvin().value(),
                       self._protocol.getEndTemperature().kelvin().value()))

            # Pressure control.
            if self._protocol.getPressure() is not None and has_box and self._has_water:
                if _gmx_version >= 2021:
                    config.append("pcoupl = C-rescale")     # C-rescale barostat.
                else:
                    config.append("pcoupl = Berendsen")     # C-rescale barostat.
                config.append("tau-p = 1.0")                # 1ps time constant for pressure coupling.
                config.append("ref-p = %.5f"                # Pressure in bar.
                    % self._protocol.getPressure().bar().value())
                config.append("compressibility = 4.5e-5")   # Compressibility of water.

            # Add any position restraints.
            self._add_position_restraints(config)

        # Add configuration variables for a production simulation.
        elif isinstance(self._protocol, _Protocol.Production):

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            config.append("integrator = sd")                    # Leap-frog stochastic dynamics.
            config.append("ld-seed = %d" % seed)                # Random number seed.
            config.append("dt = %.3f" % timestep)               # Integration time step.
            config.append("nsteps = %d" % steps)                # Number of integration steps.
            config.append("init-step = %d"
                % self._protocol.getFirstStep())                # First time step.
            config.append("nstlog = %d" % report_interval)      # Interval between writing to the log file.
            config.append("nstenergy = %d" % report_interval)   # Interval between writing to the energy file.
            config.append("nstxout = %d" % restart_interval)    # Interval between writing to the trajectory file.
            if has_box and self._has_water:
                config.append("pbc = xyz")                      # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")         # Use Verlet pair lists.
                config.append("ns-type = grid")                 # Use a grid to search for neighbours.
                config.append("nstlist = 10")                   # Rebuild neighbour list every 10 steps.
                config.append("rlist = 1.2")                    # Set short-range cutoff.
                config.append("rvdw = 1.2")                     # Set van der Waals cutoff.
                config.append("rcoulomb = 1.2")                 # Set Coulomb cutoff.
                config.append("coulombtype = PME")              # Fast smooth Particle-Mesh Ewald.
                config.append("DispCorr = EnerPres")            # Dispersion corrections for energy and pressure.
            else:
                # Perform vacuum simulations by implementing pseudo-PBC conditions,
                # i.e. run calculation in a near-infinite box (333.3 nm).
                # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588

                # Create a copy of the system.
                system = self._system.copy()

                # Create a 999.9 nm periodic box and apply to the system.
                space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
                system._sire_object.setProperty(self._property_map.get("space", "space"), space)

                # Re-write the GRO file.
                gro = _SireIO.Gro87(system._sire_object, self._property_map)
                gro.writeToFile(self._gro_file)

                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("nstlist = 1")                # Single neighbour list (all particles interact).
                config.append("rlist = 333.3")              # "Infinite" short-range cutoff.
                config.append("rvdw = 333.3")               # "Infinite" van der Waals cutoff.
                config.append("rcoulomb = 333.3")           # "Infinite" Coulomb cutoff.
                config.append("coulombtype = Cut-off")      # Plain cut-off.
            config.append("vdwtype = Cut-off")              # Twin-range van der Waals cut-off.
            config.append("constraints = h-bonds")          # Rigid water molecules.
            config.append("constraint-algorithm = LINCS")   # Linear constraint solver.

            # Temperature control.
            # No need for "berendsen" with integrator "sd".
            config.append("tc-grps = system")               # A single temperature group for the entire system.
            config.append("tau-t = 2.0")                    # 2ps time constant for temperature coupling.
                                                            # Set the reference temperature.
            config.append("ref-t = %.2f" % self._protocol.getTemperature().kelvin().value())

            # Pressure control.
            if self._protocol.getPressure() is not None and has_box and self._has_water:
                if _gmx_version >= 2021:
                    config.append("pcoupl = C-rescale")     # C-rescale barostat.
                else:
                    config.append("pcoupl = Berendsen")     # C-rescale barostat.
                config.append("tau-p = 1.0")                # 1ps time constant for pressure coupling.
                config.append("ref-p = %.5f"                # Pressure in bar.
                    % self._protocol.getPressure().bar().value())
                config.append("compressibility = 4.5e-5")   # Compressibility of water.

        elif isinstance(self._protocol, _Protocol.FreeEnergy):

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            config.append("integrator = sd")                    # Leap-frog stochastic dynamics.
            config.append("ld-seed = %d" % seed)                # Random number seed.
            config.append("dt = %.3f" % timestep)               # Integration time step.
            config.append("nsteps = %d" % steps)                # Number of integration steps.
            config.append("nstlog = %d" % report_interval)      # Interval between writing to the log file.
            config.append("nstenergy = %d" % report_interval)   # Interval between writing to the energy file.
            config.append("nstxout = %d" % restart_interval)    # Interval between writing to the trajectory file.
            if has_box and self._has_water:
                config.append("pbc = xyz")                      # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")         # Use Verlet pair lists.
                config.append("ns-type = grid")                 # Use a grid to search for neighbours.
                config.append("nstlist = 10")                   # Rebuild neighbour list every 10 steps.
                config.append("rlist = 1.2")                    # Set short-range cutoff.
                config.append("rvdw = 1.2")                     # Set van der Waals cutoff.
                config.append("rcoulomb = 1.2")                 # Set Coulomb cutoff.
                config.append("coulombtype = PME")              # Fast smooth Particle-Mesh Ewald.
                config.append("DispCorr = EnerPres")            # Dispersion corrections for energy and pressure.
            else:
                # Perform vacuum simulations by implementing pseudo-PBC conditions,
                # i.e. run calculation in a near-infinite box (333.3 nm).
                # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588

                # Create a copy of the system.
                system = self._system.copy()

                # Create a 999.9 nm periodic box and apply to the system.
                space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
                system._sire_object.setProperty(self._property_map.get("space", "space"), space)

                # Re-write the GRO file.
                gro = _SireIO.Gro87(system._sire_object, self._property_map)
                gro.writeToFile(self._gro_file)

                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("nstlist = 1")                # Single neighbour list (all particles interact).
                config.append("rlist = 333.3")              # "Infinite" short-range cutoff.
                config.append("rvdw = 333.3")               # "Infinite" van der Waals cutoff.
                config.append("rcoulomb = 333.3")           # "Infinite" Coulomb cutoff.
                config.append("coulombtype = Cut-off")      # Plain cut-off.
            config.append("vdwtype = Cut-off")              # Twin-range van der Waals cut-off.
            config.append("constraints = h-bonds")          # Rigid water molecules.
            config.append("constraint-algorithm = LINCS")   # Linear constraint solver.

            # Temperature control.
            # No need for "berendsen" with integrator "sd".
            config.append("tc-grps = system")               # A single temperature group for the entire system.
            config.append("tau-t = 2.0")                    # 2ps time constant for temperature coupling.
                                                            # Set the reference temperature.
            config.append("ref-t = %.2f" % self._protocol.getTemperature().kelvin().value())

            # Pressure control.
            if self._protocol.getPressure() is not None and has_box and self._has_water:
                if _gmx_version >= 2021:
                    config.append("pcoupl = C-rescale")     # C-rescale barostat.
                else:
                    config.append("pcoupl = Berendsen")     # C-rescale barostat.
                config.append("tau-p = 1.0")                # 1ps time constant for pressure coupling.
                config.append("ref-p = %.5f"                # Pressure in bar.
                    % self._protocol.getPressure().bar().value())
                config.append("compressibility = 4.5e-5")   # Compressibility of water.

            # Extract the lambda array.
            lam_vals = self._protocol.getLambdaValues()

            # Determine the index of the lambda value.
            idx = self._protocol.getLambdaIndex()

            # Free energy parameters.
            config.append("free-energy = yes")              # Free energy simulation.
            config.append("init-lambda-state = %d" % idx)   # Index of the lambda value.
            config.append("fep-lambdas = %s" \
                % " ".join([str(x) for x in lam_vals]))
            config.append("couple-lambda0 = vdw-q")         # All interactions on at lambda = 0
            config.append("couple-lambda1 = vdw-q")         # All interactions on at lambda = 1
            config.append("calc-lambda-neighbors = -1")     # Write all lambda values.
            config.append("nstcalcenergy = 250")            # Calculate energies every 250 steps.
            config.append("nstdhdl = 250")                  # Write gradients every 250 steps.

        # Add configuration variables for a metadynamics simulation.
        elif isinstance(self._protocol, _Protocol.Metadynamics):

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            config.append("integrator = sd")                    # Leap-frog stochastic dynamics.
            config.append("ld-seed = %d" % seed)                # Random number seed.
            config.append("dt = %.3f" % timestep)               # Integration time step.
            config.append("nsteps = %d" % steps)                # Number of integration steps.
            config.append("nstlog = %d" % report_interval)      # Interval between writing to the log file.
            config.append("nstenergy = %d" % report_interval)   # Interval between writing to the energy file.
            config.append("nstxout = %d" % restart_interval)    # Interval between writing to the trajectory file.
            if has_box and self._has_water:
                config.append("pbc = xyz")                      # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")         # Use Verlet pair lists.
                config.append("ns-type = grid")                 # Use a grid to search for neighbours.
                config.append("nstlist = 10")                   # Rebuild neighbour list every 10 steps.
                config.append("rlist = 1.2")                    # Set short-range cutoff.
                config.append("rvdw = 1.2")                     # Set van der Waals cutoff.
                config.append("rcoulomb = 1.2")                 # Set Coulomb cutoff.
                config.append("coulombtype = PME")              # Fast smooth Particle-Mesh Ewald.
                config.append("DispCorr = EnerPres")            # Dispersion corrections for energy and pressure.
            else:
                # Perform vacuum simulations by implementing pseudo-PBC conditions,
                # i.e. run calculation in a near-infinite box (333.3 nm).
                # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588

                # Create a copy of the system.
                system = self._system.copy()

                # Create a 999.9 nm periodic box and apply to the system.
                space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
                system._sire_object.setProperty(self._property_map.get("space", "space"), space)

                # Re-write the GRO file.
                gro = _SireIO.Gro87(system._sire_object, self._property_map)
                gro.writeToFile(self._gro_file)

                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("nstlist = 1")                # Single neighbour list (all particles interact).
                config.append("rlist = 333.3")              # "Infinite" short-range cutoff.
                config.append("rvdw = 333.3")               # "Infinite" van der Waals cutoff.
                config.append("rcoulomb = 333.3")           # "Infinite" Coulomb cutoff.
                config.append("coulombtype = Cut-off")      # Plain cut-off.
            config.append("vdwtype = Cut-off")              # Twin-range van der Waals cut-off.
            config.append("constraints = h-bonds")          # Rigid water molecules.
            config.append("constraint-algorithm = LINCS")   # Linear constraint solver.

            # Temperature control.
            # No need for "berendsen" with integrator "sd".
            config.append("tc-grps = system")               # A single temperature group for the entire system.
            config.append("tau-t = 2.0")                    # 2ps time constant for temperature coupling.
                                                            # Set the reference temperature.
            config.append("ref-t = %.2f" % self._protocol.getTemperature().kelvin().value())

            # Pressure control.
            if self._protocol.getPressure() is not None and has_box and self._has_water:
                if _gmx_version >= 2021:
                    config.append("pcoupl = C-rescale")     # C-rescale barostat.
                else:
                    config.append("pcoupl = Berendsen")     # C-rescale barostat.
                config.append("tau-p = 1.0")                # 1ps time constant for pressure coupling.
                config.append("ref-p = %.5f"                # Pressure in bar.
                    % self._protocol.getPressure().bar().value())
                config.append("compressibility = 4.5e-5")   # Compressibility of water.

            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(self._work_dir)
            plumed_config, auxiliary_files = self._plumed.createConfig(self._system,
                                                                       self._protocol,
                                                                       self._property_map)
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(file, self._work_dir + f"/{file_name}")
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getFreeEnergy", self._getFreeEnergy)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "sampleConfigurations", self._sampleConfigurations)
            setattr(self, "getTime", self._getTime)

        # Add configuration variables for a steered molecular dynamics protocol.
        elif isinstance(self._protocol, _Protocol.Steering):

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            config.append("integrator = sd")                    # Leap-frog stochastic dynamics.
            config.append("ld-seed = %d" % seed)                # Random number seed.
            config.append("dt = %.3f" % timestep)               # Integration time step.
            config.append("nsteps = %d" % steps)                # Number of integration steps.
            config.append("nstlog = %d" % report_interval)      # Interval between writing to the log file.
            config.append("nstenergy = %d" % report_interval)   # Interval between writing to the energy file.
            config.append("nstxout = %d" % restart_interval)    # Interval between writing to the trajectory file.
            if has_box and self._has_water:
                config.append("pbc = xyz")                      # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")         # Use Verlet pair lists.
                config.append("ns-type = grid")                 # Use a grid to search for neighbours.
                config.append("nstlist = 10")                   # Rebuild neighbour list every 10 steps.
                config.append("rlist = 1.2")                    # Set short-range cutoff.
                config.append("rvdw = 1.2")                     # Set van der Waals cutoff.
                config.append("rcoulomb = 1.2")                 # Set Coulomb cutoff.
                config.append("coulombtype = PME")              # Fast smooth Particle-Mesh Ewald.
                config.append("DispCorr = EnerPres")            # Dispersion corrections for energy and pressure.
            else:
                # Perform vacuum simulations by implementing pseudo-PBC conditions,
                # i.e. run calculation in a near-infinite box (333.3 nm).
                # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588

                # Create a copy of the system.
                system = self._system.copy()

                # Create a 999.9 nm periodic box and apply to the system.
                space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
                system._sire_object.setProperty(self._property_map.get("space", "space"), space)

                # Re-write the GRO file.
                gro = _SireIO.Gro87(system._sire_object, self._property_map)
                gro.writeToFile(self._gro_file)

                config.append("pbc = xyz")                  # Simulate a fully periodic box.
                config.append("cutoff-scheme = Verlet")     # Use Verlet pair lists.
                config.append("nstlist = 1")                # Single neighbour list (all particles interact).
                config.append("rlist = 333.3")              # "Infinite" short-range cutoff.
                config.append("rvdw = 333.3")               # "Infinite" van der Waals cutoff.
                config.append("rcoulomb = 333.3")           # "Infinite" Coulomb cutoff.
                config.append("coulombtype = Cut-off")      # Plain cut-off.
            config.append("vdwtype = Cut-off")              # Twin-range van der Waals cut-off.
            config.append("constraints = h-bonds")          # Rigid water molecules.
            config.append("constraint-algorithm = LINCS")   # Linear constraint solver.

            # Temperature control.
            # No need for "berendsen" with integrator "sd".
            config.append("tc-grps = system")               # A single temperature group for the entire system.
            config.append("tau-t = 2.0")                    # 2ps time constant for temperature coupling.
                                                            # Set the reference temperature.
            config.append("ref-t = %.2f" % self._protocol.getTemperature().kelvin().value())

            # Pressure control.
            if self._protocol.getPressure() is not None and has_box and self._has_water:
                if _gmx_version >= 2021:
                    config.append("pcoupl = C-rescale")     # C-rescale barostat.
                else:
                    config.append("pcoupl = Berendsen")     # C-rescale barostat.
                config.append("tau-p = 1.0")                # 1ps time constant for pressure coupling.
                config.append("ref-p = %.5f"                # Pressure in bar.
                    % self._protocol.getPressure().bar().value())
                config.append("compressibility = 4.5e-5")   # Compressibility of water.

            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(self._work_dir)
            plumed_config, auxiliary_files = self._plumed.createConfig(self._system,
                                                                       self._protocol,
                                                                       self._property_map)
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(file, self._work_dir + f"/{file_name}")
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "getTime", self._getTime)

        # Set the configuration.
        self.setConfig(config)

        # Flag that this isn't a custom protocol.
        self._protocol._setCustomised(False)

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("mdrun", True)          # Use mdrun.
        self.setArg("-v", True)             # Verbose output.
        self.setArg("-deffnm", self._name)  # Output file prefix.

        # Metadynamics and steered MD arguments.
        if isinstance(self._protocol, (_Protocol.Metadynamics, _Protocol.Steering)):
            self.setArg("-plumed", "plumed.dat")

    def _generate_binary_run_file(self):
        """Use grommp to generate the binary run input file."""

        # Create the name of the output mdp file.
        mdp_out = _os.path.dirname(self._config_file) + \
                  "/%s.out.mdp" % _os.path.basename(self._config_file).split(".")[0]

        # Use grompp to generate the portable binary run input file.
        command = "%s grompp -f %s -po %s -c %s -p %s -r %s -o %s" \
            % (self._exe, self._config_file, mdp_out, self._gro_file,
               self._top_file, self._gro_file, self._tpr_file)

        # Warnings don't trigger an error. Set to a suitably large number.
        if self._ignore_warnings:
            command += " --maxwarn 1000"

        # Run the command.
        proc = _subprocess.run(_shlex.split(command), shell=False, text=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Check that grompp ran successfully.
        if proc.returncode != 0:
            # Handle errors and warnings.
            if self._show_errors:
                # Capture errors and warnings from the grompp output.
                errors = []
                warnings = []
                is_error = False
                is_warn = False
                lines = proc.stderr.split("\n")
                for line in lines:
                    line = line.strip()
                    if line[0:5] == "ERROR" or is_error:
                        if line == "":
                            is_error = False
                            continue
                        errors.append(line)
                        is_error = True
                    elif line[0:7] == "WARNING" or is_warn:
                        if line == "":
                            is_warn = False
                            continue
                        warnings.append(line)
                        is_warn = True

                error_string = "\n  ".join(errors)
                warning_string = "\n  ".join(warnings)

                exception_string = "Unable to generate GROMACS binary run input file.\n"
                if len(errors) > 0:
                    exception_string += "\n'gmx grompp' reported the following errors:\n"   \
                                     + f"{error_string}\n"
                if len(warnings) > 0:
                    exception_string += "\n'gmx grompp' reported the following warnings:\n" \
                                     + f"{warning_string}\n"                                \
                                     +  "\nUse 'ignore_warnings' to ignore warnings."

                raise RuntimeError(exception_string)

            else:
                raise RuntimeError("Unable to generate GROMACS binary run input file. "
                                   "Use 'show_errors=True' to display errors/warnings.")

    def addToConfig(self, config):
        """Add a string to the configuration list.

           Parameters
           ----------

           config : str, [str]
               A configuration string, a list of configuration strings, or a
               path to a configuration file.
        """

        # Call the base class method.
        super().addToConfig(config)

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file()

    def resetConfig(self):
        """Reset the configuration parameters."""
        self._generate_config()

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file()

    def setConfig(self, config):
        """Set the list of configuration file strings.

           Parameters
           ----------

           config : str, [str]
               The list of configuration strings, or a path to a configuration
               file.
        """

        # Call the base class method.
        super().setConfig(config)

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file()

    def start(self):
        """Start the GROMACS process.

           Returns
           -------

           process : :class:`Process.Gromacs <BioSimSpace.Process.Gromacs>`
               A handle to the GROMACS process.
        """

        # The process is currently queued.
        if self.isQueued():
            return

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Clear any existing output.
        self._clear_output()

        # Run the process in the working directory.
        with _Utils.cd(self._work_dir):

            # Create the arguments string list.
            args = self.getArgStringList()

            # Write the command-line process to a README.txt file.
            with open("README.txt", "w") as f:

                # Set the command-line string.
                self._command = "%s " % self._exe + self.getArgString()

                # Write the command to file.
                f.write("# GROMACS was run with the following command:\n")
                f.write("%s\n" % self._command)

            # Start the timer.
            self._timer = _timeit.default_timer()

            # Start the simulation.
            self._process = _SireBase.Process.run(self._exe, args,
                "%s.out" % self._name, "%s.out" % self._name)

            # For historical reasons (console message aggregation with MPI), Gromacs
            # writes the majority of its output to stderr. For user convenience, we
            # redirect all output to stdout, and place a message in the stderr file
            # to highlight this.
            with open(self._stderr_file, "w") as f:
                f.write("All output has been redirected to the stdout stream!\n")

        return self

    def getSystem(self, block="AUTO"):
        """Get the latest molecular system.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The latest molecular system.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        # Minimisation trajectories have a single frame, i.e. the final state.
        if isinstance(self._protocol, _Protocol.Minimisation):
            time = 0*_Units.Time.nanosecond
        # Get the current simulation time.
        else:
            time = self.getTime()

        # Grab the most recent frame from the trajectory file.
        return self._getFrame(time)

    def getCurrentSystem(self):
        """Get the latest molecular system.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The latest molecular system.
        """
        return self.getSystem(block=False)

    def getTrajectory(self, block="AUTO"):
        """Return a trajectory object.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           trajectory : :class:`System <BioSimSpace.Trajectory.Trajectory>`
               The latest trajectory object.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        try:
            # Locate the trajectory file.
            traj_file = self._find_trajectory_file()

            if traj_file is None:
                return None
            else:
                self._traj_file = traj_file

            return _Trajectory.Trajectory(process=self)

        except:
            return None

    def getFrame(self, index):
        """Return a specific trajectory frame.

           Parameters
           ----------

           index : int
               The index of the frame.

           Returns
           -------

           frame : :class:`System <BioSimSpace._SireWrappers.System>`
               The System object of the corresponding frame.
        """

        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        max_index = int((self._protocol.getRunTime() / self._protocol.getTimeStep())
                  / self._protocol.getRestartInterval()) - 1

        if index < 0 or index > max_index:
            raise ValueError(f"'index' must be in range [0, {max_index}].")

        try:
            time = index * self._protocol.getRestartInterval() \
                 * self._protocol.getTimeStep()

            with _warnings.catch_warnings():
                system = self._getFrame(time)

            return self._getFrame(time)

        except:
            return None

    def getRecord(self, record, time_series=False, unit=None, block="AUTO"):
        """Get a record from the stdout dictionary.

           Parameters
           ----------

           record : str
               The record key.

           time_series : bool
               Whether to return a list of time series records.

           unit : :class:`Unit <BioSimSpace.Units>`
               The unit to convert the record to.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           record : :class:`Type <BioSimSpace.Types>`
               The matching record.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_stdout_dict()
        return self._get_stdout_record(record, time_series, unit)

    def getCurrentRecord(self, record, time_series=False, unit=None):
        """Get a current record from the stdout dictionary.

           Parameters
           ----------

           record : str
               The record key.

           time_series : bool
               Whether to return a list of time series records.

           unit : :class:`Unit <BioSimSpace.Units>`
               The unit to convert the record to.

           Returns
           -------

           record : :class:`Type <BioSimSpace.Types>`
               The matching record.
        """
        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_stdout_dict()
        return self._get_stdout_record(record, time_series, unit)

    def getRecords(self, block="AUTO"):
        """Return the dictionary of stdout time-series records.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
              The dictionary of time-series records.
        """
        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        return self._stdout_dict.copy()

    def getCurrentRecords(self):
        """Return the current dictionary of stdout time-series records.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
              The dictionary of time-series records.
        """
        return self.getRecords(block=False)

    def getTime(self, time_series=False, block="AUTO"):
        """Get the simulation time.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           time : :class:`Time <BioSimSpace.Types.Time>`
               The current simulation time in nanoseconds.
        """

        if isinstance(self._protocol, _Protocol.Minimisation):
            return None

        else:
            return self.getRecord("TIME", time_series, _Units.Time.picosecond, block)

    def getCurrentTime(self, time_series=False):
        """Get the current simulation time.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           time : :class:`Time <BioSimSpace.Types.Time>`
               The current simulation time in nanoseconds.
        """
        return self.getTime(time_series, block=False)

    def getStep(self, time_series=False, block="AUTO"):
        """Get the number of integration steps.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           step : int
               The current number of integration steps.
        """
        return self.getRecord("STEP", time_series, None, block)

    def getCurrentStep(self, time_series=False):
        """Get the current number of integration steps.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           step : int
               The current number of integration steps.
        """
        return self.getStep(time_series, block=False)

    def getBondEnergy(self, time_series=False, block="AUTO"):
        """Get the bond energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The bond energy.
        """
        return self.getRecord("BOND", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentBondEnergy(self, time_series=False):
        """Get the current bond energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The bond energy.
        """
        return self.getBondEnergy(time_series, block=False)

    def getAngleEnergy(self, time_series=False, block="AUTO"):
        """Get the angle energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The angle energy.
        """
        return self.getRecord("ANGLE", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentAngleEnergy(self, time_series=False):
        """Get the current angle energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The angle energy.
        """
        return self.getAngleEnergy(time_series, block=False)

    def getDihedralEnergy(self, time_series=False, block="AUTO"):
        """Get the dihedral energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The dihedral energy.
        """
        return self.getRecord("PROPERDIH", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentDihedralEnergy(self, time_series=False):
        """Get the current dihedral energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The dihedral energy.
        """
        return self.getDihedralEnergy(time_series, block=False)

    def getImproperEnergy(self, time_series=False, block="AUTO"):
        """Get the improper energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The improper energy.
        """
        return self.getRecord("IMPRPROPERDIH", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentImproperEnergy(self, time_series=False):
        """Get the current improper energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The improper energy.
        """
        return self.getImproperEnergy(time_series, block=False)

    def getLennardJones14(self, time_series=False, block="AUTO"):
        """Get the Lennard-Jones energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Lennard-Jones energy.
        """
        return self.getRecord("LJ14", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentLennardJones14(self, time_series=False):
        """Get the current Lennard-Jones energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Lennard-Jones energy.
        """
        return self.getLennardJones14(time_series, block=False)

    def getLennardJonesSR(self, time_series=False, block="AUTO"):
        """Get the short-range Lennard-Jones energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The short-range Lennard-Jones energy.
        """
        return self.getRecord("LJSR", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentLennardJonesSR(self, time_series=False):
        """Get the current short-range Lennard-Jones energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Lennard-Jones energy.
        """
        return self.getLennardJonesSR(time_series, block=False)

    def getCoulomb14(self, time_series=False, block="AUTO"):
        """Get the Coulomb energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Coulomb energy.
        """
        return self.getRecord("COULOMB14", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentCoulomb14(self, time_series=False):
        """Get the current Coulomb energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Coulomb energy.
        """
        return self.getCoulomb14(time_series, block=False)

    def getCoulombSR(self, time_series=False, block="AUTO"):
        """Get the short-range Coulomb energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Coulomb energy.
        """
        return self.getRecord("COULOMBSR", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentCoulombSR(self, time_series=False):
        """Get the current short-range Coulomb energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Coulomb energy.
        """
        return self.getCoulombSR(time_series, block=False)

    def getCoulombReciprocal(self, time_series=False, block="AUTO"):
        """Get the reciprocal space Coulomb energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Coulomb energy.
        """
        return self.getRecord("COULRECIP", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentCoulombReciprocal(self, time_series=False):
        """Get the current reciprocal space Coulomb energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The Coulomb energy.
        """
        return self.getCoulombReciprocal(time_series, block=False)

    def getDispersionCorrection(self, time_series=False, block="AUTO"):
        """Get the dispersion correction.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The dispersion correction.
        """
        return self.getRecord("DISPERCORR", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentDispersionCorrection(self, time_series=False):
        """Get the current dispersion correction.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The dispersion correction.
        """
        return self.getDispersionCorrection(time_series, block=False)

    def getRestraintEnergy(self, time_series=False, block="AUTO"):
        """Get the position restraint energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The dispersion correction.
        """
        return self.getRecord("POSITIONREST", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentRestraintEnergy(self, time_series=False):
        """Get the current position restraint energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The dispersion correction.
        """
        return self.getRestraintEnergy(time_series, block=False)

    def getPotentialEnergy(self, time_series=False, block="AUTO"):
        """Get the potential energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The potential energy.
        """
        return self.getRecord("POTENTIAL", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentPotentialEnergy(self, time_series=False):
        """Get the current potential energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The potential energy.
        """
        return self.getPotentialEnergy(time_series, block=False)

    def getKineticEnergy(self, time_series=False, block="AUTO"):
        """Get the kinetic energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The kinetic energy.
        """
        return self.getRecord("KINETICEN", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentKineticEnergy(self, time_series=False):
        """Get the current kinetic energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The kinetic energy.
        """
        return self.getKineticEnergy(time_series, block=False)

    def getTotalEnergy(self, time_series=False, block="AUTO"):
        """Get the total energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The total energy.
        """
        return self.getRecord("TOTALENERGY", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentTotalEnergy(self, time_series=False):
        """Get the current total energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The total energy.
        """
        return self.getTotalEnergy(time_series, block=False)

    def getConservedEnergy(self, time_series=False, block="AUTO"):
        """Get the conserved energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The conserved energy.
        """
        return self.getRecord("CONSERVEDEN", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentConservedEnergy(self, time_series=False):
        """Get the current conserved energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The conserved energy.
        """
        return self.getConservedEnergy(time_series, block=False)

    def getTemperature(self, time_series=False, block="AUTO"):
        """Get the temperature.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The temperature.
        """
        return self.getRecord("TEMPERATURE", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperature(self, time_series=False):
        """Get the current temperature.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The current temperature.
        """
        return self.getTemperature(time_series, block=False)

    def getPressure(self, time_series=False, block="AUTO"):
        """Get the pressure.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure.
        """
        return self.getRecord("PRESSURE", time_series, _Units.Pressure.bar, block)

    def getCurrentPressure(self, time_series=False):
        """Get the current pressure.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The current pressure.
        """
        return self.getPressure(time_series, block=False)

    def getPressureDC(self, time_series=False, block="AUTO"):
        """Get the DC pressure.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The DC pressure.
        """
        return self.getRecord("PRESDC", time_series, _Units.Pressure.bar, block)

    def getCurrentPressureDC(self, time_series=False):
        """Get the current DC pressure.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The current pressure.
        """
        return self.getPressureDC(time_series, block=False)

    def getConstraintRMSD(self, time_series=False, block="AUTO"):
        """Get the RMSD of the constrained atoms.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           length : :class:`Length <BioSimSpace.Types.Length>`
               The constrained RMSD.
        """
        return self.getRecord("CONSTRRMSD", time_series, _Units.Length.nanometer, block)

    def getCurrentConstraintRMSD(self, time_series=False):
        """Get the current RMSD of the constrained atoms.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           length : :class:`Length <BioSimSpace.Types.Length>`
               The current constrained RMSD.
        """
        return self.getConstraintRMSD(time_series, block=False)

    def stdout(self, n=10):
        """Print the last n lines of the stdout buffer.

           Parameters
           ----------

           n : int
               The number of lines to print.
        """

        # Note that thermodynamic records, e.g. energy, pressure, temperature,
        # are redirected to a log file.

        # Ensure that the number of lines is positive.
        if n < 0:
            raise ValueError("The number of lines must be positive!")

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())

        # Get the current number of lines.
        num_lines = len(self._stdout)

        # Set the line from which to start printing.
        if num_lines < n:
            start = 0
        else:
            start = num_lines - n

        # Print the lines.
        for x in range(start, num_lines):
            print(self._stdout[x])

    def _add_position_restraints(self, config):
        """Helper function to add position restraints.

           Parameters
           ----------

           config : [str]
               The list of configuration strings.
        """

        # Get the restraint type.
        restraint = self._protocol.getRestraint()

        if restraint is not None:

            # Get the force constant in units of kJ_per_mol/nanometer**2
            force_constant = self._protocol.getForceConstant()._sire_unit
            force_constant = force_constant.to(_SireUnits.kJ_per_mol/_SireUnits.nanometer2)

            # Scale reference coordinates with the scaling matrix of the pressure coupling.
            config.append("refcoord-scaling = all")

            # Copy the user property map.
            property_map = self._property_map.copy()

            # Parse the topology in serial to ensure that molecules are
            # ordered correctly. Don't sort based on name.
            property_map["parallel"] = _SireBase.wrap(False)
            property_map["sort"] = _SireBase.wrap(False)

            # Create a copy of the system.
            system = self._system.copy()

            # Convert to the lambda = 0 state if this is a perturbable system.
            system = self._checkPerturbable(system)

            # Convert the water model topology so that it matches the GROMACS naming convention.
            system._set_water_topology("GROMACS")

            # Create a GROMACS topology object.
            top = _SireIO.GroTop(system._sire_object, property_map)

            # Get the top file as a list of lines.
            top_lines = top.lines()

            # List of 'moleculetype' record indices.
            moltypes_top_idx = []

            # Store the line index for the start of each 'moleculetype' record.
            for idx, line in enumerate(top_lines):
                if "[ moleculetype ]" in line or "[ system ]" in line:
                    moltypes_top_idx.append(idx)

            # Create a dictionary to store the indices of the molecules in the
            # system for each GROMACS molecule type and the inverse.
            moltypes_sys_idx = {}
            sys_idx_moltypes = {}

            # Convert the topology to a GROMACS system.
            gro_system = top.groSystem()

            # Initialise the dictionary for each type.
            for mol_type in gro_system.uniqueTypes():
                moltypes_sys_idx[mol_type] = []

            # Now loop over each molecule and store the indices of the molecules
            # in the system that are of each type as well as a mapping from the
            # molecule index to the GROMACS molecule type.
            for idx, mol_type in enumerate(gro_system):
                moltypes_sys_idx[mol_type].append(idx)
                sys_idx_moltypes[idx] = mol_type

            # A keyword restraint.
            if isinstance(restraint, str):

                # The number of restraint files.
                num_restraint = 1

                # Loop over all of the molecule types and create a position
                # restraint file for each.
                for mol_type_idx, (mol_type, mol_idxs) in enumerate(moltypes_sys_idx.items()):

                    # Initialise a list of restrained atom indices.
                    restrained_atoms = []

                    # Loop over each molecule in the system that matches this
                    # type and append any atoms matching the restraint.
                    for idx, mol_idx in enumerate(mol_idxs):

                        # Get the indices of any restrained atoms in this molecule,
                        # making sure that indices are relative to the molecule.
                        atom_idxs = self._system.getRestraintAtoms(restraint,
                                                                   mol_index=mol_idx,
                                                                   is_absolute=False,
                                                                   allow_zero_matches=True)

                        # Store the atom index if it hasn't already been recorded.
                        for atom_idx in atom_idxs:
                            if not atom_idx in restrained_atoms:
                                restrained_atoms.append(atom_idx)

                    # Write the position restraint file for this molecule.
                    if len(restrained_atoms) > 0:
                        # Create the file name.
                        restraint_file = "%s/posre_%04d.itp" % (self._work_dir, num_restraint)

                        with open(restraint_file, "w") as file:
                            # Write the header.
                            file.write("[ position_restraints ]\n")
                            file.write(";  i funct       fcx        fcy        fcz\n")

                            # Write restraints for each atom.
                            for atom_idx in restrained_atoms:
                                file.write(f"{atom_idx+1:4}    1       {force_constant}       {force_constant}       {force_constant}\n")

                        # Include the position restraint file in the correct place within
                        # the topology file. We put the additional include directive at the
                        # end of the block so we move to the line before the next moleculetype
                        # record.
                        new_top_lines = top_lines[:moltypes_top_idx[mol_type_idx+1]-1]

                        # Append the additional information.
                        new_top_lines.append('#include "%s"' % restraint_file)
                        new_top_lines.append("")

                        # Now extend with the remainder of the file.
                        new_top_lines.extend(top_lines[moltypes_top_idx[mol_type_idx+1]:])

                        # Overwrite the topology file lines.
                        top_lines = new_top_lines

                        # Increment the number of restraint files.
                        num_restraint += 1

                        # Append the restraint file to the list of autogenerated inputs.
                        self._input_files.append(restraint_file)

                # Write the updated topology to file.
                with open(self._top_file, "w") as file:
                    for line in top_lines:
                        file.write("%s\n" % line)

            # A user-defined list of atoms indices.
            else:

                # Create an empty multi-dict for each molecule type.
                mol_atoms = {}
                for mol_type in gro_system.uniqueTypes():
                    mol_atoms[mol_type] = []

                # Now work out which MolNum corresponds to each atom in the restraint.
                for idx in restraint:
                    try:
                        # Get the molecule index and relative atom index.
                        mol_idx, atom_idx = self._system._getRelativeIndices(idx)

                        # Get the type associated with this molecule index.
                        mol_type = sys_idx_moltypes[mol_idx]

                        # Append this atom if it's not already been recorded.
                        if not atom_idx in mol_atoms[mol_type]:
                            mol_atoms[mol_type].append(atom_idx)

                    except Exception as e:
                        msg = "Unable to find restrained atom in the system?"
                        if _isVerbose():
                            raise ValueError(msg) from e
                        else:
                            raise ValueError(msg) from None

                # The number of restraint files.
                num_restraint = 1

                # Loop over all of the molecule types and create a position
                # restraint file for each.
                for mol_type_idx, (mol_type, atom_idxs) in enumerate(mol_atoms.items()):

                    # Write the position restraint file for this molecule.
                    if len(atom_idxs) > 0:
                        # Create the file name.
                        restraint_file = "%s/posre_%04d.itp" % (self._work_dir, num_restraint)

                        with open(restraint_file, "w") as file:
                            # Write the header.
                            file.write("[ position_restraints ]\n")
                            file.write(";  i funct       fcx        fcy        fcz\n")

                            # Write restraints for each atom.
                            for atom_idx in atom_idxs:
                                file.write(f"{atom_idx+1:4}    1       {force_constant}       {force_constant}       {force_constant}\n")

                        # Include the position restraint file in the correct place within
                        # the topology file. We put the additional include directive at the
                        # end of the block so we move to the line before the next moleculetype
                        # record.
                        new_top_lines = top_lines[:moltypes_top_idx[mol_type_idx+1]-1]

                        # Append the additional information.
                        new_top_lines.append('#include "%s"' % restraint_file)
                        new_top_lines.append("")

                        # Now extend with the remainder of the file.
                        new_top_lines.extend(top_lines[moltypes_top_idx[mol_type_idx+1]:])

                        # Overwrite the topology file lines.
                        top_lines = new_top_lines

                        # Increment the number of restraint files.
                        num_restraint += 1

                        # Append the restraint file to the list of autogenerated inputs.
                        self._input_files.append(restraint_file)

                # Write the updated topology to file.
                with open(self._top_file, "w") as file:
                    for line in top_lines:
                        file.write("%s\n" % line)

    def _update_stdout_dict(self):
        """Update the dictionary of thermodynamic records."""

        # Exit if log file hasn't been created.
        if not _os.path.isfile(self._log_file):
            return

        # A list of the new record lines.
        lines = []

        # Append any new lines.
        for line in _pygtail.Pygtail(self._log_file):
            lines.append(line)

        # Store the number of lines.
        num_lines = len(lines)

        # Line index counter.
        x = 0

        # Append any new records to the stdout dictionary.
        while x < num_lines:

            # We've hit any energy record section.
            if lines[x].strip() == "Energies (kJ/mol)":

                # Initialise lists to hold all of the key/value pairs.
                keys = []
                values = []

                # Loop until we reach a blank line, or the end of the lines.
                while True:

                    # End of file.
                    if x + 2 >= num_lines:
                        break

                    # Extract the lines with the keys and values.
                    k_line = lines[x+1]
                    v_line = lines[x+2]

                    # Empty line:
                    if len(k_line.strip()) == 0 or len(v_line.strip()) == 0:
                        break

                    # Add whitespace at the end so that the splitting algorithm
                    # below works properly.
                    k_line = k_line + " "
                    v_line = v_line + " "

                    # Set the starting index of a record.
                    start_idx = 0

                    # Create lists to hold the keys and values.
                    k = []
                    v = []

                    # Split the lines into the record headings and corresponding
                    # values.
                    for idx, val in enumerate(v_line):
                        # We've hit the end of the line.
                        if idx + 1 == len(v_line):
                            break

                        # This is the end of a record, i.e. we've gone from a
                        # character to whitespace. Record the key and value and
                        # update the start index for the next record.
                        if val != " " and v_line[idx+1] == " ":
                            k.append(k_line[start_idx:idx+1])
                            v.append(v_line[start_idx:idx+1])
                            start_idx=idx+1

                    # Update the keys and values, making sure the number of
                    # values matches the number of keys.
                    keys.extend(k)
                    values.extend(v[:len(k)])

                    # Update the line index.
                    x = x + 2

                # Add the records to the dictionary.
                if (len(keys) == len(values)):
                    for key, value in zip(keys, values):
                        # Replace certain characters in the key in order to make
                        # the formatting consistent.

                        # Convert to upper case.
                        key = key.upper()

                        # Strip whitespace and newlines from beginning and end.
                        key = key.strip()

                        # Remove whitespace.
                        key = key.replace(" ", "")

                        # Remove periods.
                        key = key.replace(".", "")

                        # Remove hyphens.
                        key = key.replace("-", "")

                        # Remove parentheses.
                        key = key.replace("(", "")
                        key = key.replace(")", "")

                        # Remove instances of BAR.
                        key = key.replace("BAR", "")

                        # Add the record.
                        self._stdout_dict[key] = value.strip()

            # This is a time record.
            elif "Step" in lines[x].strip():
                if x + 1 < num_lines:
                    records = lines[x+1].split()

                    # There should be two records, 'Step' and 'Time'.
                    if len(records) == 2:
                        self._stdout_dict["STEP"] = records[0].strip()
                        self._stdout_dict["TIME"] = records[1].strip()

                # Update the line index.
                x += 2

            # We've reached an averages section, abort.
            elif " A V E R A G E S" in lines[x]:
                break

            # No match, move to the next line.
            else:
                x += 1

    def _get_stdout_record(self, key, time_series=False, unit=None):
        """Helper function to get a stdout record from the dictionary.

           Parameters
           ----------

           key : str
               The record key.

           time_series : bool
               Whether to return a time series of records.

           unit : BioSimSpace.Types._type.Type
               The unit to convert the record to.

           Returns
           -------

           record :
               The matching stdout record.
        """

        # No data!
        if len(self._stdout_dict) == 0:
            return None

        if not isinstance(time_series, bool):
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Validate the unit.
        if unit is not None:
            if not isinstance(unit, _Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if key == "STEP":
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    if unit is None:
                        return [float(x) for x in self._stdout_dict[key]]
                    else:
                        return [(float(x) * unit)._to_default_unit() for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key == "STEP":
                    return int(self._stdout_dict[key][-1])
                else:
                    if unit is None:
                        return float(self._stdout_dict[key][-1])
                    else:
                        return (float(self._stdout_dict[key][-1]) * unit)._to_default_unit()

            except KeyError:
                return None

    def _getFrame(self, time):
        """Get the trajectory frame closest to a specific time value.

           Parameters
           ----------

           time : :class:`Time <BioSimSpace.Types.Time>`
               The time value.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system from the closest trajectory frame.
        """

        if not isinstance(time, _Types.Time):
            raise TypeError("'time' must be of type 'BioSimSpace.Types.Time'")

        # Grab the last frame from the current trajectory file.
        try:
            with _Utils.cd(self._work_dir):

                # Do we need to get coordinates for the lambda=1 state.
                if "is_lambda1" in self._property_map:
                    is_lambda1 = True
                else:
                    is_lambda1 = False

                # Locate the trajectory file.
                traj_file = self._find_trajectory_file()

                if traj_file is None:
                    return None
                else:
                    self._traj_file = traj_file

                # Use trjconv to get the frame closest to the current simulation time.
                command = "%s trjconv -f %s -s %s -dump %f -pbc mol -o frame.gro" \
                    % (self._exe, self._traj_file, self._tpr_file, time.picoseconds().value())

                # Run the command as a pipeline.
                proc_echo = _subprocess.Popen(["echo", "0"], shell=False, stdout=_subprocess.PIPE)
                proc = _subprocess.Popen(_shlex.split(command), shell=False,
                    stdin=proc_echo.stdout, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
                proc.wait()
                proc_echo.stdout.close()

                # Read the frame file.
                new_system = _IO.readMolecules(["frame.gro", self._top_file],
                                               property_map=self._property_map)

                # Delete the frame file.
                _os.remove("frame.gro")

                # Create a copy of the existing system object.
                old_system = self._system.copy()

                # Update the coordinates and velocities and return a mapping between
                # the molecule indices in the two systems.
                sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                        old_system._sire_object,
                        new_system._sire_object,
                        self._mapping,
                        is_lambda1,
                        self._property_map,
                        self._property_map)

                # Update the underlying Sire object.
                old_system._sire_object = sire_system

                # Store the mapping between the MolIdx in both systems so we don't
                # need to recompute it next time.
                self._mapping = mapping

                # Update the box information in the original system.
                if "space" in new_system._sire_object.propertyKeys():
                    box = new_system._sire_object.property("space")
                    old_system._sire_object.setProperty(self._property_map.get("space", "space"), box)

                return old_system

        except:
            _warnings.warn("Failed to extract trajectory frame with trjconv. "
                           "Try running 'getSystem' again.")
            frame = "%s/frame.gro" % self._work_dir
            if _os.path.isfile(frame):
                _os.remove(frame)
            return None

    def _find_trajectory_file(self):
        """Helper function to find the trajectory file associated with the
           process.

           Returns
           -------

           traj_file : str
               The path to the trajectory file.
        """

        # Check that the current trajectory file is found.
        if not _os.path.isfile(self._traj_file):
            # If not, first check for any trr extension.
            traj_file = _IO.glob("%s/*.trr" % self._work_dir)

            # Store the number of trr files.
            num_trr = len(traj_file)

            # Only accept if a single trajectory file is present.
            if num_trr == 1:
                return traj_file[0]
            else:
                # Now check for any xtc files.
                traj_file = _IO.glob("%s/*.xtc" % self._work_dir)

                if len(traj_file) == 1:
                    return traj_file[0]
                else:
                    _warnings.warn("Invalid trajectory file! "
                                   "%d trr files found, %d xtc files found."
                                   % (num_trr, len(traj_file)))
                    return None
        else:
            return self._traj_file

def _is_minimisation(config):
    """Helper function to check whether a custom configuration
       is a minimisation.

       Parameters
       ----------

       config : [str]
           A list of configuration strings.

       Returns
       -------

       is_minimisation : bool
           Whether this is a minimisation configuration.
    """

    for line in config:
        # Convert to lower-case and remove any whitespace.
        line = line.lower().replace(" ", "")

        # Check for integrators used for minimisation.
        if "integrator=steep" in line or \
           "integrator=cg" in line    or \
           "integrator=l-bfgs" in line:
               return True

    return False

def _is_vacuum(config):
    """Helper function to check whether a configuration corresponds to a
       vacuum simulation.

       Parameters
       ----------

       config : [str]
           A list of configuration strings.

       Returns
       -------

       is_vacuum : bool
           Whether this is (likely) a vacuum configuration.
    """

    # Join all of the config strings together.
    config = "".join(config)

    # Convert to lower-case and remove any whitespace.
    config = config.lower().replace(" ", "")

    # Check for likely indicators of a vacuum simulation.
    # These can be adapted as we think of more, or options
    # change.
    if "pbc=no"     in config and \
       "nstlist=0"  in config and \
       "rlist=0"    in config and \
       "rvdw=0"     in config and \
       "rcoulomb=0" in config:
           return True
    else:
        return False

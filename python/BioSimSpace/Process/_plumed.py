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
Utility class for interfacing with PLUMED.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Plumed"]

import glob as _glob
import os as _os
import pygtail as _pygtail
import shutil as _shutil
import subprocess as _subprocess
import warnings as _warnings

from Sire.Base import findExe as _findExe
from Sire.Mol import MolNum as _MolNum

from BioSimSpace._SireWrappers import System as _System
from BioSimSpace.Metadynamics import CollectiveVariable as _CollectiveVariable
from BioSimSpace.Protocol import Metadynamics as _Metadynamics
from BioSimSpace.Types import Coordinate as _Coordinate

from BioSimSpace import _Exceptions as _Exceptions
from BioSimSpace import Types as _Types
from BioSimSpace import _Utils as _Utils
from BioSimSpace import Units as _Units

from ._process import _MultiDict

class Plumed():
    def __init__(self, work_dir):
        """Constructor.

           Parameters
           ----------

           work_dir : str
               The working directory of the process that is interfacing
               with PLUMED.
        """

        # Check that the working directory is valid.
        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")
        else:
            if not _os.path.isdir(work_dir):
                raise ValueError("'work_dir' doesn't exist: %s" % work_dir)

        # Try to locate the PLUMED executable.
        try:
            self._exe= _findExe("plumed").absoluteFilePath()
        except:
            raise _Exceptions.MissingSoftwareError("Metadynamics simulations required "
                "that PLUMED is installed: www.plumed.org")

        # Run a PLUMED as a background process to query the version number.
        process = _subprocess.run("%s info --version" % self._exe, shell=True, stdout=_subprocess.PIPE)
        plumed_version = float(process.stdout.decode("ascii").strip())

        if plumed_version < 2.5:
            raise _Exceptions.IncompatibleError("PLUMED version >= 2.5 is required.")

        # Set the working directory of the process.
        self._work_dir = work_dir

        # Set the location of the HILLS and COLVAR files.
        self._hills_file = "%s/HILLS" % self._work_dir
        self._colvar_file = "%s/COLVAR" % self._work_dir

        # The number of collective variables.
        self._num_colvar = 0

        # The number of lower/upper walls.
        self._num_lower_walls = 0
        self._num_upper_walls = 0

        # Initialise a list of the collective variable argument names.
        self._colvar_name = []

        # Initalise a dictionary to map the collective variable names
        # to their unit. This can be used when returning time series
        # records from the log files.
        self._colvar_unit = {}

        # Initalise the list of configuration file strings.
        self._config = []

        # Initialise dictionaries to hold COLVAR and HILLS time-series records.
        self._colvar_dict = _MultiDict()
        self._hills_dict = _MultiDict()

        # Initalise lists to store the keys used to index the above dictionary.
        self._colvar_keys = []
        self._hills_keys = []

    def createConfig(self, system, protocol, is_restart=False):
        """Create a PLUMED configuration file.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               A BioSimSpace system object.

           protocol : :class:`Protocol.Metadynamics <BioSimSpace.Protocol.Metadynamics>`
               The metadynamics protocol.

           Returns
           -------

           config : [str]
               The list of PLUMED configuration strings.
        """

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

        if type(protocol) is not _Metadynamics:
            raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.Metadynamics'")

        # Clear data.
        self._num_colvar = 0
        self._num_lower_walls = 0
        self._num_upper_walls = 0
        self._colvar_name = []
        self._colvar_unit = {}
        self._config = []

        # Check for restart files.
        is_restart = False
        if protocol.getHillsFile() is not None:
            is_restart = True
            _shutil.copyfile(protocol.getHillsFile(), "%s/HILLS" % self._work_dir)
        if protocol.getGridFile() is not None:
            is_restart = True
            _shutil.copyfile(protocol.getGridFile(), "%s/GRID" % self._work_dir)
        if protocol.getColvarFile() is not None:
            is_restart = True
            _shutil.copyfile(protocol.getColvarFile(), "%s/COLVAR" % self._work_dir)

        # Is the simulation a restart?
        if is_restart:
            self._config.append("RESTART")
        else:
            self._config.append("RESTART NO")

            # Delete any existing COLVAR, HILLS, and GRID files.
            try:
                _os.remove("%s/COLVAR" % self._work_dir)
                _os.remove("%s/COLVAR.offset" % self._work_dir)
                _os.remove("%s/HILLS" % self._work_dir)
                _os.remove("%s/HILLS.offset" % self._work_dir)
                _os.remove("%s/GRID" % self._work_dir)
            except:
                pass

        # Intialise molecule number to atom tally lookup dictionary in the system.
        try:
            system.getIndex(system[0].getAtoms()[0])
        except:
            raise ValueError("The system contains no molecules?")

        # Store the collective variable(s).
        colvars = protocol.getCollectiveVariable()
        self._num_colvar = len(colvars)

        # Loop over each collective variable and create WHOLEMOLECULES entities
        # for any molecule that involve atoms in a collective variable. We only
        # want to record each molecule once, so keep a list of the molecules
        # that we've already seen.

        molecules = []

        for colvar in colvars:

            # Store all of the atoms to which the collective variable applies.
            atoms = []

            # Distance.
            if type(colvar) is _CollectiveVariable.Distance:
                atom0 = colvar.getAtom0()
                atom1 = colvar.getAtom1()

                if type(atom0) is int:
                    atoms.append(atom0)
                elif type(atom0) is list:
                    atoms.extend(atom0)

                if type(atom1) is int:
                    atoms.append(atom1)
                elif type(atom1) is list:
                    atoms.extend(atom1)

            # Torsion.
            elif type(colvar) is _CollectiveVariable.Torsion:
                atoms = colvar.getAtoms()

            # Loop over all of the atoms. Make sure the index is valid and
            # check if we need to create an entity for the molecule containing
            # the atom.

            for idx in atoms:
                # The atom index is invalid.
                if idx >= system.nAtoms():
                    raise __Exceptions.IncompatibleError("The collective variable is incompatible with the "
                        "system. Contains atom index %d, number of atoms in system is %d " % (idx, system.nAtoms()))

                # Get the molecule numbers in this system.
                mol_nums = system._sire_object.molNums()

                # Loop over each molecule and find the one that contains this atom.
                for x, num in enumerate(mol_nums):
                    # The atom was in the previous molecule.
                    if system._atom_index_tally[num] > idx:
                        num = mol_nums[x-1]
                        break

                # This is a new molecule.
                if num not in molecules:
                    molecules.append(num)

        # Initialise the configuration string.
        string = "WHOLEMOLECULES"

        # Create an entity for each unique molecule.
        for x, molecule in enumerate(molecules):
            # Get the start index.
            idx = system._atom_index_tally[molecule]

            # Get the number of atoms in the molecule.
            num_atoms = system._sire_object.molecule(molecule).nAtoms()

            # Create the entity record. Rember to one-index the atoms.
            string += " ENTITY%d=%d-%d" % (x, idx+1, idx+num_atoms)

        # Append the string to the configuration list.
        self._config.append(string)

        # Intialise tally counters.
        num_distance = 0
        num_torsion = 0
        num_center = 0
        num_fixed = 0

        # Initialise a list to store the grid data for each variable.
        grid_data = []

        # Initialise the METAD string.
        metad_string = "metad: METAD ARG="

        for idx, colvar in enumerate(colvars):

            # Get the lower/upper bounds and the grid data.
            lower_wall = colvar.getLowerBound()
            upper_wall = colvar.getUpperBound()
            grid = colvar.getGrid()

            # Whether the collective variable is a torsion.
            is_torsion = False

            # Distance.
            if type(colvar) is _CollectiveVariable.Distance:
                num_distance += 1

                # Create the argument name.
                arg_name = "d%d" % num_distance

                # Get the atoms between which the distance is measured.
                # Also get the weights, and the center of mass flag.

                # Start.
                atom0 = colvar.getAtom0()
                weights0 = colvar.getWeights0()
                is_com0 = colvar.getCoM0()

                # End.
                atom1 = colvar.getAtom1()
                weights1 = colvar.getWeights1()
                is_com1 = colvar.getCoM1()

                # Initialise the collective variable string.
                colvar_string = "d%d: DISTANCE ATOMS=" % num_distance

                # Process the first atom(s) or fixed coordinate.

                # A single atom.
                if type(atom0) is int:
                    colvar_string += "%d" % (atom0 + 1)

                # A list of atom indices.
                elif type(atom0) is list:
                    num_center += 1
                    colvar_string += "c%d" % num_center

                    center_string = "c%d: CENTER ATOMS=%s" \
                        % (num_center, ",".join([str(x+1) for x in atom0]))

                    # Center of mass weighting takes precendence.
                    if is_com0:
                        center_string += " MASS"

                    # User weights.
                    elif weights0 is not None:
                        center_string += " WEIGHTS=%s" % ",".join([str(x) for x in weights0])

                    self._config.append(center_string)

                # A coordinate of a fixed point.
                elif type(atom0) is _Coordinate:
                    # Convert to nanometers.
                    x = atom0.x().nanometers().magnitude()
                    y = atom0.y().nanometers().magnitude()
                    z = atom0.z().nanometers().magnitude()

                    num_fixed += 1
                    self._config.append("f%d: FIXEDATOM AT=%s,%s,%s" % (num_fixed, x, y, z))
                    colvar_string += "f%d" % num_fixed

                # Process the second atom(s) or fixed coordinate.

                # A single atom.
                if type(atom1) is int:
                    colvar_string += ",%d" % (atom1 + 1)

                # A list of atom indices.
                elif type(atom1) is list:
                    num_center += 1
                    colvar_string += ",c%d" % num_center

                    center_string = "c%d: CENTER ATOMS=%s" \
                        % (num_center, ",".join([str(x+1) for x in atom1]))

                    # Center of mass weighting takes precendence.
                    if is_com1:
                        center_string += " MASS"

                    # User weights.
                    elif weights1 is not None:
                        center_string += " WEIGHTS=%s" % ",".join([str(x) for x in weights1])

                    self._config.append(center_string)

                # A coordinate of a fixed point.
                elif type(atom1) is _Coordinate:
                    # Convert to nanometers.
                    x = atom1.x().nanometers().magnitude()
                    y = atom1.y().nanometers().magnitude()
                    z = atom1.z().nanometers().magnitude()

                    num_fixed += 1
                    self._config.append("f%d: FIXEDATOM AT=%s,%s,%s" % (num_fixed, x, y, z))
                    colvar_string += ",f%d" % num_fixed

                # Disable periodic boundaries.
                if not colvar.getPeriodicBoundaries():
                    colvar_string += " NOPBC"

                # Measure the x, y, and z distance components.
                if colvar.getComponent() is not None:
                    colvar_string += " COMPONENT"
                    arg_name += ".%s" % colvar.getComponent()

                # Append the collective variable record.
                self._config.append(colvar_string)

                # Store the collective variable name and its unit.
                self._colvar_name.append(arg_name)
                self._colvar_unit[arg_name] = _Units.Length.nanometer

            # Torsion.
            elif type(colvar) is _CollectiveVariable.Torsion:
                is_torsion = True
                num_torsion += 1
                arg_name = "t%d" % num_torsion
                colvar_string = "%s: TORSION ATOMS=%s" \
                    % (arg_name, ",".join([str(x+1) for x in colvar.getAtoms()]))

                # Store the collective variable name and its unit.
                self._colvar_name.append(arg_name)
                self._colvar_unit[arg_name] = _Units.Angle.radian

                # Disable periodic boundaries.
                if not colvar.getPeriodicBoundaries():
                    colvar_string += " NOPBC"

                # Append the collective variable record.
                self._config.append(colvar_string)

            # Check for lower and upper bounds on the collective variable.
            if lower_wall is not None:
                self._num_lower_walls += 1
                lower_wall_string = "lwall%d: LOWER_WALLS ARG=%s" % (self._num_lower_walls, arg_name)
                try:
                    # Unit based.
                    lower_wall_string += ", AT=%s" % lower_wall.getValue().magnitude()
                except:
                    # Dimensionless.
                    lower_wall_string += ", AT=%s" % lower_wall.getValue()
                lower_wall_string += ", KAPPA=%s" % lower_wall.getForceConstant()
                lower_wall_string += ", EXP=%s" % lower_wall.getExponent()
                lower_wall_string += ", EPS=%s" % lower_wall.getEpsilon()
                self._config.append(lower_wall_string)

            # Check for lower and upper bounds on the collective variable.
            if upper_wall is not None:
                self._num_upper_walls += 1
                upper_wall_string = "uwall%d: UPPER_WALLS ARG=%s" % (self._num_upper_walls, arg_name)
                try:
                    # Unit based.
                    upper_wall_string += ", AT=%s" % upper_wall.getValue().magnitude()
                except:
                    # Dimensionless.
                    upper_wall_string += ", AT=%s" % upper_wall.getValue()
                upper_wall_string += ", KAPPA=%s" % upper_wall.getForceConstant()
                upper_wall_string += ", EXP=%s" % upper_wall.getExponent()
                upper_wall_string += ", EPS=%s" % upper_wall.getEpsilon()
                self._config.append(upper_wall_string)

            # Store grid data.
            if grid is not None:
                if is_torsion:
                        grid_data.append(("-pi", "pi", grid.getBins()))
                else:
                    try:
                        # Unit based.
                        grid_data.append((grid.getMinimum().magnitude(),
                                          grid.getMaximum().magnitude(),
                                          grid.getBins()))
                    except:
                        # Dimensionless.
                        grid_data.append((grid.getMinimum(),
                                          grid.getMaximum(),
                                          grid.getBins()))

            # Add the argument to the METAD record.
            metad_string += "%s" % arg_name

            # Update the METAD record to separate the collective variable arguments.
            if idx < self._num_colvar - 1:
                metad_string += ","

        # Now complete the METAD record string.

        # Hill width.
        metad_string += " SIGMA="
        for idx, colvar in enumerate(colvars):
            metad_string += "%s" % colvar.getHillWidth().magnitude()
            if idx < self._num_colvar - 1:
                metad_string += ","

        # Hill height.
        metad_string += " HEIGHT=%s" % protocol.getHillHeight().magnitude()

        # Hill frequency.
        metad_string += " PACE=%s" % protocol.getHillFrequency()

        # Grid parameters.
        if len(grid_data) > 0:
            grid_min_string = " GRID_MIN="
            grid_max_string = " GRID_MAX="
            grid_bin_string = " GRID_BIN="

            for idx, grid in enumerate(grid_data):
                grid_min_string += str(grid[0])
                grid_max_string += str(grid[1])
                grid_bin_string += str(grid[2])
                if idx < self._num_colvar - 1:
                    grid_min_string += ","
                    grid_max_string += ","
                    grid_bin_string += ","

            metad_string += grid_min_string + grid_max_string + grid_bin_string
            metad_string += " GRID_WFILE=GRID GRID_WSTRIDE=%s" % protocol.getHillFrequency()
            if protocol.getGridFile() is not None:
                metad_string += " GRID_RFILE=GRID"
            metad_string += " CALC_RCT"

        # Temperature and bias parameters.
        metad_string += " TEMP=%s" % protocol.getTemperature().kelvin().magnitude()
        if protocol.getBiasFactor() is not None:
            metad_string += " BIASFACTOR=%s" % protocol.getBiasFactor()

        # Append the METAD record to the config.
        self._config.append(metad_string)

        # Print all record data to the COLVAR file.
        print_string = "PRINT STRIDE=%s ARG=* FILE=COLVAR" % protocol.getHillFrequency()
        self._config.append(print_string)

        return self._config

    def getTime(self, time_series=False):
        """Get the simulation run time.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           time : :class:`Time <BioSimSpace.Types.Time>`
               The simulation run time.
        """
        # Get the latest records.
        self._update_colvar_dict()

        # Get the corresponding data from the dictionary and return.
        return self._get_colvar_record(key="time",
                                       time_series=time_series,
                                       unit=_Units.Time.picosecond)

    def getCollectiveVariable(self, index, time_series=False):
        """Get the value of a collective variable.

           Parameters
           ----------

           index : int
               The index of the collective variable.

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           collective_variable : :class:`Type <BioSimSpace.Types>`
               The value of the collective variable.
        """
        if type(index) is not int:
            raise TypeError("'index' must be of type 'int'")
        if index > self._num_colvar - 1 or index < -self._num_colvar:
            raise IndexError("'index' must be in range -%d to %d" % (self._num_colvar, self._num_colvar-1))

        # Get the latest records.
        self._update_colvar_dict()

        # Get the corresponding data from the dictionary and return.
        return self._get_colvar_record(key=self._colvar_name[index],
                                       time_series=time_series,
                                       unit=self._colvar_unit[self._colvar_name[index]])

    def getFreeEnergy(self, index=None, stride=None, kt=None):
        """Get the current free energy estimate.

           Parameters
           ----------

           index : int
               The index of the collective variable. If None, then all variables
               will be considered.

           stride : int
               The stride for integrating the free energy. This can be used to
               check for convergence.

           kt : BioSimSpace.Types.Energy
               The temperature in energy units for intergrating out variables.

           free_energies : [:class:`Type <BioSimSpace.Types>`, ...], \
                           [[:class:`Type <BioSimSpace.Types>`, :class:`Type <BioSimSpace.Types>`, ...], ...]
               The free energy estimate for the chosen collective variables.
        """
        if index is not None:
            if type(index) is not int:
                raise TypeError("'index' must be of type 'int'")
            if index > self._num_colvar - 1 or index < -self._num_colvar:
                raise IndexError("'index' must be in range -%d to %d" % (self._num_colvar, self._num_colvar-1))

        if stride is not None:
            if type(stride) is not int:
                raise TypeError("'stride' must be of type 'int'")
            if stride < 0:
                raise ValueError("'stride' must be >= 0")

        if kt is not None:
            if type(kt) is not _Types.Energy:
                raise TypeError("'kt' must be of type 'BioSimSpace.Type.Energy'")

            # Convert to kt and get the magnitude.
            kt = kt.kt().magnitude()

            if kt <= 0:
                raise ValueError("'kt' must have magnitude > 0")
        else:
            if index is not None:
                raise ValueError("You must specify 'kt' when making a dimensionality reduction.")

        # Create the command string.
        command = "%s sum_hills --hills HILLS --mintozero" % self._exe

        # Append additional arguments.
        if index is not None:
            command += " --idw %s" % self._colvar_name[index]
            command += " --kt %s" % kt
        if stride is not None:
            command += " --stride %s" % stride

        # Initialise a list to hold the free energy estimates.
        free_energies = []

        # Move to the working directory.
        with _Utils.cd(self._work_dir):

            # Run the sum_hills command as a background process.
            proc = _subprocess.run(command, shell=True,
                stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

            if proc.returncode != 0:
                raise RuntimeError("Failed to generate free energy estimate.\n"
                                   "Error: %s" % proc.stderr.decode("utf-8"))

            # Get a sorted list of all the fes*.dat files.
            fes_files = _glob.glob("fes*.dat")
            fes_files.sort()

            # Process each of the files.
            for fes in fes_files:
                # Create a list to store the free energy estimate for this file.
                free_energy = []

                if index is None:
                    # Create lists for each the collective variable data point.
                    for x in range(0, self._num_colvar):
                        free_energy.append([])
                else:
                    free_energy.append([])
                # Create a list to store the free energy.
                free_energy.append([])

                # Read the file.
                with open(fes, "r") as file:

                    # Loop over all lines in the file.
                    for line in file:

                        # Ignore comments and blank lines.
                        if line[0] != "#":

                            # Extract the data.
                            # This is: colvar1, colvar2, ..., fes
                            data = [float(x) for x in line.split()]

                            # The line contains data.
                            if len(data) > 0:

                                # Store data for each of the collective variables.
                                if index is None:
                                    for x in range(0, self._num_colvar):
                                        name = self._colvar_name[x]
                                        free_energy[x].append(data[x] * self._colvar_unit[name])
                                    free_energy[self._num_colvar].append(data[self._num_colvar] * _Units.Energy.kj_per_mol)
                                else:
                                    name = self._colvar_name[0]
                                    free_energy[0].append(data[0] * self._colvar_unit[name])
                                    free_energy[1].append(data[1] * _Units.Energy.kj_per_mol)

                if len(fes_files) == 1:
                    free_energies = free_energy
                else:
                    free_energies.append(tuple(free_energy))

                # Remove the file.
                _os.remove(fes)

        return tuple(free_energies)

    def _update_colvar_dict(self):
        """Read the COLVAR file and update any records."""

        # Exit if the COLVAR file hasn't been created.
        if not _os.path.isfile(self._colvar_file):
            return

        # Loop over all new lines in the file.
        for line in _pygtail.Pygtail(self._colvar_file):

            # Is this a header line. If so, store the keys.
            if line[3:9] == "FIELDS":
                self._colvar_keys = line[10:].split()

            # This is an actual data record. Update the multi-dictionary.
            elif line[0] != "#":
                data = [float(x) for x in line.split()]
                for key, value in zip(self._colvar_keys, data):
                    self._colvar_dict[key] = value

    def _update_hills_dict(self):
        """Read the HILLS file and update any records."""

        # Exit if the HILLS file hasn't been created.
        if not _os.path.isfile(self._hills_file):
            return

        # Loop over all new lines in the file.
        for line in _pygtail.Pygtail(self._hills_file):

            # Is this a header line. If so, store the keys.
            if line[3:9] == "FIELDS":
                self._hills_keys = line[10:].split()

            # This is an actual data record. Update the multi-dictionary.
            elif line[0] != "#":
                data = [float(x) for x in line.split()]
                for key, value in zip(self._hills_keys, data):
                    self._hills_dict[key] = value

    def _get_colvar_record(self, key, time_series=False, unit=None):
        """Helper function to get a COLVAR record from the dictionary.

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
        if len(self._colvar_dict) is 0:
            return None

        if type(time_series) is not bool:
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Valdate the unit.
        if unit is not None:
            if not isinstance(unit, _Types._type.Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if unit is None:
                    return self._colvar_dict[key]
                else:
                    if key == "time":
                        return [(x * unit).nanoseconds() for x in self._colvar_dict[key]]
                    else:
                        return [x * unit for x in self._colvar_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if unit is None:
                    return self._colvar_dict[key][-1]
                else:
                    if key == "time":
                        return (self._colvar_dict[key][-1] * unit).nanoseconds()
                    else:
                        return self._colvar_dict[key][-1] * unit

            except KeyError:
                return None

    def _get_hills_record(self, key, time_series=False, unit=None):
        """Helper function to get a HILLS record from the dictionary.

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
        if len(self._hills_dict) is 0:
            return None

        if type(time_series) is not bool:
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Valdate the unit.
        if unit is not None:
            if not isinstance(unit, _Types._type.Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if unit is None:
                    return self._hills_dict[key]
                else:
                    return [x * unit for x in self._hills_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if unit is None:
                    return self._hills_dict[key][-1]
                else:
                    return self._hills_dict[key][-1] * unit

            except KeyError:
                return None

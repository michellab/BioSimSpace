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

import os as _os
import subprocess as _subprocess

from Sire.Base import findExe as _findExe
from Sire.Mol import MolNum as _MolNum

from .._SireWrappers import System as _System
from ..Metadynamics import CollectiveVariable as _CollectiveVariable
from ..Protocol import Metadynamics as _Metadynamics
from ..Types import Coordinate as _Coordinate

import BioSimSpace._Exceptions as _Exceptions
import BioSimSpace.Units as _Units

class Plumed():
    def __init__(self):
        """Constructor."""

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

        # The number of collective variables.
        self._num_colvar = 0

        # The number of lower/upper walls.
        self._num_lower_walls = 0
        self._num_upper_walls = 0

        # Initialise a list of the collective variable argument names.
        self._colvar_names = []

        # Initalise a dictionary to map the collective variable names
        # to their unit. This can be used when returning time series
        # records from the log files.
        self._colvar_unit = {}

        # Initalise the list of configuration file strings.
        self._config = []

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
        self._colvar_names = []
        self._colvar_unit = {}
        self._config = []

        # Is the simulation a restart?
        if protocol.isRestart():
            self._config.append("RESTART")
        else:
            self._config.append("RESTART NO")

        # Intialise molecule number to atom tally lookup dictionary in the system.
        try:
            system.getIndex(system.getMolecules()[0].getAtoms()[0])
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
                self._colvar_names.append(arg_name)
                self._colvar_unit[arg_name] = _Units.Length.nanometer

                # Check for lower and upper bounds on the collective variable.
                if lower_wall is not None:
                    self._num_lower_walls += 1
                    lower_wall_string = "lwall%d: LOWER_WALLS ARG=%s" % (self._sum_lower_walls, arg_name)
                    lower_wall_string += ", AT=%s" % lower_wall.getValue().nanometers().magnitude()
                    lower_wall_string += ", KAPPA=%s" % lower_wall.getForceConstant()
                    lower_wall_string += ", EXP=%s" % lower_wall.getExponent()
                    lower_wall_string += ", EPS=%s" % lower_wall.getEpsilon()
                    self._config.append(lower_wall_string)

                # Check for lower and upper bounds on the collective variable.
                if upper_wall is not None:
                    self._num_upper_walls += 1
                    upper_wall_string = "uwall%d: UPPER_WALLS ARG=%s" % (self._num_upper_walls, arg_name)
                    upper_wall_string += ", AT=%s" % upper_wall.getValue().nanometers().magnitude()
                    upper_wall_string += ", KAPPA=%s" % upper_wall.getForceConstant()
                    upper_wall_string += ", EXP=%s" % upper_wall.getExponent()
                    upper_wall_string += ", EPS=%s" % upper_wall.getEpsilon()
                    self._config.append(upper_wall_string)

                # Store grid data.
                if grid is not None:
                    grid_data.append((grid.getMinimum().nanometers().magnitude(),
                                      grid.getMaximum().nanometers().magnitude(),
                                      grid.getBins()))

            # Torsion.
            elif type(colvar) is _CollectiveVariable.Torsion:
                num_torsion += 1
                arg_name = "t%d" % num_torsion
                colvar_string = "%s: TORSION ATOMS=%s" \
                    % (arg_name, ",".join([str(x+1) for x in colvar.getAtoms()]))

                # Store the collective variable name and its unit.
                self._colvar_names.append(arg_name)
                self._colvar_unit[arg_name] = _Units.Length.nanometer

                # Disable periodic boundaries.
                if not colvar.getPeriodicBoundaries():
                    colvar_string += " NOPBC"

                # Append the collective variable record.
                self._config.append(colvar_string)

                # Check for lower and upper bounds on the collective variable.
                if lower_wall is not None:
                    self._num_lower_walls += 1
                    lower_wall_string = "lwall%d: LOWER_WALLS ARG=%s" % (self._num_lower_walls, arg_name)
                    lower_wall_string += ", AT=%s" % lower_wall.getValue().radians().magnitude()
                    lower_wall_string += ", KAPPA=%s" % lower_wall.getForceConstant()
                    lower_wall_string += ", EXP=%s" % lower_wall.getExponent()
                    lower_wall_string += ", EPS=%s" % lower_wall.getEpsilon()
                    self._config.append(lower_wall_string)

                # Check for lower and upper bounds on the collective variable.
                if upper_wall is not None:
                    self._num_upper_walls += 1
                    upper_wall_string = "uwall%d: UPPER_WALLS ARG=%s" % (self._num_upper_walls, arg_name)
                    upper_wall_string += ", AT=%s" % upper_wall.getValue().radians().magnitude()
                    upper_wall_string += ", KAPPA=%s" % upper_wall.getForceConstant()
                    upper_wall_string += ", EXP=%s" % upper_wall.getExponent()
                    upper_wall_string += ", EPS=%s" % upper_wall.getEpsilon()
                    self._config.append(upper_wall_string)

                # Store grid data.
                if grid is not None:
                    grid_data.append((grid.getMinimum().radians().magnitude(),
                                      grid.getMaximum().radians().magnitude(),
                                      grid.getBins()))

            # Add the argument to the METAD record.
            metad_string += "%s" % arg_name

            # Update the METAD record to separate the collective variable arguments.
            if idx < self._num_colvar - 1:
                metad_string += ","

        # Now complete the METAD record string.

        # Hill width.
        metad_string += " SIGMA="
        for idx, width in enumerate(protocol.getHillWidth()):
            metad_string += "%s" % width
            if idx < self._num_colvar - 1:
                metad_string += ","

        # Hill height.
        metad_string += " HEIGHT=%s" % protocol.getHillHeight()

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

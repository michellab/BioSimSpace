######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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
Functionality for handling parameterisation protocols.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from BioSimSpace import _amber_home, _gmx_exe
from ..._SireWrappers import Molecule as _Molecule

import BioSimSpace.IO as _IO

import os as _os
import queue as _queue
import subprocess as _subprocess

__all__ = ["Protocol"]

# Set the tLEaP cmd directory.
_cmd_dir = _os.path.dirname(_amber_home) + "/dat/leap/cmd"

# Search for all of the required executables.

# Search for the tLEaP exe.

_tleap_exe = None

if _amber_home is not None:
    _exe = "%s/bin/tleap" % _amber_home
    if _os.path.isfile(_exe):
        _tleap_exe = _exe
    else:
        raise IOError("Missing tLEaP executable: '%s'" % _exe)

# Search for the Antechamber exe.

_antechamber_exe = None

if _amber_home is not None:
    _exe = "%s/bin/antechamber" % _amber_home
    if _os.path.isfile(_exe):
        _antechamber_exe = _exe
    else:
        raise IOError("Missing Antechamber executable: '%s'" % _exe)

# Search for the parmchk exe.

_parmchk_exe = None

if _amber_home is not None:
    _exe = "%s/bin/parmchk2" % _amber_home
    if _os.path.isfile(_exe):
        _parmchk_exe = _exe
    else:
        raise IOError("Missing parmchk executable: '%s'" % _exe)

class Protocol():
    """A base class for parameterisation protocols."""

    def __init__(self, forcefield, map={}):
        """Constructor.

           Positional arguments:

           forcefield -- The name of the force field.

           Keyword arguments:

           map        -- A dictionary that maps system "properties" to their user defined
                         values. This allows the user to refer to properties with their
                         own naming scheme, e.g. { "charge" : "my-charge" }
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Validate and set the force field name.
        if type(forcefield) is not str:
            raise TypeError("'forcefield' must be of type 'str'")
        else:
            self._forcefield = forcefield

        # Validate and set the property map.
        if type(map) is not dict:
            raise TypeError("'map' must be of type 'dict'")
        else:
            self._map = map.copy()

        # Set default compatability flags for the utility programs.
        # These must be overridden in the constructor of any derived classes.
        self._tleap = False
        self._pdb2gmx = False

    def run(self, molecule, work_dir=None, queue=None):
        """Run the parameterisation protocol.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
           work_dir -- The working directory.
           queue    -- The thread queue is which this method has been run.
        """

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

        if type(work_dir) is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        if type(queue) is not None and type(queue) is not _queue.Queue:
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Store the current working directory.
        dir = _os.getcwd()

        # Set up the working directory.
        if work_dir is not None:
            # Create the working directory, if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir)

            # Change to the working directory for the process.
            # This avoid problems with relative paths.
            _os.chdir(work_dir)

        # Create a new molecule using a deep copy of the internal Sire Molecule.
        new_mol = _Molecule(molecule._getSireMolecule().__deepcopy__())

        # Choose the program to run with depending on the force field compatibility.
        # If tLEaP and pdb2gmx are supported, default to tLEaP, then use pdb2gmx if
        # tLEaP fails to produce output.

        # Parameterise using tLEaP.
        if self._tleap:
            if _tleap_exe is not None:
                output = self._run_tleap(molecule)
            else:
                output = None

            # If there was no output, try using pdbgmx
            if output is None:
                if _gmx_exe:
                    output = self._run_pdb2gmx(molecule)
                else:
                    output = None

        # Parameterise using pdb2gmx.
        else:
            if _gmx_exe:
                output = self._run_pdb2gmx(molecule)
            else:
                return None

        # Change back to the original directory.
        if work_dir is not None:
            _os.chdir(dir)

        # No output, parameterisation failed.
        if output is None:
            if queue is not None:
                queue.put(None)
            return None

        else:
            # Prepend the working directory to the output file names.
            if work_dir is not None:
                output = ["%s/%s" % (work_dir, output[0]),
                          "%s/%s" % (work_dir, output[1])]

            # Load the parameterised molecule.
            par_mol = _Molecule(_IO.readMolecules(output)._getSireSystem()[_Sire.Mol.MolIdx(0)])

            # Make the molecule 'mol' compatible with 'par_mol'. This will create
            # a mapping between atom indices in the two molecules and add all of
            # the new properties from 'par_mol' to 'mol'.
            new_mol._makeCompatibleWith(par_mol, map=self._map, overwrite=True, verbose=False)

            if queue is not None:
                queue.put(new_mol)
            return new_mol

    def _run_tleap(self, molecule):
        """Run using tLEaP.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
        """

        # Create a new system and molecule group.
        s = _Sire.System.System("BioSimSpace System")
        m = _Sire.Mol.MoleculeGroup("all")

        # Add the molecule.
        m.add(molecule._getSireMolecule())
        s.add(m)

        # Write the system to a PDB file.
        pdb = _Sire.IO.PDB2(s)
        pdb.writeToFile("leap.pdb")

        # Try to find a force field file.
        ff = _find_force_field(self._forcefield)

        # Write the LEaP input file.
        with open("leap.txt", "w") as f:
            f.write("source %s\n" % ff)
            f.write("mol = loadPdb leap.pdb\n")
            f.write("saveAmberParm mol leap.top leap.crd\n")
            f.write("quit")

        # Generate the tLEaP command.
        command = "%s -f leap.txt" % _tleap_exe

        with open("README.txt", "w") as f:
            # Write the command to file.
            f.write("# tLEaP was run with the following command:\n")
            f.write("%s\n" % command)

        # Run tLEaP as a subprocess.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # tLEaP doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile("leap.top") and _os.path.isfile("leap.crd"):
            return ["leap.top", "leap.crd"]
        else:
            return None

    def _run_pdb2gmx(self, molecule):
        """Run using pdb2gmx.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
        """

        # A list of supported force fields, mapping to their GROMACS ID string.
        # GROMACS supports a sub-set of the AMBER force fields.
        supported_ff = { "ff99"   : "amber99",
                         "ff99SB" : "amber99sb",
                         "ff03"   : "amber03"
                       }

        if self._forcefield not in supported_ff:
            raise ValueError("'pdb2gmx' does not support the '%s' force field." % self._forcefield)

        # Create a new system and molecule group.
        s = _Sire.System.System("BioSimSpace System")
        m = _Sire.Mol.MoleculeGroup("all")

        # Add the molecule.
        m.add(molecule._getSireMolecule())
        s.add(m)

        # Write the system to a PDB file.
        pdb = _Sire.IO.PDB2(s)
        pdb.writeToFile("input.pdb")

        # Generate the pdb2gmx command.
        command = "%s pdb2gmx -f input.pdb -o output.gro -p output.top -ignh -ff %s -water none" \
            % (_gmx_exe, supported_ff[self._forcefield])

        with open("README.txt", "w") as f:
            # Write the command to file.
            f.write("# pdb2gmx was run with the following command:\n")
            f.write("%s\n" % command)

        # Run pdb2gmx as a subprocess.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Check for the expected output.
        if _os.path.isfile("output.gro") and _os.path.isfile("output.top"):
            return ["output.gro", "output.top"]
        else:
            return None

def _find_force_field(forcefield):
    """Internal function to search LEaP compatible force field files."""

    # Whether the force field is found.
    is_found = False

    # Whether the force field is old.
    is_old = False

    # Search for a compatible force field file.
    if _amber_home is not None:
        ff = _IO.glob("%s/dat/leap/cmd/*.%s" % (_amber_home, forcefield))

        # Search the old force fields. First try a specific match.
        if len(ff) == 0:
            ff = _IO.glob("%s/dat/leap/cmd/oldff/leaprc.s" % (_amber_home, forcefield))

            # No matches, try globbing all files with matching extension.
            if len(ff) == 0:
                ff = _IO.glob("%s/oldff/*.%s" % (_cmd_dir, forcefield))

            if len(ff) > 0:
                is_found = True
                is_old = True

    if not is_found:
        ff = _IO.glob("%s/*.%s" % (_cmd_dir, forcefield))

        # Search the old force fields. First try a specific match.
        if len(ff) == 0:
            ff = _IO.glob("%s/oldff/leaprc.%s" % (_cmd_dir, forcefield))
            is_old = True

            # No matches, try globbing all files with matching extension.
            if len(ff) == 0:
                ff = _IO.glob("%s/oldff/*.%s" % (_cmd_dir, self._forcefield))
                # No force field found!
                if len(ff) == 0:
                    raise ValueError("No force field file found for '%s'" % self._forcefield)

    # Multiple force fields found.
    if len(ff) > 1:
        raise ValueError("Multiple force fields found for '%s': %s" % (self._forcefield, ff))

    # Create the force field name.
    ff = _os.path.basename(ff[0])
    if is_old:
        ff = "oldff/" + ff

    # Return the force field.
    return ff

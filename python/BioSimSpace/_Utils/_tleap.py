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
Functionality for parameterising molecules using LEaP.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from BioSimSpace import _amber_home

import BioSimSpace.IO as _IO

import os as _os
import subprocess as _subprocess
import tempfile as _tempfile

__all__ = ["Tleap"]

# Set the bundled tLEaP cmd directory.
_cmd_dir = _os.path.dirname(_Sire.Base.getBinDir()) + "/dat/leap/cmd"

# Search for the tLEaP exe.

_exe = None

if _amber_home is not None:
    if _os.path.isfile("%s/bin/tleap" % _amber_home):
        _exe = "%s/bin/tleap" % _amber_home

if _exe is None:
    # Search Sire bin directory.
    bin_dir = _Sire.Base.getBinDir()
    exe = "%s/tleap" % bin_dir

    if _os.path.isfile(exe):
        _exe = exe

    # Search system PATH.
    else:
        _exe = _Sire.Base.findExe("tleap").absoluteFilePath()

class Tleap():
    """A simple class for parameterising molecules using LEaP."""

    def run(molecule, forcefield):
        """Positional arguments:

           molecule   -- The molecule to parameterise.
           forcefield -- The forcefield to use.
        """

        if type(molecule) is not _Sire.Mol._Mol.Molecule:
            raise TypeError("'molecule' must be of type 'Sire.Mol._Mol.Molecule'")

        if type(forcefield) is not str:
            raise TypeError("'forcefield' must be of type 'str'")

        # Whether the force field was found.
        is_found = False

        # Whether the force field is old.
        is_old = False

        # Search for a compatible force field file.
        if _amber_home is not None:
            ff = _IO.glob("%s/dat/leap/cmd/*.%s" % (_amber_home, forcefield))

            # Search the old force fields
            if len(ff) == 0:
                ff = _IO.glob("%s/dat/leap/cmd/oldff/*.%s" % (_amber_home, forcefield))

                if len(ff) > 0:
                    is_found = True
                    is_old = True

        if not is_found:
            ff = _IO.glob("%s/*.%s" % (_cmd_dir, forcefield))

            # Search the old force fields
            if len(ff) == 0:
                ff = _IO.glob("%s/oldff/*.%s" % (_cmd_dir, forcefield))

                # No force field found!
                if len(ff) == 0:
                    raise ValueError("No force field file found for '%s'" % forcefield)

        # Multiple force fields found.
        if len(ff) > 1:
            raise ValueError("Multiple force fields found for '%s': %s" % (forcefield, ff))

        # Create the force field name.
        ff = _os.path.basename(ff[0])
        if is_old:
            ff = "oldff/" + ff

        # Create a new system and molecule group.
        s = _Sire.System.System("BioSimSpace molecule")
        m = _Sire.Mol.MoleculeGroup("all")

        # Add the molecule.
        m.add(molecule)
        s.add(m)

        # Generate a temporary directory and store its name.
        tmp_dir = _tempfile.TemporaryDirectory()
        work_dir = tmp_dir.name

        # Store the current working directory.
        dir = _os.getcwd()

        # Change to the working directory.
        _os.chdir(work_dir)

        # Create a list of the expected output files.
        output = ["%s/tmp.top" % work_dir, "%s/tmp.crd" % work_dir]

        # Write the system to a PDB file.
        pdb = _Sire.IO.PDB2(s)
        pdb.writeToFile("tmp.pdb")

        # Write the LEaP input file.
        with open("tmp.txt", "w") as f:
            f.write("source %s\n" % ff)
            f.write("tmp = loadPdb tmp.pdb\n")
            f.write("saveAmberParm tmp tmp.top tmp.crd\n")
            f.write("quit")

        # Generate the tLEaP command.
        command = "%s -f tmp.txt" % exe

        # Run tLEaP as a subprocess.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Change back to the current directory.
        _os.chdir(dir)

        # tLEaP doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile(output[0]) and _os.path.isfile(output[1]):

            # Return the parameterised molecule.
            return _Sire.IO.MoleculeParser.read(output)[_Sire.Mol.MolIdx(0)]

            # Things to check:
            #    1) Number of atoms is consistent.
            #    2) Leap log file does not specify addition or removal of atoms.
            #    3) Funky work necessary for HIS, HIE, HID, HIP (histidene).
            #    4) More funky logic for CYS/CYX, including having to add addBond lines.

        else:
            raise RuntimeError("tLEaP failed to generate the expected output!")

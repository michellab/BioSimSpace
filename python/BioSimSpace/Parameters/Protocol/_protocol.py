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

from BioSimSpace import _amber_home, _bin_dir
from ..._SireWrappers import Molecule as _Molecule

import BioSimSpace.IO as _IO

import os as _os
import subprocess as _subprocess

__all__ = ["Protocol"]

# Set the bundled tLEaP cmd directory.
_cmd_dir = _os.path.dirname(_bin_dir) + "/dat/leap/cmd"

# Search for all of the required executables.

# Search for the tLEaP exe.

_tleap_exe = None

if _amber_home is not None:
    if _os.path.isfile("%s/bin/tleap" % _amber_home):
        _tleap_exe = "%s/bin/tleap" % _amber_home

if _tleap_exe is None:
    _tleap_exe = "%s/tleap" % _bin_dir

    if not _os.path.isfile(_tleap_exe):
        _tleap_exe = _Sire.Base.findExe("tleap").absoluteFilePath()

# Search for the Antechamber exe.

_antechamber_exe = None

if _amber_home is not None:
    if _os.path.isfile("%s/bin/antechamber" % _amber_home):
        _antechamber_exe = "%s/bin/antechamber" % _amber_home

if _antechamber_exe is None:
    _antechamber_exe = "%s/antechamber" % _bin_dir

    if not _os.path.isfile(_antechamber_exe):
        _antechamber_exe = _Sire.Base.findExe("antechamber").absoluteFilePath()

# Search for the parmchk exe.

_parmchk_exe = None

if _amber_home is not None:
    if _os.path.isfile("%s/bin/parmchk2" % _amber_home):
        _parmchk_exe = "%s/bin/parmchk2" % _amber_home

if _parmchk_exe is None:
    _parmchk_exe = "%s/parmchk2" % _bin_dir

    if not _os.path.isfile(_parmchk_exe):
        _parmchk_exe = _Sire.Base.findExe("parmchk2").absoluteFilePath()

# Search for the gmx exe.

_gmx_exe = "%s/gmx" % _bin_dir
if not _os.path.isfile(_gmx_exe):
    _gmx_exe = _Sire.Base.findExe("gmx").absoluteFilePath()

class Protocol():
    """A base class for parameterisation protocols."""

    def __init__(self, forcefield):
        """Constructor."""

	# Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Validate and set the force field name.
        if type(forcefield) is not str:
            raise TypeError("'forcefield' must be of type 'str'")
        else:
            self._forcefield = forcefield

        # Set default compatability flags for the utility programs.
        # These must be overridden in the constructor of any derived classes.
        self._tleap = False
        self._pdb2gmx = False

    def run(self, molecule):
        """Run the parameterisation protocol.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
        """

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

        # Create a new molecule using a deep copy of the internal Sire Molecule.
        new_mol = _Molecule(molecule._molecule.__deepcopy__())

        # Choose the program to run with depending on the force field compatibility.
        # If tLEaP and pdb2gmx are supported, default to tLEaP, then use pdb2gmx if
        # tLEaP fails to produce output.

        # Parameterise using tLEaP.
        if self._tleap:
            output = self._run_tleap(molecule)

            # If there was no output, try using pdbgmx
            if output is None:
                output = self._run_pdb2gmx(molecule)

        # Parameterise using pdb2gmx.
        else:
            output = self._run_pdb2gmx(molecule)

        # No output, parameterisation failed.
        if output is None:
            return None

        else:
            # Load the parameterised molecule.
            par_mol = _Molecule(_Sire.IO.MoleculeParser.read(output)[_Sire.Mol.MolIdx(0)])

            # Make the molecule 'mol' compatible with 'par_mol'. This will create
            # a mapping between atom indices in the two molecules and add all of
            # the new properties from 'par_mol' to 'mol'.
            new_mol._makeCompatibleWith(par_mol, overwrite=True, verbose=False)

            return new_mol

    def _run_tleap(self, molecule):
        """Run using tLEaP.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
        """

        # Whether the force field is found.
        is_found = False

        # Whether the force field is old.
        is_old = False

        # Search for a compatible force field file.
        if _amber_home is not None:
            ff = _IO.glob("%s/dat/leap/cmd/*.%s" % (_amber_home, self._forcefield))

            # Search the old force fields. First try a specific match.
            if len(ff) == 0:
                ff = _IO.glob("%s/dat/leap/cmd/oldff/leaprc.s" % (_amber_home, self._forcefield))

                # No matches, try globbing all files with matching extension.
                if len(ff) == 0:
                    ff = _IO.glob("%s/oldff/*.%s" % (_cmd_dir, self._forcefield))

                if len(ff) > 0:
                    is_found = True
                    is_old = True

        if not is_found:
            ff = _IO.glob("%s/*.%s" % (_cmd_dir, self._forcefield))

            # Search the old force fields. First try a specific match.
            if len(ff) == 0:
                ff = _IO.glob("%s/oldff/leaprc.%s" % (_cmd_dir, self._forcefield))
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

        # Create a new system and molecule group.
        s = _Sire.System.System("BioSimSpace System")
        m = _Sire.Mol.MoleculeGroup("all")

        # Add the molecule.
        m.add(molecule._molecule)
        s.add(m)

        # Write the system to a PDB file.
        pdb = _Sire.IO.PDB2(s)
        pdb.writeToFile("input.pdb")

        # Write the LEaP input file.
        with open("leap.txt", "w") as f:
            f.write("source %s\n" % ff)
            f.write("mol = loadPdb input.pdb\n")
            f.write("saveAmberParm mol output.top output.crd\n")
            f.write("quit")

        # Generate the tLEaP command.
        command = "%s -f leap.txt" % _tleap_exe

        # Run tLEaP as a subprocess.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # tLEaP doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile("output.top") and _os.path.isfile("output.crd"):
            return ["output.top", "output.crd"]
        else:
            return None

    def _run_pdb2gmx(self, molecule):
        """Run using pdb2gmx.

           Positional arguments:

           molecule -- The molecule to apply the parameterisation protocol to.
        """
        raise NotImplementedError("Parameterisation using pdb2gmx is not yet supported.")

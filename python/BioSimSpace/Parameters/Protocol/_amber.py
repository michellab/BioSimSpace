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
Functionality for handling parameterisation protocols
for AMBER force field models.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

# To override any protocols, just implement a custom "run" method in any
# of the classes.

import Sire as _Sire

from . import _protocol
from ..._SireWrappers import Molecule as _Molecule

import BioSimSpace.IO as _IO

import os as _os
import queue as _queue
import subprocess as _subprocess

__all__ = ["FF03", "FF99", "FF99SB", "FF14SB", "GAFF", "GAFF2"]

class FF03(_protocol.Protocol):
    """A class for handling protocols for the FF03 force field model."""

    def __init__(self, map={}):
        """Constructor.

           Keyword arguments:

           map -- A dictionary that maps system "properties" to their user defined
                  values. This allows the user to refer to properties with their
                  own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff03", map=map)

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF99(_protocol.Protocol):
    """A class for handling protocols for the FF99 force field model."""

    def __init__(self, map={}):
        """Constructor.

           Keyword arguments:

           map -- A dictionary that maps system "properties" to their user defined
                  values. This allows the user to refer to properties with their
                  own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99", map=map)

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF99SB(_protocol.Protocol):
    """A class for handling protocols for the FF99SB force field model."""

    def __init__(self, map={}):
        """Constructor.

           Keyword arguments:

           map -- A dictionary that maps system "properties" to their user defined
                  values. This allows the user to refer to properties with their
                  own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99SB", map=map)

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF14SB(_protocol.Protocol):
    """A class for handling protocols for the FF14SB force field model."""

    def __init__(self):
        """Constructor.

           Keyword arguments:

           map -- A dictionary that maps system "properties" to their user defined
                  values. This allows the user to refer to properties with their
                  own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff14SB", map=map)

        # Set the compatibility flags.
        self._tleap = True

class GAFF(_protocol.Protocol):
    """A class for handling protocols for the GAFF force field model."""

    # A list of supported charge methods.
    _charge_methods = [ "RESP",
                        "CM2",
                        "MUL",
                        "BCC",
                        "ESP",
                        "GAS" ]

    def __init__(self, charge_method="BCC", map={}):
        """Constructor.

           Keyword arguments:

           charge_method -- The method to use when calculating atomic charges.
           map           -- A dictionary that maps system "properties" to their user defined
                            values. This allows the user to refer to properties with their
                            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        if type(charge_method) is not str:
            raise TypeError("'charge_method' must be of type 'str'")

        # Strip whitespace and convert to upper case.
        charge_method = charge_method.replace(" ", "").upper()

        # Check that the charge method is valid.
        if not charge_method in self._charge_methods:
            raise ValueError("Unsupported charge method: '%s'. Supported methods are: %s"
                % (charge_method, self.charge_methods))

        # Set the charge method.
        self._charge_method = charge_method

        # Call the base class constructor.
        super().__init__(forcefield="gaff", map=map)

        self._version = 1

    @classmethod
    def chargeMethods(cls):
        """Return a list of the supported charge methods."""
        return cls._charge_methods

    def chargeMethod(self):
        """Return the chosen charge method."""
        return self._charge_method

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

        # Create a new system and molecule group.
        s = _Sire.System.System("BioSimSpace System")
        m = _Sire.Mol.MoleculeGroup("all")

        # Add the molecule.
        m.add(molecule._getSireMolecule())
        s.add(m)

        # Write the system to a PDB file.
        pdb = _Sire.IO.PDB2(s)
        pdb.writeToFile("antechamber.pdb")

        # Generate the Antechamber command.
        command = ("%s -at %d -i antechamber.pdb -fi pdb " +
                   "-o antechamber.mol2 -fo mol2 -c %s -s 2"
                  ) % (_protocol._antechamber_exe,
                        self._version, self._charge_method.lower())

        with open("README.txt", "w") as f:
            # Write the command to file.
            f.write("# Antechamber was run with the following command:\n")
            f.write("%s\n" % command)

        # Run Antechamber as a subprocess.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Antechamber doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile("antechamber.mol2"):

            # Run parmchk to check for missing parameters.
            command = ("%s -s %d -i antechamber.mol2 -f mol2 " +
                       "-o antechamber.frcmod"
                      ) % (_protocol._parmchk_exe, self._version)

            with open("README.txt", "a") as f:
                # Write the command to file.
                f.write("\n# ParmChk was run with the following command:\n")
                f.write("%s\n" % command)

            # Run parmchk as a subprocess.
            proc = _subprocess.run(command, shell=True,
                stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

            # The frcmod file was created.
            if _os.path.isfile("antechamber.frcmod"):

                # Now call tLEaP using the partially parameterised molecule and the frcmod file.
                # tLEap will run in the same working directory, using the Mol2 file generated by
                # Antechamber.

                # Try to find a force field file.
                if self._version == 1:
                    ff = _protocol._find_force_field("gaff")
                else:
                    ff = _protocol._find_force_field("gaff2")

                # Write the LEaP input file.
                with open("leap.txt", "w") as f:
                    f.write("source %s\n" % ff)
                    f.write("mol = loadMol2 antechamber.mol2\n")
                    f.write("loadAmberParams antechamber.frcmod\n")
                    f.write("saveAmberParm mol leap.top leap.crd\n")
                    f.write("quit")

                # Generate the tLEaP command.
                command = "%s -f leap.txt" % _protocol._tleap_exe

                with open("README.txt", "a") as f:
                    # Write the command to file.
                    f.write("\n# tLEaP was run with the following command:\n")
                    f.write("%s\n" % command)

                # Run tLEaP as a subprocess.
                proc = _subprocess.run(command, shell=True,
                    stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

                # tLEaP doesn't return sensible error codes, so we need to check that
                # the expected output was generated.
                if _os.path.isfile("leap.top") and _os.path.isfile("leap.crd"):
                    # Load the parameterised molecule.
                    par_mol = _Molecule(_IO.readMolecules(["leap.top", "leap.crd"])
                              ._getSireSystem()[_Sire.Mol.MolIdx(0)])

                    # Make the molecule 'mol' compatible with 'par_mol'. This will create
                    # a mapping between atom indices in the two molecules and add all of
                    # the new properties from 'par_mol' to 'mol'.
                    new_mol._makeCompatibleWith(par_mol, map=self._map, overwrite=True, verbose=False)

                else:
                    new_mol = None
            else:
                new_mol = None
        else:
            new_mol = None

        # Change back to the current directory.
        _os.chdir(dir)

        if queue is not None:
            queue.put(new_mol)
        return new_mol

class GAFF2(_protocol.Protocol):
    """A class for handling protocols for the GAFF2 force field model."""

    # Copy the GAFF run method.
    run = GAFF.run

    # Copy the supported charge methods.
    _charge_methods = GAFF._charge_methods
    chargeMethods = GAFF.chargeMethods
    chargeMethod = GAFF.chargeMethod

    def __init__(self, charge_method="BCC", map={}):
        """Constructor.

           Keyword arguments:

           charge_method -- The method to use when calculating atomic charges.
           map           -- A dictionary that maps system "properties" to their user defined
                            values. This allows the user to refer to properties with their
                            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Strip whitespace and convert to upper case.
        charge_method = charge_method.replace(" ", "").upper()

        # Check that the charge method is valid.
        if not charge_method in self._charge_methods:
            raise ValueError("Unsupported charge method: '%s'. Supported methods are: %s"
                % (charge_method, self._charge_methods))

        # Set the charge method.
        self._charge_method = charge_method

        # Call the base class constructor.
        super().__init__(forcefield="gaff2")

        self._version = 2

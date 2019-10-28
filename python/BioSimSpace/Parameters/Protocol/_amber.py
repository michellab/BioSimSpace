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
Functionality for handling parameterisation protocols
for AMBER force field models.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["FF03", "FF99", "FF99SB", "FF99SBILDN", "FF14SB", "GAFF", "GAFF2"]

# To override any protocols, just implement a custom "run" method in any
# of the classes.

import os as _os
import queue as _queue
import subprocess as _subprocess
import warnings as _warnings

from Sire import IO as _SireIO
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace import _isVerbose
from BioSimSpace import IO as _IO
from BioSimSpace._Exceptions import ParameterisationError as _ParameterisationError
from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace.Parameters._utils import formalCharge as _formalCharge
from BioSimSpace.Types import Charge as _Charge

from . import _protocol

class FF03(_protocol.Protocol):
    """A class for handling protocols for the FF03 force field model."""

    def __init__(self, property_map={}):
        """Constructor.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff03", property_map=property_map)

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF99(_protocol.Protocol):
    """A class for handling protocols for the FF99 force field model."""

    def __init__(self, property_map={}):
        """Constructor.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99", property_map=property_map)

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = True

class FF99SB(_protocol.Protocol):
    """A class for handling protocols for the FF99SB force field model."""

    def __init__(self, property_map={}):
        """Constructor.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99SB", property_map=property_map)

        # Set the compatibility flags.
        self._tleap = True

class FF99SBILDN(_protocol.Protocol):
    """A class for handling protocols for the FF99SBILDN force field model."""

    def __init__(self, property_map={}):
        """Constructor.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99SBildn", property_map=property_map)

        # Set the compatibility flags.
        self._tleap = True
        self._pdb2gmx = False

class FF14SB(_protocol.Protocol):
    """A class for handling protocols for the FF14SB force field model."""

    def __init__(self, property_map={}):
        """Constructor.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff14SB", property_map=property_map)

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

    def __init__(self, charge_method="BCC", net_charge=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           charge_method : str
               The method to use when calculating atomic charges:
               "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

           net_charge: int
               The net charge on the molecule.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
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

        if net_charge is not None:
            # Get the magnitude of the charge.
            if type(net_charge) is _Charge:
                net_charge = net_charge.magnitude()

            if type(net_charge) is float:
                if net_charge % 1 != 0:
                    raise ValueError("'net_charge' must be integer valued.")

            # Try to convert to int.
            try:
                net_charge = int(net_charge)
            except:
                raise TypeError("'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'")

        # Set the net molecular charge.
        self._net_charge = net_charge

        # Set the charge method.
        self._charge_method = charge_method

        # Call the base class constructor.
        super().__init__(forcefield="gaff", property_map=property_map)

        self._version = 1

    @classmethod
    def chargeMethods(cls):
        """Return a list of the supported charge methods.

           Returns
           -------

           charge_methods : [str]
               A list of supported charge methods.
        """
        return cls._charge_methods

    def chargeMethod(self):
        """Return the chosen charge method.

           Returns
           -------

           charge_method : str
               The chosen charge method.
        """
        return self._charge_method

    def netCharge(self):
        """Return the net molecular charge.

           Returns
           -------

           net_charge : int
               The net molecular charge.
        """
        return self._net_charge

    def run(self, molecule, work_dir=None, queue=None):
        """Run the parameterisation protocol.

           Parameters
           ----------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule to apply the parameterisation protocol to.

           work_dir : str
               The working directory.

           queue : queue.Queue
               The thread queue is which this method has been run.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The parameterised molecule.
        """

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

        if type(work_dir) is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        if type(queue) is not None and type(queue) is not _queue.Queue:
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Create the file prefix.
        prefix = work_dir + "/"

        # Create a copy of the molecule.
        new_mol = molecule.copy()

        # Use the net molecular charge passed as an option.
        if self._net_charge is not None:
            charge = self._net_charge
        else:
            # The user will likely have passed a bare PDB or Mol2 file.
            # Antechamber expects the molecule to be uncharged, or integer
            # charged (where the charge, or number of electrons, is passed with
            # the -nc flag).

            # Get the total charge on the molecule.
            if "charge" in self._property_map:
                _property_map = { "charge": self._property_map["charge"] }
                prop = self._property_map["charge"]
            else:
                _property_map = { "charge": "charge" }
                prop = "charge"

            # The molecule has a charge property.
            if new_mol._getSireObject().hasProperty(prop):
                charge = new_mol.charge(property_map=_property_map).magnitude()

                # Charge is non-integer, try to fix it.
                if abs(round(charge) - charge) > 0:
                    new_mol._fixCharge(property_map=_property_map)
                    charge = round(charge)
            else:
                charge = None

            # Only try "formal_charge" when "charge" is missing. Unlikely to have
            # both if this is a bare molecule, but the user could be re-parameterising
            # an existing molecule.
            if charge is None:
                # Get the total formal charge on the molecule.
                if "formal_charge" in self._property_map:
                    _property_map = { "charge": self._property_map["formal_charge"] }
                    prop = self._property_map["charge"]
                else:
                    _property_map = { "charge": "formal_charge" }
                    prop = "formal_charge"

                if new_mol._getSireObject().hasProperty(prop):
                    charge = new_mol.charge(property_map=_property_map).magnitude()

                    # Compute the formal charge ourselves to check that it is consistent.
                    formal_charge = _formalCharge(molecule).magnitude()

                    if charge != formal_charge:
                        _warnings.warn("The formal charge on the molecule is %d "
                                       "but we estimate it to be %d" % (charge, formal_charge))
                else:
                    msg = ("The molecule has no 'charge' or 'formal_charge' information, and "
                           "no 'net_charge' option has been passed. You can use the "
                           "'BioSimSpace.Parameters.formalCharge' function to compute the "
                           "formal charge")
                    raise _ParameterisationError(msg)

        # Create a new system and molecule group.
        s = _SireSystem.System("BioSimSpace System")
        m = _SireMol.MoleculeGroup("all")

        # Add the molecule.
        m.add(new_mol._getSireObject())
        s.add(m)

        # Write the system to a PDB file.
        try:
            pdb = _SireIO.PDB2(s)
            pdb.writeToFile(prefix + "antechamber.pdb")
        except Exception as e:
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Generate the Antechamber command.
        command = ("%s -at %d -i antechamber.pdb -fi pdb " +
                "-o antechamber.mol2 -fo mol2 -c %s -s 2 -nc %d"
                ) % (_protocol._antechamber_exe, self._version,
                        self._charge_method.lower(), charge)

        with open(prefix + "README.txt", "w") as file:
            # Write the command to file.
            file.write("# Antechamber was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open(prefix + "antechamber.out", "w")
        stderr = open(prefix + "antechamber.err", "w")

        # Run Antechamber as a subprocess.
        proc = _subprocess.run(command, cwd=work_dir, shell=True, stdout=stdout, stderr=stderr)
        stdout.close()
        stderr.close()

        # Antechamber doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile(prefix + "antechamber.mol2"):

            # Run parmchk to check for missing parameters.
            command = ("%s -s %d -i antechamber.mol2 -f mol2 " +
                    "-o antechamber.frcmod"
                    ) % (_protocol._parmchk_exe, self._version)

            with open(prefix + "README.txt", "a") as file:
                # Write the command to file.
                file.write("\n# ParmChk was run with the following command:\n")
                file.write("%s\n" % command)

            # Create files for stdout/stderr.
            stdout = open(prefix + "parmchk.out", "w")
            stderr = open(prefix + "parmchk.err", "w")

            # Run parmchk as a subprocess.
            proc = _subprocess.run(command, cwd=work_dir, shell=True, stdout=stdout, stderr=stderr)
            stdout.close()
            stderr.close()

            # The frcmod file was created.
            if _os.path.isfile(prefix + "antechamber.frcmod"):

                # Now call tLEaP using the partially parameterised molecule and the frcmod file.
                # tLEap will run in the same working directory, using the Mol2 file generated by
                # Antechamber.

                # Try to find a force field file.
                if self._version == 1:
                    ff = _protocol._find_force_field("gaff")
                else:
                    ff = _protocol._find_force_field("gaff2")

                # Write the LEaP input file.
                with open(prefix + "leap.txt", "w") as file:
                    file.write("source %s\n" % ff)
                    file.write("mol = loadMol2 antechamber.mol2\n")
                    file.write("loadAmberParams antechamber.frcmod\n")
                    file.write("saveAmberParm mol leap.top leap.crd\n")
                    file.write("quit")

                # Generate the tLEaP command.
                command = "%s -f leap.txt" % _protocol._tleap_exe

                with open(prefix + "README.txt", "a") as file:
                    # Write the command to file.
                    file.write("\n# tLEaP was run with the following command:\n")
                    file.write("%s\n" % command)

                # Create files for stdout/stderr.
                stdout = open(prefix + "tleap.out", "w")
                stderr = open(prefix + "tleap.err", "w")

                # Run tLEaP as a subprocess.
                proc = _subprocess.run(command, cwd=work_dir, shell=True, stdout=stdout, stderr=stderr)
                stdout.close()
                stderr.close()

                # tLEaP doesn't return sensible error codes, so we need to check that
                # the expected output was generated.
                if _os.path.isfile(prefix + "leap.top") and _os.path.isfile(prefix + "leap.crd"):
                    # Load the parameterised molecule.
                    try:
                        par_mol = _Molecule(_IO.readMolecules([prefix + "leap.top", prefix + "leap.crd"])
                                ._getSireObject()[_SireMol.MolIdx(0)])
                    except Exception as e:
                        msg = "Failed to read molecule from: 'leap.top', 'leap.crd'"
                        if _isVerbose():
                            raise IOError(msg) from e
                        else:
                            raise IOError(msg) from None

                    # Make the molecule 'mol' compatible with 'par_mol'. This will create
                    # a mapping between atom indices in the two molecules and add all of
                    # the new properties from 'par_mol' to 'mol'.
                    new_mol._makeCompatibleWith(par_mol, property_map=self._property_map, overwrite=True, verbose=False)

                    # Record the forcefield used to parameterise the molecule.
                    new_mol._forcefield = ff

                else:
                    raise _ParameterisationError("tLEaP failed!")
            else:
                raise _ParameterisationError("Parmchk failed!")
        else:
            raise _ParameterisationError("Antechamber failed!")

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

    def __init__(self, charge_method="BCC", net_charge=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           charge_method : str
               The method to use when calculating atomic charges:
               "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

           net_charge: int
               The net charge on the molecule.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Strip whitespace and convert to upper case.
        charge_method = charge_method.replace(" ", "").upper()

        # Check that the charge method is valid.
        if not charge_method in self._charge_methods:
            raise ValueError("Unsupported charge method: '%s'. Supported methods are: %s"
                % (charge_method, self._charge_methods))

        if net_charge is not None:
            # Get the magnitude of the charge.
            if type(net_charge) is _Charge:
                net_charge = net_charge.magnitude()

            if type(net_charge) is float:
                if net_charge % 1 != 0:
                    raise ValueError("'net_charge' must be integer valued.")

            # Try to convert to int.
            try:
                net_charge = int(net_charge)
            except:
                raise TypeError("'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'")

        # Set the net molecular charge.
        self._net_charge = net_charge

        # Set the charge method.
        self._charge_method = charge_method

        # Call the base class constructor.
        super().__init__(forcefield="gaff2")

        self._version = 2

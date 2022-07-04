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
__email__ = "lester.hedges@gmail.com"

__all__ = ["FF03", "FF99", "FF99SB", "FF99SBILDN", "FF14SB", "GAFF", "GAFF2"]

# To override any protocols, just implement a custom "run" method in any
# of the classes.

import os as _os
import queue as _queue
import shlex as _shlex
import subprocess as _subprocess
import warnings as _warnings

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from ... import _isVerbose
from ... import IO as _IO
from ..._Exceptions import ParameterisationError as _ParameterisationError
from ..._Exceptions import ThirdPartyError as _ThirdPartyError
from ..._SireWrappers import Molecule as _Molecule
from ...Parameters._utils import formalCharge as _formalCharge
from ...Types import Charge as _Charge
from ... import _Utils

from . import _protocol

class FF03(_protocol.Protocol):
    """A class for handling protocols for the FF03 force field model."""

    def __init__(self, water_model=None, leap_commands=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           water_model : str
               The water model used to parameterise any structural ions. This
               will be ignored when it is not supported by the chosen force field.
               Run 'BioSimSpace.Solvent.waterModels()' to see the supported
               water models.

           leap_commands : [str]
               An optional list of extra commands for the LEaP program. These
               will be added after any default commands and can be used to, e.g.,
               load additional parameter files. When this option is set, we can no
               longer fall back on GROMACS's pdb2gmx.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff03", property_map=property_map)

        # Validate the water model.
        if water_model is not None and not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'")
        else:
            self._water_model = water_model

        # Validate the additional leap commands.
        if leap_commands is not None:
            if not isinstance(leap_commands, (list, tuple)):
                raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in leap_commands):
                    raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
        self._leap_commands = leap_commands

        # Set the compatibility flags.
        self._tleap = True
        if self._leap_commands is None:
            self._pdb2gmx = True

class FF99(_protocol.Protocol):
    """A class for handling protocols for the FF99 force field model."""

    def __init__(self, water_model=None, leap_commands=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           water_model : str
               The water model used to parameterise any structural ions. This
               will be ignored when it is not supported by the chosen force field.
               Run 'BioSimSpace.Solvent.waterModels()' to see the supported
               water models.

           leap_commands : [str]
               An optional list of extra commands for the LEaP program. These
               will be added after any default commands and can be used to, e.g.,
               load additional parameter files. When this option is set, we can no
               longer fall back on GROMACS's pdb2gmx.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99", property_map=property_map)

        # Validate the water model.
        if water_model is not None and not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'")
        else:
            self._water_model = water_model

        # Validate the additional leap commands.
        if leap_commands is not None:
            if not isinstance(leap_commands, (list, tuple)):
                raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in leap_commands):
                    raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
        self._leap_commands = leap_commands

        # Set the compatibility flags.
        self._tleap = True
        if self._leap_commands is None:
            self._pdb2gmx = True

class FF99SB(_protocol.Protocol):
    """A class for handling protocols for the FF99SB force field model."""

    def __init__(self, water_model=None, leap_commands=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           water_model : str
               The water model used to parameterise any structural ions. This
               will be ignored when it is not supported by the chosen force field.
               Run 'BioSimSpace.Solvent.waterModels()' to see the supported
               water models.

           leap_commands : [str]
               An optional list of extra commands for the LEaP program. These
               will be added after any default commands and can be used to, e.g.,
               load additional parameter files. When this option is set, we can no
               longer fall back on GROMACS's pdb2gmx.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99SB", property_map=property_map)

        # Validate the water model.
        if water_model is not None and not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'")
        else:
            self._water_model = water_model

        # Validate the additional leap commands.
        if leap_commands is not None:
            if not isinstance(leap_commands, (list, tuple)):
                raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in leap_commands):
                    raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
        self._leap_commands = leap_commands

        # Set the compatibility flags.
        self._tleap = True

class FF99SBILDN(_protocol.Protocol):
    """A class for handling protocols for the FF99SBILDN force field model."""

    def __init__(self, water_model=None, leap_commands=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           water_model : str
               The water model used to parameterise any structural ions. This
               will be ignored when it is not supported by the chosen force field.
               Run 'BioSimSpace.Solvent.waterModels()' to see the supported
               water models.

           leap_commands : [str]
               An optional list of extra commands for the LEaP program. These
               will be added after any default commands and can be used to, e.g.,
               load additional parameter files. When this option is set, we can no
               longer fall back on GROMACS's pdb2gmx.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff99SBildn", property_map=property_map)

        # Validate the water model.
        if water_model is not None and not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'")
        else:
            self._water_model = water_model

        # Validate the additional leap commands.
        if leap_commands is not None:
            if not isinstance(leap_commands, (list, tuple)):
                raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in leap_commands):
                    raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
        self._leap_commands = leap_commands

        # Set the compatibility flags.
        self._tleap = True

class FF14SB(_protocol.Protocol):
    """A class for handling protocols for the FF14SB force field model."""

    def __init__(self, water_model=None, leap_commands=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           water_model : str
               The water model used to parameterise any structural ions. This
               will be ignored when it is not supported by the chosen force field.
               Run 'BioSimSpace.Solvent.waterModels()' to see the supported
               water models.

           leap_commands : [str]
               An optional list of extra commands for the LEaP program. These
               will be added after any default commands and can be used to, e.g.,
               load additional parameter files. When this option is set, we can no
               longer fall back on GROMACS's pdb2gmx.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield="ff14SB", property_map=property_map)

        # Validate the water model.
        if water_model is not None and not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'")
        else:
            self._water_model = water_model

        # Validate the additional leap commands.
        if leap_commands is not None:
            if not isinstance(leap_commands, (list, tuple)):
                raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in leap_commands):
                    raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
        self._leap_commands = leap_commands

        # Set the compatibility flags.
        self._tleap = True

class GAFF(_protocol.Protocol):
    """A class for handling protocols for the GAFF force field model."""

    # A list of supported charge methods.
    _charge_methods = ["RESP",
                       "CM2",
                       "MUL",
                       "BCC",
                       "ESP",
                       "GAS"]

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

        if not isinstance(charge_method, str):
            raise TypeError("'charge_method' must be of type 'str'")

        # Strip whitespace and convert to upper case.
        charge_method = charge_method.replace(" ", "").upper()

        # Check that the charge method is valid.
        if not charge_method in self._charge_methods:
            raise ValueError("Unsupported charge method: '%s'. Supported methods are: %s"
                % (charge_method, self._charge_methods))

        if net_charge is not None:
            # Get the value of the charge.
            if isinstance(net_charge, _Charge):
                net_charge = net_charge.value()

            if isinstance(net_charge, float):
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

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
               The molecule to parameterise, either as a Molecule object or SMILES
               string.

           work_dir : str
               The working directory.

           queue : queue.Queue
               The thread queue is which this method has been run.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The parameterised molecule.
        """

        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'")

        if work_dir is not None and not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

        if queue is not None and not isinstance(queue, _queue.Queue):
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Create the file prefix.
        prefix = work_dir + "/"

        # Convert SMILES to a molecule.
        if isinstance(molecule, str):
            is_smiles = True
            try:
                new_mol = self._smiles_to_molecule(molecule, work_dir)
            except Exception as e:
                msg = "Unable to convert SMILES to Molecule using Open Force Field."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None
        else:
            is_smiles = False
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
                charge = new_mol.charge(property_map=_property_map).value()

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
                    charge = new_mol.charge(property_map=_property_map).value()

                    # Compute the formal charge ourselves to check that it is consistent.
                    formal_charge = _formalCharge(new_mol).value()

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
                msg += ": " + getattr(e, "message", repr(e))
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
        proc = _subprocess.run(_Utils.command_split(command), cwd=work_dir,
            shell=False, stdout=stdout, stderr=stderr)
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
            proc = _subprocess.run(_Utils.command_split(command), cwd=work_dir,
                shell=False, stdout=stdout, stderr=stderr)
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
                stdout = open(prefix + "leap.out", "w")
                stderr = open(prefix + "leap.err", "w")

                # Run tLEaP as a subprocess.
                proc = _subprocess.run(_Utils.command_split(command), cwd=work_dir,
                    shell=False, stdout=stdout, stderr=stderr)
                stdout.close()
                stderr.close()

                # tLEaP doesn't return sensible error codes, so we need to check that
                # the expected output was generated.
                if _os.path.isfile(prefix + "leap.top") and _os.path.isfile(prefix + "leap.crd"):
                    # Check the output of tLEaP for missing atoms.
                    if self._has_missing_atoms(prefix + "leap.out"):
                        raise _ParameterisationError("tLEaP added missing atoms. The topology is now "
                                                     "inconsistent with the original molecule. Please "
                                                     "make sure that your initial molecule has a "
                                                     "complete topology.")

                    # Load the parameterised molecule. (This could be a system of molecules.)
                    try:
                        par_mol = _IO.readMolecules([prefix + "leap.top", prefix + "leap.crd"])
                        # Extract single molecules.
                        if par_mol.nMolecules() == 1:
                            par_mol = par_mol.getMolecules()[0]
                    except Exception as e:
                        msg = "Failed to read molecule from: 'leap.top', 'leap.crd'"
                        if _isVerbose():
                            msg += ": " + getattr(e, "message", repr(e))
                            raise IOError(msg) from e
                        else:
                            raise IOError(msg) from None

                    # Make the molecule 'mol' compatible with 'par_mol'. This will create
                    # a mapping between atom indices in the two molecules and add all of
                    # the new properties from 'par_mol' to 'mol'.
                    if is_smiles:
                        new_mol = par_mol

                        # We'll now add MolName and ResName info to the molecule, since
                        # this will be missing.

                        # Rename the molecule with the original SMILES string.
                        edit_mol = new_mol._sire_object.edit()
                        edit_mol = edit_mol.rename(molecule).molecule()

                        # Rename the residue LIG.
                        resname = _SireMol.ResName("LIG")
                        edit_mol = edit_mol.residue(_SireMol.ResIdx(0)).rename(resname).molecule()

                        # Commit the changes.
                        new_mol._sire_object = edit_mol.commit()

                    else:
                        new_mol.makeCompatibleWith(par_mol, property_map=self._property_map, overwrite=True, verbose=False)

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
            # Get the value of the charge.
            if isinstance(net_charge, _Charge):
                net_charge = net_charge.value()

            if isinstance(net_charge, float):
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

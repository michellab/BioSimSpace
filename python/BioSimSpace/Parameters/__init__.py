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
Functionality for parameterising molecules.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule

from .._Utils import Antechamber as _Antechamber
from .._Utils import Tleap as _Tleap

import multiprocessing as _multiprocessing

def ff99(molecules, work_dir=None, verbose=False):
    """Parameterise using the ff99 force field.

       Positional arguments:

       molecules -- The molecules to parameterise.

       Keyword arguments:

       work_dir  -- The working directory for external processes.
       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "ff99", work_dir, verbose)

def ff99SB(molecules, work_dir=None, verbose=False):
    """Parameterise using the ff99SB force field.

       Positional arguments:

       molecules -- The molecules to parameterise.

       Keyword arguments:

       work_dir  -- The working directory for external processes.
       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "ff99SB", work_dir, verbose)

def ff03(molecules, work_dir=None, verbose=False):
    """Parameterise using the ff03 force field.

       Positional arguments:

       molecules -- The molecules to parameterise.

       Keyword arguments:

       work_dir  -- The working directory for external processes.
       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "ff03", work_dir, verbose)

def ff14(molecules, work_dir=None, verbose=False):
    """Parameterise using the ff14 force field.

       Positional arguments:

       molecules -- The molecules to parameterise.

       Keyword arguments:

       work_dir  -- The working directory for external processes.
       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "ff14", work_dir, verbose)

def ff14SB(molecules, work_dir=None, verbose=False):
    """Parameterise using the ff14SB force field.

       Positional arguments:

       molecules -- The molecules to parameterise.
       work_dir  -- The working directory for external processes.

       Keyword arguments:

       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "ff14SB", work_dir, verbose)

def gaff(molecules, work_dir=None, verbose=False):
    """Parameterise using the gaff force field.

       Positional arguments:

       molecules -- The molecules to parameterise.

       Keyword arguments:

       work_dir  -- The working directory for external processes.
       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "gaff", work_dir, verbose)

def gaff2(molecules, work_dir=None, verbose=False):
    """Parameterise using the gaff force field.

       Positional arguments:

       molecules -- The molecules to parameterise.

       Keyword arguments:

       work_dir  -- The working directory for external processes.
       verbose   -- Whether to report stdout/stderr of external processes.
    """
    return _parameterise(molecules, "gaff2", work_dir, verbose)

def _parameterise(molecules, forcefield, work_dir=None, verbose=False):
    """Internal function to parameterise a set of molecules using a given force field.

       Positional arguments:

       molecules  -- The molecules to parameterise.
       forcefield -- The force field to use.

       Keyword arguments:

       work_dir   -- The working directory for external processes.
       verbose    -- Whether to report stdout/stderr of external processes.
    """

    if work_dir is not None and type(work_dir) is not str:
        raise TypeError("'work_dir' must be of type 'str'")

    if type(verbose) is not bool:
        raise TypeError("'verbose' must be of type 'bool'")

    # A list of Molecule objects.
    if type(molecules) is list:
        if all(isinstance(x, _Molecule) for x in molecules):
            is_list = True
            is_system = False

    # A single Molecule object.
    elif type(molecules) is _Molecule:
        molecules = [molecules]
        is_list = False
        is_system = False

    # Unsupported.
    else:
        raise TypeError("'molecules' must be of type 'BioSimSpace.Molecule', "
            + "or a list of objects of this type.")

    # Create a list of new molecules.
    new_mols = []

    # Loop over all of the molecules.
    for mol in molecules:
        new_mols.append(_parameterise_molecule(mol, forcefield,
            work_dir=work_dir, verbose=verbose))

    if is_list:
        return new_mols
    else:
        return new_mols[0]

# Create a list of the force field names.
# This needs to come after all of the force field functions.
_forcefields = []
for _var in dir():
    if _var[0] != "_":
        _forcefields.append(_var)

def forceFields():
    "Print a list of the supported force fields"
    print(", ".join(_forcefields))

def _parameterise_molecule(molecule, forcefield, work_dir=None, verbose=False):
    """Internal function to parameterise a single molecule using a given force field.

       Positional arguments:

       molecule   -- The molecules to parameterise.
       forcefield -- The force field to use.

       Keyword arguments:

       work_dir   -- The working directory for external processes.
       verbose    -- Whether to report stdout/stderr of external processes.
    """

    # A single Molecule object.
    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    # Create a new molecule using a deep copy of the internal Sire Molecule.
    new_mol = _Molecule(molecule._molecule.__deepcopy__())

    # Parameterise the molecule.
    if forcefield == "gaff" or forcefield == "gaff2":
        par_mol = _Molecule(_Antechamber.parameterise(molecule._getSireMolecule(),
            forcefield, work_dir=work_dir, verbose=verbose))
    else:
        par_mol = _Molecule(_Tleap.parameterise(molecule._getSireMolecule(),
            forcefield, work_dir=work_dir, verbose=verbose))

    # Make the molecule 'mol' compatible with 'par_mol'. This will create
    # a mapping between atom indices in the two molecules and add all of
    # the new properties from 'par_mol' to 'mol'.
    new_mol._makeCompatibleWith(par_mol, overwrite=True, verbose=verbose)

    # Return the updated molecule.
    return new_mol

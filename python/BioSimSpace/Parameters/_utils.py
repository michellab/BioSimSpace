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
Utility functions.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["formalCharge"]

import tempfile as _tempfile

from BioSimSpace import _is_notebook
from BioSimSpace import IO as _IO
from BioSimSpace import _Utils as _Utils
from BioSimSpace.Units.Charge import electron_charge as _electron_charge
from BioSimSpace._SireWrappers import Molecule as _Molecule

def formalCharge(molecule):
    """Compute the formal charge on a molecule. This function requires that
       the molecule has explicit hydrogen atoms.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           A molecule object.

       Returns
       -------

       formal_charge : :class:`Charge <BioSimSpace.Types.Charge>`
           The total formal charge on the molecule.
    """

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    from rdkit import Chem as _Chem

    # Create a temporary working directory.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Zero the total formal charge.
    formal_charge = 0

    # Stdout/stderr redirection doesn't work from within Jupyter.
    if _is_notebook:
        # Run in the working directory.
        with _Utils.cd(work_dir):

            # Save the molecule to a PDB file.
            _IO.saveMolecules("tmp", molecule, "PDB")

            # Read the ligand PDB into an RDKit molecule.
            mol = _Chem.MolFromPDBFile("tmp.pdb")

            # Compute the formal charge.
            formal_charge = _Chem.rdmolops.GetFormalCharge(mol)

    else:
        # Run in the working directory and redirect stderr from RDKit.
        with _Utils.cd(work_dir), _Utils.stderr_redirected():

            # Save the molecule to a PDB file.
            _IO.saveMolecules("tmp", molecule, "PDB")

            # Read the ligand PDB into an RDKit molecule.
            mol = _Chem.MolFromPDBFile("tmp.pdb")

            # Compute the formal charge.
            formal_charge = _Chem.rdmolops.GetFormalCharge(mol)

    return formal_charge * _electron_charge

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

import BioSimSpace.IO as _IO
import BioSimSpace._Utils as _Utils

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

    # Create a copy of the molecule.
    mol = molecule.copy()

    # Delete any existing formal charge information.
    try:
        mol._sire_object = mol.edit().removeProperty("formal_charge").molecule().commit()
    except:
        pass

    import openbabel as ob

    # Create a temporary working directory.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Zero the total formal charge.
    formal_charge = 0

    # Run in the working directory.
    with _Utils.cd(work_dir):

        # Save the molecule to a PDB file.
        _IO.saveMolecules("tmp", molecule, "PDB")

        # Read the ligand PDB into an OpenBabel molecule.
        conv = ob.OBConversion()
        mol = ob.OBMol()
        conv.ReadFile(mol, "tmp.pdb")

        # Adapted from FESetup.

        phospho = '[#15D4](~[OD1])(~[OD1])(~[OD2])~[OD2]'
        sulfonyl = '[#16](~[OD1])(~[OD1])'

        # This compensates for missing charge determination by Openbabel, so
        # will fail to detect charges in cases not covered here. Relies on
        # explicit hydrogens!
        CHARGE_SMARTS = (
            ('[#16;$([#16D4](~[OD1])(~[OD1])~[OD1])]', -1),             # sulfates, Sac
            ('[#16;$([#16D3](~[OD1])~[OD1])]', -1),                     # sulfonates, Sac
            ('[#15;$([#15D4](~[OD1])(~[OD1])~[OD1])]', -2),             # phosphates, Pac
            ('[#15;$([#15D3](~[OD1])~[OD1])]', -2),                     # phosphonates, Pac
            ('[#15;$(%s)]' % phospho, -1),                              # polyphosphates or phosphate diesters
            ('[CD3](~[NH2])(~[NH2])~[NH]', +1),                         # guanidinium
            ('O=C~[CD3,ND2]~C=O', -1),                                  # beta-dicarbonyls, imides
            ('%s~[ND2]~C(=O)~[NH]' % sulfonyl, -1),                     # sulfonylureas
            ('%s~[ND2]~C(=O)~[#6]' % sulfonyl, -1),                     # (sulfon)imides
            ('%s~[ND2]~c' % sulfonyl, -1),                              # N-aryl sulfonamides
            ('{0}~[ND2]~{0}'.format(sulfonyl), -1),                     # sulfonimides 2
            ('[#6][ND3]([#6])~[#6H]~[ND3]([#6])[#6]', +1),              # formamidiniums
            ('C(~[OD1])~[OD1]', -1),                                    # carbonates, Cac in atomtyp.txt
            ('[#6][#16D1]', -1),                                        # thiolates
            ('[#6][#8D1]', -1),                                         # alkoxides (alcoholates)
            ('[#7D4]', +1),                                             # atomtyp.txt checks for [#7X4]
            ('{0}~[ND3](~{0})~[CD3](~{0})~{0}'.format('[#6,#1]'), +1),  # imines
            )

        # Initialise the smarts pattern object.
        sm = ob.OBSmartsPattern()

        # Loop over all smarts strings / formal charge pairs.
        for smarts, fchg in CHARGE_SMARTS:
            sm.Init(smarts)

            # Does the molecule match this string.
            # If so, update the formal charge.
            if sm.Match(mol):
                m = list(sm.GetUMapList() )
                formal_charge += fchg * len(m)

    return formal_charge * _electron_charge

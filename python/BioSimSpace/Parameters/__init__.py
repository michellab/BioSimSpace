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
.. currentmodule:: BioSimSpace.Parameters

Functions
=========

.. autosummary::
    :toctree: generated/

    parameterise
    ff99
    ff99SB
    ff99SBildn
    ff14SB
    gaff
    gaff2
    forceFields
    formalCharge

Examples
========

Print the list of supported force fields.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Parameters.forceFields())

Parameterise a molecule using GAFF.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Initialise the parameterisation process. This will run the parameterisation
   # in the background. This is useful when you are working interactively and
   # wish to continue doing other things while the parameterisation runs.
   # (Parameterisation can be slow.)
   process = BSS.Parameters.gaff(molecule)

   # Get the parameterised molecule. This will now block until the
   # parameterisation is finished.
   molecule = process.getMolecule()

The same as above, but immediately getting the parameterised molecule from
the process object, i.e. not saving it as an intermediate variable.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Initialise the parameterisation process and block until the molecule is
   # ready to be returned.
   molecule = BSS.Parameters.gaff(molecule).getMolecule()

The same as above, but instead passing "GAFF" as an argument to the
:class:`parameterise <BioSimSpace.Parameters.parameterise>` function. This
function should be used in any interoperable workflow
:class:`Node <BioSimSpace.Gateway.Node>` where the force field is specified
as an input requirement by the user.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Initialise the parameterisation process and block until the molecule is
   # ready to be returned.
   molecule = BSS.Parameters.parameterise(molecule, "gaff").getMolecule()

Parameterise a molecule using GAFF, passing the net formal charge computed
by BioSimSpace. This is useful when your input PDB file has missing or
incorrect formal charge information.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Compute the net formal charge on the molecule.
   formal_charge = BSS.Parameters.formalCharge(molecule)

   # Initialise the parameterisation process. This will run the parameterisation
   # in the background. This is useful when you are working interactively and
   # wish to continue doing other things while the parameterisation runs.
   # (Parameterisation can be slow.)
   process = BSS.Parameters.gaff(molecule, net_charge=formal_charge)

   # Get the parameterised molecule. This will now block until the
   # parameterisation is finished.
   molecule = process.getMolecule()
"""

from ._parameters import *
from ._utils import *

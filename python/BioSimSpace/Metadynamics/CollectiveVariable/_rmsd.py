######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
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
Functionality for a root-mean-square deviation collective variable.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["RMSD"]

from math import ceil as _ceil

from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace._SireWrappers import System as _System

from ._collective_variable import CollectiveVariable as _CollectiveVariable
from .._bound import Bound as _Bound
from .._grid import Grid as _Grid

class RMSD(_CollectiveVariable):
    """A class for a root-mean-square deviation (RMSD) collective variable."""

    def __init__(self, system, reference, reference_index=None, rmsd_indices=None,
            hill_width=0.1, lower_bound=None, upper_bound=None, grid=None,
            alignment_type="optimal", pbc=True):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system of interest.

           reference : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The reference molecule, against which the RMSD will be measured.
               This molecule should match with a single molecule from the
               system, i.e. contain the same residues as the matching molecule
               in the same order.

           reference_index : int
               The index of the molecule in the system that matches the
               reference. If none, then we will assume that the molecule with
               the closest number of residues is the match. If not all atoms
               from the reference are matched in the system, then an exception
               will be thrown.

           rmsd_indices : None
               The indices of the atoms within the reference that should be
               used when calculating the RMSD. If None, then all atoms from
               the reference will be used. Atoms from the reference that are
               not in rmsd_indices will instead be used for alignment.

           hill_width : float
               The width of the Gaussian hill used to sample this variable.

           lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               A lower bound on the value of the collective variable.

           upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               An upper bound on the value of the collective variable.

           grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
               The grid on which the collective variable will be sampled.
               This can help speed up long metadynamics simulations where
               the number of Gaussian kernels can become prohibitive.

           alignment_type : str
               The mannier in which RMSD alignment is performed. Options are
               "optimal" or "simple".

           pbc : bool
               Whether to use periodic boundary conditions when computing the
               collective variable.
        """

        # Call the base class constructor.
        super().__init__()

        # Validate input.

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'.")

        if type(reference) is not _Molecule:
            raise TypeError("'reference' must be of type 'BioSimSpace._SireWrappers.Molecule'.")
        self._reference = reference.copy()

        # Extract the corresponding molecule from the system by index.
        if reference_index is not None:
            if type(reference_index) is not int:
                raise TypeError("'reference_index' must be of type 'int'.")
            molecule = system[reference_index]

        # Match to the molecule with the closest number of residues.
        else:
            # Store the number of residues in the target.
            num_res_ref = reference.nResidues()

            # Create a list to hold the number of residues for each molecule in the system.
            num_res_system = []
            for mol in system:
                num_res_system.append(mol.nResidues())

            # Find the index of the molecule with the closest number of residues.
            idx = min(enumerate(num_res_system), key=lambda x: abs(x[1]-num_res_ref))[0]

            # Extract the molecule from the system.
            molecule = system[idx]

        if rmsd_indices is not None:
            # Convert tuple to list.
            if type(rmsd_indices) is tuple:
                rmsd_indices = list(rmsd_indices)

            for idx in rmsd_indices:
                if type(idx) is not int:
                    raise TypeError("'rmsd_indices' must be a list of 'int' types.")
                if idx < 0:
                    idx += reference.nAtoms()
                if idx < 0 or idx >= reference.nAtoms():
                    raise ValueError("'rmsd_indices' is out of range: "
                                    f"{-reference.nAtoms()}:{reference.nAtoms()-1}")

            # Convert to AtomIdx.
            rmsd_indices = [_SireMol.AtomIdx(x) for x in rmsd_indices]

        else:
            # All indices are used for RMSD by default.
            rmsd_indices = [_SireMol.AtomIdx(x) for x in range(0, reference.nAtoms())]

        # Create the matcher.
        matcher = _SireMol.ResNumAtomNameMatcher()

        # Match the atoms by residue number and name.
        matches = matcher.match(molecule._sire_object, reference._sire_object)

        # We need to match all of the atoms in the reference.
        if len(matches) < reference.nAtoms():
            raise _IncompatibleError("Didn't match all of the atoms in the reference molecule: "
                                    f"Found {len(matches)}, expected {reference.nAtoms()}.")

        # Get the indices of the matching atoms.
        idx_matches = matches.keys()

        # Extract the molecule of interest and make it editable.
        edit_mol = molecule._sire_object.edit()

        # A list of atoms to select.
        selected = []

        # Loop over all atoms in the molecule.
        for x in range(0, molecule.nAtoms()):
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)
            # This atom was matched to an atom in the reference.
            if idx in idx_matches:
                # Append to the list of atoms to select.
                selected.append(idx)
                # Copy across the coordinates.
                edit_mol = edit_mol.atom(idx).setProperty("coordinates",
                    reference._sire_object.atom(matches[idx]).property("coordinates")).molecule()
                # Set occupancy and beta factor.
                if matches[idx] in rmsd_indices:
                    # This atom is used to compute the RMSD.
                    edit_mol = edit_mol.atom(idx).setProperty("occupancy", 1.0).molecule()
                    edit_mol = edit_mol.atom(idx).setProperty("beta_factor", 0.0).molecule()
                else:
                    # This atom is used for alignment.
                    edit_mol = edit_mol.atom(idx).setProperty("occupancy", 0.0).molecule()
                    edit_mol = edit_mol.atom(idx).setProperty("beta_factor", 1.0).molecule()

        # Commit the changes.
        new_molecule = edit_mol.commit()

        # Create an AtomSelection.
        selection = new_molecule.selection()

        # Unselect all of the atoms.
        selection.selectNone()

        # Now add all of the atoms that appear in the reference.
        for idx in selected:
            selection.select(idx)

        # Create a partial molecule and extract the atoms.
        partial_molecule = _SireMol.PartialMolecule(new_molecule, selection).extract().molecule()

        # Save the changes to the molecule.
        molecule = _Molecule(partial_molecule)

        # Parse as a PDB file and store the lines.
        pdb = _SireIO.PDB2(molecule.toSystem()._sire_object)
        lines = pdb.toLines()

        # Format for PLUMED.
        self._reference_pdb = lines[1:-2]
        self._reference_pdb.append(lines[-1])

        # Set the "settable" parameters.
        self.setHillWidth(hill_width)
        self.setAlignmentType(alignment_type)
        self.setPeriodicBoundaries(pbc)

        # Set defaults for optional values.
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None

        # Set the optional parameters.
        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)
        if grid is not None:
            self.setGrid(grid)

        # Validate that the state is self-consistent.
        self._validate()

        # Flag that the object has been instantiated, i.e. it is no longer "new".
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.RMSD: "
        string += " reference=%s" % self._reference
        string += ", hill_width=%s" % self._hill_width
        if self._lower_bound is not None:
            string += ", lower_bound=%s" % self._lower_bound
        if self._upper_bound is not None:
            string += ", upper_bound=%s" % self._upper_bound
        if self._grid is not None:
            string += ", grid=%s" % self._grid
        string += ", alignment_type=%s"% self._alignment_type
        string += ", pbc=%s"% self._pbc
        string += ">"
        return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def getReferencePDB(self):
        """Return the reference PDB file as a list of strings.

           Returns
           -------

           pdb : [str]
               The reference PDB file as list of strings.
        """
        return self._reference_pdb

    def setHillWidth(self, hill_width):
        """Set the width of the Gaussian hills used to bias this collective
           variable.

           hill_width : :class:`Length <BioSimSpace.Types.Length>`
               The width of the Gaussian hill.
        """
        if type(hill_width) is not float:
            raise TypeError("'hill_width' must be of type 'float'")

        if hill_width < 0:
            raise ValueError("'hill_width' must be > 0")

        self._hill_width = hill_width

    def getHillWidth(self):
        """Return the width of the Gaussian hill used to bias this collective
           variable.

           Returns
           -------

           hill_width : float
               The width of the Gaussian hill.
        """
        return self._hill_width

    def setAlignmentType(self, alignment_type):
        """Set the RMSD alignment type. Options are "optimal" or "simple".

           Parameters
           ----------

           alignment_type : str
               The RMSD alignment type.
        """
        if type(alignment_type) is not str:
            raise TypeError("'alignment_type' must be of type 'str'")

        alignment_type = alignment_type.lower().replace(" ", "")
        if alignment_type not in ["optimal", "simple"]:
            raise ValueError("'alignment_type' must be 'optimal' or 'simple'.")

        self._alignment_type = alignment_type

    def getAlignmentType(self):
        """Return the RMSD alignment type.

           Returns
           -------

           alignment_type : str
               The RMSD alignment type.
        """
        return self._alignment_type

    def setPeriodicBoundaries(self, pbc):
        """Set whether to use periodic_boundaries when calculating the
           collective variable.

           Parameters
           ----------

           pbc : bool
               Whether to use periodic boundaries conditions.
        """
        if type(pbc) is not bool:
            raise TypeError("'pbc' must be of type 'bool'")
        self._pbc = pbc

    def getPeriodicBoundaries(self):
        """Return whether to take account of periodic boundary conditions
           when computing the collective variable.

           Returns
           -------

           pbc : bool
               Whether to use periodic boundaries conditions.
        """
        return self._pbc

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._lower_bound is not None:
            if type(self._lower_bound.getValue()) is not int and \
               type(self._lower_bound.getValue()) is not float:
                raise TypeError("'lower_bound' must be of type 'int' or 'float'")
        if self._upper_bound is not None:
            if type(self._upper_bound.getValue()) is not int and \
               type(self._upper_bound.getValue()) is not float:
                raise TypeError("'upper_bound' must be of type 'int' or 'float'")
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound.getValue() >= self._upper_bound.getValue():
                raise TypeError("'lower_bound' must less than 'upper_bound'")

        if self._grid is not None:
            if type(self._grid.getMinimum()) is not int and \
               type(self._grid.getMinimum()) is not float:
                raise TypeError("'grid' minimum must be of type 'int' or 'float'")
            if type(self._grid.getMaximum()) is not int and \
               type(self._grid.getMaximum()) is not float:
                raise TypeError("Grid 'maximum' must be of type 'int' or 'float'")
            if self._lower_bound is not None and self._grid.getMinimum() > self._lower_bound.getValue():
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if self._upper_bound is not None and self._grid.getMaximum() < self._upper_bound.getValue():
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")

            # If the number of bins isn't specified, estimate it out from the hill width.
            if self._grid.getBins() is None:
                grid_range = self._grid.getMaximum() - self._grid.getMinimum()
                num_bins = _ceil(5.0 * (grid_range / self._hill_width))
                self._grid.setBins(num_bins)

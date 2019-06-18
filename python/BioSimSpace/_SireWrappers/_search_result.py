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
A thin wrapper around Sire.Mol.SelectResult. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["SearchResult"]

import Sire.Mol as _SireMol

class SearchResult():
    """A thin wrapper around Sire.Mol.SelectResult."""

    def __init__(self, select_result):
        """Constructor.

           Parameters
           ----------

           select_result : Sire.Mol.SelectResult
               The Sire select result object.
        """

        # Check that the select_result is valid.

        if type(select_result) is SearchResult:
            select_result = select_result._sire_object
        elif type(select_result) is _SireMol.SelectResult:
            pass
        else:
            raise TypeError("'select_result' must be of type 'BioSimSpace._SireWrappers.SearchResult' "
                            "or 'Sire.Mol.SelectResult'")

        # Store the Sire select result.
        self._sire_object = select_result.__deepcopy__()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.SearchResult: nResults=%d>" % self.nResults()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.SearchResult: nResults=%d>" % self.nResults()

    def __len__(self):
        """Return the number of results.

           Returns
           -------

           num_results : int
               The number of search results.
        """
        return self.nResults()

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object == other._sire_object
        else:
            return False

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object != other._sire_object
        else:
            return False

    def copy(self):
        """Create a copy of this object.

           Returns
           -------

           search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
               A copy of the object.
        """
        return SearchResult(self)

    def nResults(self):
        """Return the number of results.

           Returns
           -------

           num_results : int
               The number of search results.
        """
        return len(self._sire_object)

    def getResults(self):
        """Return a list containing the results.

           Returns
           -------

           results : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                      :class:`Residue <BioSimSpace._SireWrappers.Residue>`, \
                      :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, ...]
               A list of objects matching the search query.
        """

        results = []

        # Loop over each result and determine the type.
        for x in self._sire_object:
            # Atom.
            if type(x) is _SireMol.Atom:
                results.append(_Atom(x))
            # Residue.
            if type(x) is _SireMol.Residue:
                results.append(_Residue(x))
            # Molecule.
            if type(x) is _SireMol.Molecule:
                # If the molecule contains a single atom, then convert to an atom.
                if x.nAtoms() == 1:
                    results.append(_Atom(x.atom()))
                # If there's a single residue, the convert to a residue.
                elif x.nResidues() == 1:
                    results.append(_Residue(x.residue()))
                # Otherwise, append the molecule.
                else:
                    results.append(_Molecule(x))
            # Unsupported.
            else:
                pass

        return results

    def _getSireObject(self):
        """Return the underlying Sire object.

           Returns
           -------

           object : Sire.Mol.SelectResult
        """
        return self._sire_object

# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecule import Molecule as _Molecule
from ._residue import Residue as _Residue

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

        if type(select_result) is not _SireMol.SelectResult:
            raise TypeError("'select_result' must be of type 'Sire.Mol.SelectResult'")

        # Store the Sire select result.
        self._select_result = select_result

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

    def nResults(self):
        """Return the number of results.

           Returns
           -------

           num_results : int
               The number of search results.
        """
        return len(self._select_result)

    def getResults(self):
        """Return a list containing the results.

           Returns
           -------

           results : [:class:`Atom <BioSimSpace._SireWrappers.Atom`,
                      :class:`Residue <BioSimSpace._SireWrappers.Residue`,
                      :class:`Molecule <BioSimSpace._SireWrappers.Molecule`, ...]
               A list of objects matching the search query.
        """

        results = []

        # Loop over each result and determine the type.
        for x in self._select_result:
            # Atom.
            if type(x) is _SireMol.Atom:
                results.append(_Atom(x))
            # Residue.
            if type(x) is _SireMol.Residue:
                results.append(_Residue(x))
            # Molecule.
            if type(x) is _SireMol.Molecule:
                # Sometimes residue based searches are converted to molecules.
                # Attempt to cast as a residue. This will work if the molecule
                # contains a single residue.
                try:
                    results.append(_Residue(x.residue()))
                except:
                    results.append(_Molecule(x))
            # Unsupported.
            else:
                pass

        return results

# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecule import Molecule as _Molecule
from ._residue import Residue as _Residue

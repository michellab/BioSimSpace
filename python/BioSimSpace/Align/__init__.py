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
Functionality for aligning molecules.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Mol as _SireMol

from .._SireWrappers import Molecule as _Molecule

import BioSimSpace.Units as _Units

def matchAtoms(molecule0,
               molecule1,
               scoring_function="RMSD",
               matches=1,
               prematch={},
               timeout=5*_Units.Time.second,
               match_light=True,
               map0={},
               map1={},
               verbose=False):
    """Find mappings from atoms in molecule0 to those in molecule1.

       Positional arguments:

       molecule0        -- The reference molecule.
       molecule1        -- The target molecule.

       Keyword arguments:

       scoring_function -- The scoring function used to match atoms. Available
                           options are: "something", "something else", ...
       matches          -- The maximum number of matches to return. (Sorted in order of score).
       prematch         -- A pre-match to use as the basis of the search.
       timeout          -- The timeout for the matching algorithm.
       match_light      -- Whether to match light atoms.
       map0             -- A dictionary that maps "properties" in molecule0 to their user
                           defined values. This allows the user to refer to properties
                           with their own naming scheme, e.g. { "charge" : "my-charge" }
       map1             -- A dictionary that maps "properties" in molecule1 to their user
                           defined values.
       verbose          -- Whether to print status information from the matcher.
    """

    # A list of supported scoring functions.
    scoring_functions = ["RMSD"]

    # Validate input.

    if type(molecule0) is not _Molecule:
        raise TypeError("'molecule0' must be of type 'BioSimSpace.Molecule'")

    if type(molecule1) is not _Molecule:
        raise TypeError("'molecule1' must be of type 'BioSimSpace.Molecule'")

    if type(scoring_function) is not str:
        raise TypeError("'scoring_function' must be of type 'str'")
    else:
        if not scoring_function.replace(" ", "").upper() in scoring_functions:
            raise ValueError("Unsupported scoring function '%s'. Options are: %s"
                % (scoring_function, scoring_functions))

    if type(matches) is not int:
        raise TypeError("'matches' must be of type 'int'")
    else:
        if matches < 0:
            raise ValueError("'matches' must be positive!")

    if type(prematch) is not dict:
        raise TypeError("'prematch' must be of type 'dict'")
    else:
        for idx0, idx1 in prematch.items():
            if type(idx0) is not _SireMol.AtomIdx or type(idx1) is not _SireMol.AtomIdx:
                raise TypeError("'prematch' dictionary key:value pairs must be of type 'Sire.Mol.AtomIdx'")

    if type(timeout) is not _Units.Time._Time:
        raise TypeError("'timeout' must be of type 'BioSimSpace.Types.Time'")

    if type(match_light) is not bool:
        raise TypeError("'match_light' must be of type 'bool'")

    if type(map0) is not dict:
        raise TypeError("'map0' must be of type 'dict'")

    if type(map1) is not dict:
        raise TypeError("'map1' must be of type 'dict'")

    if type(verbose) is not bool:
        raise TypeError("'verbose' must be of type 'bool'")

    # Extract the Sire molecule from each BioSimSpace molecule.
    mol0 = molecule0._getSireMolecule()
    mol1 = molecule1._getSireMolecule()

    # Convert the timeout to a Sire unit.
    timeout = timeout.magnitude() * timeout._supported_units[timeout.unit()]

    # Find all of the best maximum common substructure matches.
    mappings = ( mol0.evaluate()
                     .findMCSmatches(mol1, _SireMol.AtomResultMatcher(prematch),
                         timeout, match_light, map0, map1, verbose)
               )

    # No matches!
    if len(mappings) == 0:
        return None

    # Score the mappings and return them in sorted order (best to worst).
    # For now we default to RMSD scoring, since it's the only option.
    else:
        # Return the best match.
        if matches == 1:
            return _score_rmsd(mol0, mol1, mappings)[0]
        # Return a list of matches from best to worst.
        else:
            return _score_rmsd(mol0, mol1, mappings)[0:matches]

def rmsdAlign(molecule0, molecule1, mapping):
    """Align atoms in molecule0 to those in molecule1 using the mapping
       between matched atom indices. The molecule is aligned based on
       a root mean squared displacement (RMSD) fit to find the optimal
       translation vector (as opposed to merely taking the difference of
       centroids).

       Positional arguments:

       molecule0 -- The reference molecule.
       molecule1 -- The target molecule.
       mapping   -- A dictionary mapping atoms in molecule0 to those in molecule1.
    """

    if type(molecule0) is not _Molecule:
        raise TypeError("'molecule0' must be of type 'BioSimSpace.Molecule'")

    if type(molecule1) is not _Molecule:
        raise TypeError("'molecule1' must be of type 'BioSimSpace.Molecule'")

    if type(mapping) is not dict:
        raise TypeError("'mapping' must be of type 'dict'.")
    else:
        # Make sure all key/value pairs are of type AtomIdx.
        for idx0, idx1 in mapping.items():
            if type(idx0) is not _SireMol.AtomIdx or type(idx1) is not _SireMol.AtomIdx:
                raise TypeError("key:value pairs in 'mapping' must be of type 'Sire.Mol.AtomIdx'")

    # Extract the Sire molecule from each BioSimSpace molecule.
    mol0 = molecule0._getSireMolecule()
    mol1 = molecule1._getSireMolecule()

    # Perform the alignment, mol0 to mol1.
    mol0 = mol0.move().align(mol1, _SireMol.AtomResultMatcher(mapping)).molecule()

    # Return the aligned molecule.
    return _Molecule(mol0)

def _score_rmsd(molecule0, molecule1, mappings):
    """Internal function to score atom mappings based on the root mean squared
       displacement between matched atoms that are aligned based on each mapping.
       Returns the mappings sorted based on their score from best to worst.

       Positional arguments:

       molecule0 -- The reference molecule.
       molecule1 -- The target molecule.
       mappings  -- A list of dictionaries mapping atoms in molecule0 to
                    those in molecule1.
    """

    # Initialise a list of scores.
    scores = []

    # Loop over all mappings.
    for map in mappings:
        # Align molecule0 to molecule1 based on the mapping.
        aligned_mol = molecule0.move().align(molecule1, _SireMol.AtomResultMatcher(map))

        # Compute the RMSD between the two molecules and add to the scores.
        scores.append(molecule0.evaluate().rmsd(aligned_mol).value())

    # Sort the scores and return the sorted keys.
    keys = sorted(range(len(scores)), key=lambda k: scores[k])

    # Return the sorted mappings.
    return [mappings[x] for x in keys]

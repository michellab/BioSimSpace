from pytest import approx as _approx

import os.path as _path
import random as _random
import string as _string

from Sire import Base as _SireBase
from Sire import CAS as _SireCAS
from Sire import IO as _SireIO
from Sire import MM as _SireMM
from Sire import Mol as _SireMol
from Sire import System as _SireSystem
from Sire import Units as _SireUnits

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace.Types import Coordinate as _Coordinate
from BioSimSpace.Types import Length as _Length

def _random_suffix(basename, size=4, chars=_string.ascii_uppercase + _string.digits):
    """Internal helper function to generate a random atom name suffix to avoid
       naming clashes.

       Adapted from:
       https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python

       Parameters
       ----------

       basename : str
           The base string to which a suffix will be appended.

       size : int
           The maximum width of the string, i.e. len(basename + suffix).

       chars : str
           The set of characters to include in the suffix.

       Returns
       -------

       suffix : str
           The randomly generated suffix.
    """

    basename_size = len(basename)
    if basename_size >= size:
        raise ValueError("Cannot generate suffix for basename '%s'. " % basename
                       + "AMBER atom names can only be 4 characters wide.")
    return "".join(_random.choice(chars) for _ in range(size-basename_size))
    
def _has_pert_atom(idxs, pert_idxs):
    """Internal function to check whether a potential contains perturbed atoms.

       Parameters
       ----------

       idxs : [AtomIdx]
           A list of atom indices involved in the potential.

       pert_idxs : [AtomIdx]
           A list of atom indices that are perturbed.

       Returns
       -------

       has_pert_atom : bool
           Whether the potential includes a perturbed atom.
    """

    for idx in idxs:
        if idx in pert_idxs:
            return True

    return False

def _has_dummy(mol, idxs, is_lambda1=False):
    """Internal function to check whether any atom is a dummy.

       Parameters
       ----------

       mol : Sire.Mol.Molecule
           The molecule.

       idxs : [AtomIdx]
           A list of atom indices.

       is_lambda1 : bool
           Whether to check the lambda = 1 state.

       Returns
       -------

       has_dummy : bool
           Whether a dummy atom is present.
    """

    # Set the element property associated with the end state.
    if is_lambda1:
        prop = "element1"
    else:
        prop = "element0"

    dummy = _SireMol.Element(0)

    # Check whether an of the atoms is a dummy.
    for idx in idxs:
        if mol.atom(idx).property(prop) == dummy:
            return True

    return False

def _is_dummy(mol, idxs, is_lambda1=False):
    """Internal function to return whether each atom is a dummy.

       Parameters
       ----------

       mol : Sire.Mol.Molecule
           The molecule.

       idxs : [AtomIdx]
           A list of atom indices.

       is_lambda1 : bool
           Whether to check the lambda = 1 state.

       Returns
       -------

       is_dummy : [bool]
           Whether each atom is a dummy.
    """

    # Set the element property associated with the end state.
    if is_lambda1:
        prop = "element1"
    else:
        prop = "element0"

    # Store a dummy element.
    dummy = _SireMol.Element(0)

    # Initialise a list to store the state of each atom.
    is_dummy = []

    # Check whether each of the atoms is a dummy.
    for idx in idxs:
        is_dummy.append(mol.atom(idx).property(prop) == dummy)

    return is_dummy

def _random_suffix(basename, size=4, chars=_string.ascii_uppercase + _string.digits):
    """Internal helper function to generate a random atom name suffix to avoid
       naming clashes.

       Adapted from:
       https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python

       Parameters
       ----------

       basename : str
           The base string to which a suffix will be appended.

       size : int
           The maximum width of the string, i.e. len(basename + suffix).

       chars : str
           The set of characters to include in the suffix.

       Returns
       -------

       suffix : str
           The randomly generated suffix.
    """

    basename_size = len(basename)
    if basename_size >= size:
        raise ValueError("Cannot generate suffix for basename '%s'. " % basename
                       + "AMBER atom names can only be 4 characters wide.")
    return "".join(_random.choice(chars) for _ in range(size-basename_size))

def _is_ring_broken(conn0, conn1, idx0, idy0, idx1, idy1):
    """Internal function to test whether a perturbation changes the connectivity
       around two atoms such that a ring is broken.

       Parameters
       ----------

       conn0 : Sire.Mol.Connectivity
           The connectivity object for the first end state.

       conn1 : Sire.Mol.Connectivity
           The connectivity object for the second end state.

       idx0 : Sire.Mol.AtomIdx
           The index of the first atom in the first state.

       idy0 : Sire.Mol.AtomIdx
           The index of the second atom in the first state.

       idx1 : Sire.Mol.AtomIdx
           The index of the first atom in the second state.

       idy1 : Sire.Mol.AtomIdx
           The index of the second atom in the second state.
    """

    # Have we opened/closed a ring? This means that both atoms are part of a
    # ring in one end state (either in it, or on it), whereas at least one
    # are the result of changes in ring size, where atoms remain in or on a
    # ring in both end states.

    # Whether each atom is in a ring in both end states.
    in_ring_idx0 = conn0.inRing(idx0)
    in_ring_idy0 = conn0.inRing(idy0)
    in_ring_idx1 = conn1.inRing(idx1)
    in_ring_idy1 = conn1.inRing(idy1)

    # Whether each atom is on a ring in both end states.
    on_ring_idx0 = _onRing(idx0, conn0)
    on_ring_idy0 = _onRing(idy0, conn0)
    on_ring_idx1 = _onRing(idx1, conn1)
    on_ring_idy1 = _onRing(idy1, conn1)

    # Both atoms are in a ring in one end state and at least one isn't in the other.
    if (in_ring_idx0 & in_ring_idy0) ^ (in_ring_idx1 & in_ring_idy1):
        return True

    # Both atoms are on a ring in one end state and at least one isn't in the other.
    if ((on_ring_idx0 & on_ring_idy0 & (conn0.connectionType(idx0, idy0) == 4))
        ^ (on_ring_idx1 & on_ring_idy1 & (conn1.connectionType(idx1, idy1) == 4))):
        return True

    # Both atoms are in or on a ring in one state and at least one isn't in the other.
    if (((in_ring_idx0 | on_ring_idx0) & (in_ring_idy0 | on_ring_idy0) & (conn0.connectionType(idx0, idy0) == 3)) ^
        ((in_ring_idx1 | on_ring_idx1) & (in_ring_idy1 | on_ring_idy1) & (conn1.connectionType(idx1, idy1) == 3))):
        iscn0 = set(conn0.connectionsTo(idx0)).intersection(set(conn0.connectionsTo(idy0)))
        if (len(iscn0) != 1):
            return True
        common_idx = iscn0.pop()
        in_ring_bond0 = (conn0.inRing(idx0, common_idx) | conn0.inRing(idy0, common_idx))
        iscn1 = set(conn1.connectionsTo(idx1)).intersection(set(conn1.connectionsTo(idy1)))
        if (len(iscn1) != 1):
            return True
        common_idx = iscn1.pop()
        in_ring_bond1 = (conn1.inRing(idx1, common_idx) | conn1.inRing(idy1, common_idx))
        if (in_ring_bond0 ^ in_ring_bond1):
            return True

    # If we get this far, then a ring wasn't broken.
    return False

def _is_ring_size_changed(conn0, conn1, idx0, idy0, idx1, idy1, max_ring_size=12):
    """Internal function to test whether a perturbation changes the connectivity
       around two atoms such that a ring changes size.

       Parameters
       ----------

       conn0 : Sire.Mol.Connectivity
           The connectivity object for the first end state.

       conn1 : Sire.Mol.Connectivity
           The connectivity object for the second end state.

       idx0 : Sire.Mol.AtomIdx
           The index of the first atom in the first state.

       idy0 : Sire.Mol.AtomIdx
           The index of the second atom in the first state.

       idx1 : Sire.Mol.AtomIdx
           The index of the first atom in the second state.

       idy1 : Sire.Mol.AtomIdx
           The index of the second atom in the second state.

       max_ring_size : int
           The maximum size of what is considered to be a ring.
    """

    # Have a ring changed size? If so, then the minimum path size between
    # two atoms will have changed.

    # Work out the paths connecting the atoms in the two end states.
    paths0 = conn0.findPaths(idx0, idy0, max_ring_size)
    paths1 = conn1.findPaths(idx1, idy1, max_ring_size)

    # Initalise the ring size in each end state.
    ring0 = None
    ring1 = None

    # Determine the minimum path in the lambda = 0 state.
    if len(paths0) > 1:
        path_lengths0 = []
        for path in paths0:
            path_lengths0.append(len(path))
        ring0 = min(path_lengths0)

    if ring0 is None:
        return False

    # Determine the minimum path in the lambda = 1 state.
    if len(paths1) > 1:
        path_lengths1 = []
        for path in paths1:
            path_lengths1.append(len(path))
        ring1 = min(path_lengths1)

    # Return whether the ring has changed size.
    if ring1:
        return ring0 != ring1
    else:
        return False

def _onRing(idx, conn):
    """Internal function to test whether an atom is adjacent to a ring.

       Parameters
       ----------

       idx : Sire.Mol.AtomIdx
           The index of the atom

       conn : Sire.Mol.Connectivity
           The connectivity object.

       Returns
       -------

       is_on_ring : bool
           Whether the atom is adjacent to a ring.
    """

    # Loop over all atoms connected to this atom.
    for x in conn.connectionsTo(idx):
        # The neighbour is in a ring.
        if conn.inRing(x) and (not conn.inRing(x, idx)):
            return True

    # If we get this far, then the atom is not adjacent to a ring.
    return False
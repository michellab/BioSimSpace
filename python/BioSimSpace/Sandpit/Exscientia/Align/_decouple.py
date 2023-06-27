import warnings

from sire.legacy import Base as _SireBase

from .._SireWrappers import Molecule as _Molecule
from .._Exceptions import IncompatibleError as _IncompatibleError


__all__ = ["decouple"]


def decouple(molecule, charge=(True, False), LJ=(True, False), intramol=True):
    """
    Make the molecule as being decoupled, where the interactions with the
    rest of the environment are removed, or annihilate, where the interactions
    within the molecule are removed as well (choose this mode with
    intramol=False).

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The molecule to be decoupled or annihilated.

    charge : Tuple[bool, bool]
        A Tuple of length two defines the charge of the molecule at the
        start and the end of the transformation. This allows user to
        selectively turn on or off the charge interactions
        e.g. (True, False)

    LJ : Tuple[bool, bool]
        A Tuple of length two defines the LJ of the molecule at the
        start and the end of the transformation. This allows user to
        selectively turn on or off the LJ interactions
        e.g. (True, False)

    intramol : bool
        Whether to couple the intra-molecule forces. Setting to ``False``
        means the intra-molecule forces will *not* be changed by the lambda
        thereby decouple the molecule instead of annihilate it.

    Returns
    -------

    decoupled : Sire.Mol.Molecule
        The molecule marked as being decoupled.
    """
    # Validate input.

    if not isinstance(molecule, _Molecule):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    # Cannot decouple a perturbable molecule.
    if molecule.isDecoupled():
        raise _IncompatibleError("'molecule' has already been decoupled!")

    for field, name in zip((charge, LJ), ("charge", "LJ")):
        try:
            if len(field) != 2:
                raise ValueError(f"{name} must only have two values.")
        except TypeError:
            raise TypeError(f"{name} has to be a iterable.")

        for value in field:
            if not isinstance(value, bool):
                raise ValueError(f"{value} in {name} must be bool.")

    if not isinstance(intramol, bool):
        raise TypeError("'intramol' must be of type 'bool'")

    # Create a copy of this molecule.
    mol = _Molecule(molecule)
    mol_sire = mol._sire_object

    # Edit the molecule
    mol_edit = mol_sire.edit()

    mol_edit.setProperty("decouple", {"charge": charge, "LJ": LJ, "intramol": intramol})

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()

    return mol

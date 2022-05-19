import warnings

from .._SireWrappers import Molecule as _Molecule
from .._Exceptions import IncompatibleError as _IncompatibleError


def decouple(molecule, property_map0={}, property_map1={}, intramol=True):
    """Make the molecule as being decoupled, where the interactions with the
    rest of the environment are removed, or annihilate, where the interactions
    within the molecule are removed as well (choose this mode with
    intramol=False).

        Parameters
        ----------

        molecule0 : BioSimSpace._SireWrappers.Molecule
            The molecule to be decoupled or annihilated.

        property_map0 : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values at the start of the transformation. This allows
            user to selectively turn on or off the charge or the vdw
            interactions e.g. { "charge" : True, "vdw" : True}

        property_map1 : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values at the end of the transformation. This allows
            user to selectively turn on or off the charge or the vdw
            interactions e.g. { "charge" : False, "vdw" : False}

        intramol : bool
            Whether to remove the intra-molecule forces, thereby annihilate the
            molecule instead of decouple it.

        Returns
        -------

        decoupled : Sire.Mol.Molecule
            The molecule marked as being decoupled.
    """
    # Validate input.

    if not isinstance(molecule, _Molecule):
        raise TypeError("'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    # Cannot decouple a perturbable molecule.
    if molecule._is_perturbable:
        raise _IncompatibleError("'molecule0' has already been merged!")

    # Check the keys in the property_map
    for property_map in (property_map0, property_map1):
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        for key in property_map:
            if not key in ('charge', 'vdw'):
                warnings.warn(f'Key {key} not supported for decouple, will be '
                              f'ignored. The recognised keys are '
                              f'("charge", "vdw")')

    if not isinstance(intramol, bool):
        raise TypeError("'intramol' must be of type 'bool'")

    # Create a copy of this molecule.
    mol = _Molecule(molecule)


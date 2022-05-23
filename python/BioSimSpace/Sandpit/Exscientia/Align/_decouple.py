import warnings

from Sire import Base as _SireBase

from .._SireWrappers import Molecule as _Molecule
from .._Exceptions import IncompatibleError as _IncompatibleError


def decouple(molecule, property_map0=None, property_map1=None, intramol=True):
    """Make the molecule as being decoupled, where the interactions with the
    rest of the environment are removed, or annihilate, where the interactions
    within the molecule are removed as well (choose this mode with
    intramol=False).

        Parameters
        ----------

        molecule : BioSimSpace._SireWrappers.Molecule
            The molecule to be decoupled or annihilated.

        property_map0 : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values at the start of the transformation. This allows
            user to selectively turn on or off the charge or the vdw
            interactions e.g. { "charge" : True, "LJ" : True}

        property_map1 : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values at the end of the transformation. This allows
            user to selectively turn on or off the charge or the vdw
            interactions e.g. { "charge" : False, "LJ" : False}

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
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    # Cannot decouple a perturbable molecule.
    if molecule.isDecoupled():
        raise _IncompatibleError("'molecule' has already been decoupled!")

    if property_map0 is None and property_map1 is None:
        property_map0 = {"charge": True, "LJ": True}
        property_map1 = {"charge": False, "LJ": False}
    if property_map0 is None:
        property_map0 = {}
    if property_map1 is None:
        property_map1 = {}

    # Check the keys in the property_map
    for property_map in (property_map0, property_map1):
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        for key in property_map:
            if not key in ('charge', 'LJ'):
                warnings.warn(f'Key {key} not supported for decouple, will be '
                              f'ignored. The recognised keys are '
                              f'("charge", "LJ")')

    if not isinstance(intramol, bool):
        raise TypeError("'intramol' must be of type 'bool'")

    # Create a copy of this molecule.
    mol = _Molecule(molecule)
    mol_sire = mol._sire_object

    # Edit the molecule
    mol_edit = mol_sire.edit()

    # Set the start and the end state
    mol_edit.setProperty("charge0", _SireBase.wrap(property_map0.get("charge",True)))
    mol_edit.setProperty("charge1", _SireBase.wrap(property_map1.get("charge", True)))
    mol_edit.setProperty("LJ0", _SireBase.wrap(property_map0.get("LJ", True)))
    mol_edit.setProperty("LJ1", _SireBase.wrap(property_map1.get("LJ", True)))

    mol_edit.setProperty("annihilated", _SireBase.wrap(not intramol))

    # Flag that this molecule is decoupled.
    mol_edit.setProperty("is_decoupled", _SireBase.wrap(True))

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()

    return mol


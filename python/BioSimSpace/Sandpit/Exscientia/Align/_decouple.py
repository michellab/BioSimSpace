import warnings

from Sire import Base as _SireBase
from Sire import Units as _SireUnits
from Sire import MM as _SireMM

from .._SireWrappers import Molecule as _Molecule
from .._Exceptions import IncompatibleError as _IncompatibleError

__all__ = ["decouple"]

def decouple(molecule, property_map={}, intramol=True):
    """Make the molecule as being decoupled, where the interactions with the
    rest of the environment are removed, or annihilate, where the interactions
    within the molecule are removed as well (choose this mode with
    intramol=False).

        Parameters
        ----------

        molecule : BioSimSpace._SireWrappers.Molecule
            The molecule to be decoupled or annihilated.

        property_map : dict
            A dictionary that maps "properties" to their user defined values. 
            This allows the user to refer to properties with their own naming 
            scheme, e.g. { "charge" : "my-charge" }

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
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    # Cannot decouple a perturbable molecule.
    if molecule.isDecoupled():
        raise _IncompatibleError("'molecule' has already been decoupled!")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Invert the user property mappings.
    inv_property_map = {v: k for k, v in property_map.items()}

    # Create a copy of this molecule and Sire object to check properties
    mol = _Molecule(molecule)
    mol_sire = mol._sire_object

    # Get the user name for the required properties.
    ff = inv_property_map.get("forcefield", "forcefield")
    LJ = inv_property_map.get("LJ", "LJ")
    charge = inv_property_map.get("charge", "charge")
    ambertype = inv_property_map.get("ambertype", "ambertype")

    # Check for missing information
    if not mol_sire.hasProperty(ff):
        raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule'!")
    if not mol_sire.hasProperty(LJ):
        raise _IncompatibleError("Cannot determine LJ terms for molecule")
    if not mol_sire.hasProperty(charge):
        raise _IncompatibleError("Cannot determine charges for molecule")

    # Check for ambertype property (optional)
    has_ambertype = True
    if not mol_sire.hasProperty(ambertype):
        has_ambertype = False

    if not isinstance(intramol, bool):
        raise TypeError("'intramol' must be of type 'bool'")

    # Edit the molecule
    mol_edit = mol_sire.edit()

    # Set starting properties based on fully-interacting molecule
    mol_edit.setProperty("charge0", molecule._sire_object.property(charge))
    mol_edit.setProperty("LJ0", molecule._sire_object.property(LJ))
    if has_ambertype:
        mol_edit.setProperty("ambertype0", molecule._sire_object.property(ambertype))

    # Set final charges and LJ terms to 0 and (if required) ambertypes to du
    for atom in mol_sire.atoms():
            mol_edit = mol_edit.atom(atom.index()).setProperty("charge1", 0*_SireUnits.e_charge).molecule()
            mol_edit = mol_edit.atom(atom.index()).setProperty("LJ1", _SireMM.LJParameter()).molecule()
            if has_ambertype:
                mol_edit = mol_edit.atom(atom.index()).setProperty("ambertype1", "du").molecule()

    mol_edit.setProperty("annihilated", _SireBase.wrap(intramol))

    # Flag that this molecule is decoupled.
    mol_edit.setProperty("is_decoupled", _SireBase.wrap(True))

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()

    return mol


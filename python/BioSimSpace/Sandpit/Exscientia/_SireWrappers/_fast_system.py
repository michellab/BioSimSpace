# TODO: remove this file when BSS performance is sufficient for large systems

from Sire import IO as _SireIO
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace._SireWrappers import Molecules as _Molecules
from BioSimSpace._SireWrappers import System as _System


def _fastSystemInit(molecules):
    """Add a molecule, or list of molecules to the system.

       Parameters
       ----------

       molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                   :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                   [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`], \
                   :class:`System <BioSimSpace._SireWrappers.System>`
          A Molecule, Molecules object, a list of Molecule objects, or a System containing molecules.
    """

    system = _SireSystem.System("BioSimSpace System")
    molgrp = _SireMol.MoleculeGroup("all")

    # Add the molecule groups in a quick way.
    # A BioSimSpace System object.
    if type(molecules) is _System:
        molecules = molecules.getMolecules()

    # Convert tuple to a list.
    if type(molecules) is tuple:
        molecules = list(molecules)

    # A Molecule object.
    if type(molecules) is _Molecule:
        molecules = [molecules]

    if type(molecules) is _SireMol.Molecule:
        molecules = [_Molecule(molecules)]

    if type(molecules) is not _Molecules:
        molecules = _Molecules(molecules)

    molgrp.add(molecules._sire_object)

    # Initialise the new system and renumber its molecules in a quick way.
    system.add(molgrp)
    _SireIO.renumberConstituents(system, system.nMolecules())
    system = _System(system)

    if type(molecules) is _System:
        # Copy all the properties from the old system.
        for prop in molecules._sire_object.propertyKeys():
            val = system._sire_object.property(prop)
            system._sire_object.setProperty(prop, val)

    return system

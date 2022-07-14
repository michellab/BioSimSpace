from Sire.MM import InternalFF, IntraCLJFF, IntraFF
from Sire.Mol import AtomIdx, PartialMolecule
from Sire.Maths import Vector

import BioSimSpace as BSS

import pytest

# Parameterise the function with a set of valid atom pre-matches.
@pytest.mark.parametrize("prematch", [{3 : 1},
                                      {5 : 9},
                                      {4 : 5},
                                      {1 : 0}])
def test_prematch(prematch):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand01*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand02*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules that contains the prematch.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second,
                                   prematch=prematch)

    # Check that the prematch key:value pair is in the mapping.
    for key, value in prematch.items():
        assert mapping[key] == value

# Parameterise the function with a set of invalid atom pre-matches.
@pytest.mark.parametrize("prematch", [{-1 :  1},
                                      {50 :  9},
                                      { 4 : 48},
                                      { 1 : -1}])
def test_invalid_prematch(prematch):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand01*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand02*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Assert that the invalid prematch raises a ValueError.
    with pytest.raises(ValueError):
        mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second,
                                       prematch=prematch)

def test_merge():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand38*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Create the merged molecule.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

    # Store the number of atoms in m0.
    n0 = m0._sire_object.nAtoms()

    # Test that the intramolecular energies area the same.

    # IntraCLJFF:
    #  Old interface. Uses the "intrascale" matrix. Validate that this
    #  is consistent.
    # IntraFF:
    #  New interface. Uses atom "connectivity". Validate that the bonding
    #  is consistent.

    intraclj0 = IntraCLJFF("intraclj")
    intraclj0.add(m0._sire_object)

    intraff0 = IntraFF("intraclj")
    intraff0.add(m0._sire_object)

    intraclj1 = IntraCLJFF("intraclj")
    intraclj1.add(m1._sire_object)

    intraff1 = IntraFF("intraclj")
    intraff1.add(m1._sire_object)

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    # Create maps between property names: { "prop" : "prop0" }, { "prop" : "prop1" }
    pmap0 = {}
    pmap1 = {}
    for prop in m2._sire_object.propertyKeys():
        if prop[-1] == "0":
            pmap0[prop[:-1]] = prop
        elif prop[-1] == "1":
            pmap1[prop[:-1]] = prop

    intraclj2.add(m2._sire_object, pmap0)
    intraff2.add(m2._sire_object, pmap0)

    assert intraclj0.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff0.energy().value() == pytest.approx(intraff2.energy().value())

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    intraclj2.add(m2._sire_object, pmap1)
    intraff2.add(m2._sire_object, pmap1)

    assert intraclj1.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff1.energy().value() == pytest.approx(intraff2.energy().value())

    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.

    internalff0 = InternalFF("internal")
    internalff0.setStrict(True)
    internalff0.add(m0._sire_object)

    internalff1 = InternalFF("internal")
    internalff1.setStrict(True)
    internalff1.add(m1._sire_object)

    # First extract a partial molecule using the atoms from molecule0 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for atom in m0._sire_object.atoms():
        selection.select(atom.index())
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap0)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=0 end state.
    amber_mol, _ = m2._extractMolecule()

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Now extract a partial molecule using the atoms from molecule1 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for idx in mapping.keys():
        selection.select(AtomIdx(idx))
    for idx in range(n0, m2._sire_object.nAtoms()):
        selection.select(AtomIdx(idx))
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap1)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=1 end state.
    amber_mol, _ = m2._extractMolecule(is_lambda1=True)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

def test_roi_merge():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/wild*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/mutated*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Extract the residue name and index
    def get_res_info(mol):
        mol_res_name = []
        mol_res_idx = []
        for bss_res in mol.getResidues():
            mol_res_name.append(bss_res.name())
            mol_res_idx.append([x.index() for x in bss_res.getAtoms()])
        return mol_res_name, mol_res_idx

    def recover_mapping(mut1_idx, mut2_idx, mapping):
        new_mapping = {}
        for k, v in mapping.items():
            new_mapping[mut1_idx[k]] = mut2_idx[v]
        return new_mapping
    
    def update_coordinate(mol, coord_dict):
        edit_mol = mol._sire_object.edit()
        for idx in range(0, mol.nAtoms()):
            if idx in coord_dict:
                vec = coord_dict[idx]
                vec = Vector(vec.x().angstroms().value(), 
                             vec.y().angstroms().value(),
                             vec.z().angstroms().value())
                edit_mol = edit_mol.atom(AtomIdx(idx)).setProperty("coordinates", vec).molecule()
        mol._sire_object = edit_mol.commit()

    
    mol0_res_name, mol0_res_idx = get_res_info(m0) 
    mol1_res_name, mol1_res_idx = get_res_info(m1)
     
    # Get the best mapping between the molecules.
    mapping = {}
    for i in range(len(mol0_res_name)):
        # The mutated residue would have different name
        if mol0_res_name[i] == mol1_res_name[i]:
            for k, v in zip(mol0_res_idx[i], mol1_res_idx[i]):
                mapping[k] = v
        else:
            # Assume the atom order is kept during conversion
            mut0 = m0.getResidues()[i].toMolecule()
            mut1 = m1.getResidues()[i].toMolecule()
            mut0_idx = [x.index() for x in m0.getResidues()[i].getAtoms()]
            mut1_idx = [x.index() for x in m1.getResidues()[i].getAtoms()]

            mapping_mut = BSS.Align.matchAtoms(mut0, mut1)

            # We should also translate the moelcule according to the mapping
            # or there will be odd dummy atom in the merged topology
            mut0 = BSS.Align.rmsdAlign(mut0, mut1, mapping=mapping_mut)
            aligned_coordinates = mut0.coordinates()
            coord_dict = dict(zip(mut0_idx, aligned_coordinates))

            mapping_recovered = recover_mapping(mut0_idx, mut1_idx, mapping_mut)
            mut_idx = i
    mapping.update(mapping_recovered)
    update_coordinate(m0, coord_dict)
    
    # Create the merged molecule.
    m2 = BSS.Align.merge(m0, m1, mapping=mapping, force=True, roi=[mut0_idx, mut1_idx], mut_idx=mut_idx)

    # Store the number of atoms in m0.
    n0 = m0._sire_object.nAtoms()

    # Test that the intramolecular energies are the same.

    # IntraCLJFF:
    #  Old interface. Uses the "intrascale" matrix. Validate that this
    #  is consistent.
    # IntraFF:
    #  New interface. Uses atom "connectivity". Validate that the bonding
    #  is consistent.

    intraclj0 = IntraCLJFF("intraclj")
    intraclj0.add(m0._sire_object)

    intraff0 = IntraFF("intraclj")
    intraff0.add(m0._sire_object)

    intraclj1 = IntraCLJFF("intraclj")
    intraclj1.add(m1._sire_object)

    intraff1 = IntraFF("intraclj")
    intraff1.add(m1._sire_object)

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    # Create maps between property names: { "prop" : "prop0" }, { "prop" : "prop1" }
    pmap0 = {}
    pmap1 = {}
    for prop in m2._sire_object.propertyKeys():
        if prop[-1] == "0":
            pmap0[prop[:-1]] = prop
        elif prop[-1] == "1":
            pmap1[prop[:-1]] = prop

    intraclj2.add(m2._sire_object, pmap0)
    intraff2.add(m2._sire_object, pmap0)

    assert intraclj0.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff0.energy().value() == pytest.approx(intraff2.energy().value())

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    intraclj2.add(m2._sire_object, pmap1)
    intraff2.add(m2._sire_object, pmap1)

    assert intraclj1.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff1.energy().value() == pytest.approx(intraff2.energy().value())

    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.

    internalff0 = InternalFF("internal")
    internalff0.setStrict(True)
    internalff0.add(m0._sire_object)

    internalff1 = InternalFF("internal")
    internalff1.setStrict(True)
    internalff1.add(m1._sire_object)

    # First extract a partial molecule using the atoms from molecule0 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for atom in m0._sire_object.atoms():
        selection.select(atom.index())
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap0)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=0 end state.
    amber_mol, _ = m2._extractMolecule()

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Now extract a partial molecule using the atoms from molecule1 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for idx in mapping.keys():
        selection.select(AtomIdx(idx))
    for idx in range(n0, m2._sire_object.nAtoms()):
        selection.select(AtomIdx(idx))
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap1)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=1 end state.
    amber_mol, _ = m2._extractMolecule(is_lambda1=True)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

@pytest.mark.xfail(reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception")
def test_ring_breaking_three_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/CAT-13a*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/CAT-17g*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Generate the mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

@pytest.mark.xfail(reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception")
def test_ring_breaking_five_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand04*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Load the pre-defined mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

@pytest.mark.xfail(reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception")
def test_ring_breaking_six_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand38*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Load the pre-defined mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

@pytest.mark.parametrize("ligands", [ pytest.param(["CAT-13c", "CAT-17i"], marks=pytest.mark.xfail(reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception")),
                                     pytest.param(["CAT-13e", "CAT-17g"], marks=pytest.mark.xfail(reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception"))
                                     ],
                        )
def test_ring_size_change(ligands):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/%s.*" % ligands[0]))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/%s.*" % ligands[1]))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Generate the mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True, allow_ring_size_change=True)

# Parameterise the function with a valid mapping.
@pytest.mark.parametrize("ligands, mapping", [(("grow1", "grow2"),
                                               {2 : 21,
                                                4 : 23,
                                                6 : 25,
                                                8 : 27,
                                                10 : 18,
                                                1 : 19,
                                                0 : 20,
                                                11 : 16,
                                                12 : 17,
                                                13 : 14,
                                                15 : 13,
                                                18 : 11,
                                                20 : 9,
                                                22 : 8,
                                                23 : 5,
                                                16 : 6,
                                                24 : 3,
                                                26 : 1,
                                                27 : 0,
                                                9 : 28,
                                                5 : 24,
                                                3 : 22,
                                                7 : 26,
                                                14 : 15,
                                                19 : 12,
                                                21 : 10,
                                                17 : 7,
                                                25 : 4}),
                                               (("grow3", "grow4"),
                                               {1: 6,
                                                2: 7,
                                                3: 8,
                                                4: 9,
                                                5: 10,
                                                6: 11,
                                                14: 21,
                                                13: 20,
                                                12: 19,
                                                11: 18,
                                                10: 17})])
def test_grow_whole_ring(ligands, mapping):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob(f"test/input/ligands/{ligands[0]}*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob(f"test/input/ligands/{ligands[1]}*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Check that we can merge without allowing ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping)

def test_hydrogen_mass_repartitioning():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand38*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Create the merged molecule.
    merged = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

    # Work out the initial mass of the system.
    initial_mass0 = 0
    for mass in merged._sire_object.property("mass0").toVector():
        initial_mass0 += mass.value()
    initial_mass1 = 0
    for mass in merged._sire_object.property("mass1").toVector():
        initial_mass1 += mass.value()

    # Repartition the hydrogen mass.
    merged.repartitionHydrogenMass()

    # Work out the final mass of the system.
    final_mass0 = 0
    for mass in merged._sire_object.property("mass0").toVector():
        final_mass0 += mass.value()
    final_mass1 = 0
    for mass in merged._sire_object.property("mass1").toVector():
        final_mass1 += mass.value()

    # Assert the the masses are approximately the same.
    assert final_mass0 == pytest.approx(initial_mass0)
    assert final_mass1 == pytest.approx(initial_mass1)

if __name__ == '__main__':
    test_roi_merge()

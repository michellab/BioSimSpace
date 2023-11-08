import BioSimSpace.Sandpit.Exscientia as BSS
from tests.conftest import root_fp


def test_sire_properties():
    s = BSS.IO.readMolecules([f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"])

    m = s[0]._sire_object

    charges = m.property("charge")

    c = charges.toVector()

    assert len(c) == m.nAtoms()

    masses = m.property("mass")

    v = masses.toVector()

    assert len(v) == m.nAtoms()

    a = m.property("ambertype")

    a = a.toVector()

    assert len(a) == m.nAtoms()

    g = m.property("gb_radii")

    g = g.toVector()

    assert len(g) == m.nAtoms()

    s = m.property("gb_screening")

    s = s.toVector()

    assert len(s) == m.nAtoms()

    t = m.property("treechain")

    t = t.toVector()

    assert len(t) == m.nAtoms()


if __name__ == "__main__":
    test_sire_properties()

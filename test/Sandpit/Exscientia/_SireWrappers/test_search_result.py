import BioSimSpace.Sandpit.Exscientia as BSS

import pytest


@pytest.fixture
def system(scope="session"):
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules("test/input/amber/ala/*")


@pytest.fixture
def molecule(system, scope="session"):
    """Re-use the same molecule for each test."""
    return system[0]


@pytest.fixture
def residue(molecule, scope="session"):
    """Re-use the same residue for each test."""
    return molecule.getResidues()[0]


# Parameterise the function with a set of searchable objects.
@pytest.mark.parametrize("fixture", ["system", "molecule", "residue"])
# Ensure that valid searches return a result.
def test_valid_search(fixture, request):
    obj = request.getfixturevalue(fixture)
    obj.search("element H")


# Parameterise the function with a set of searchable objects.
@pytest.mark.parametrize("fixture", ["system", "molecule", "residue"])
# Ensure that invalid searches raise a ValueError exception.
def test_invalid_search(fixture, request):
    obj = request.getfixturevalue(fixture)

    # Assert that a search for an invalid element raises a ValueError.
    with pytest.raises(ValueError):
        obj.search("element Z")

    # Assert that a search for a valid element using an invalid property
    # map raises a ValueError.
    with pytest.raises(ValueError):
        obj.search("element H", property_map={"element": "cheese"})

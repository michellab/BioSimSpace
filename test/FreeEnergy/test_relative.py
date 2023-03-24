import BioSimSpace as BSS

import pytest

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def perturbable_system():
    """Re-use the same perturbable system for each test."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/perturbable_system0.prm7",
        f"{url}/perturbable_system0.rst7",
        f"{url}/perturbable_system1.prm7",
        f"{url}/perturbable_system1.rst7",
    )


@pytest.mark.parametrize("engine", BSS.FreeEnergy.engines())
def test_setup(perturbable_system, engine):
    """Test setup for a relative alchemical free energy leg."""

    # Try to instantiate a free energy object for the simulation leg.
    free_nrg = BSS.FreeEnergy.Relative(perturbable_system)

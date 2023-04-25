import BioSimSpace as BSS

import pytest

# Store the tutorial URL.
url = BSS.tutorialUrl()

# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None


@pytest.fixture(scope="session")
def perturbable_system():
    """Re-use the same perturbable system for each test."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/perturbable_system0.prm7",
        f"{url}/perturbable_system0.rst7",
        f"{url}/perturbable_system1.prm7",
        f"{url}/perturbable_system1.rst7",
    )


def test_setup_somd(perturbable_system):
    """Test setup for a relative alchemical free energy leg using SOMD."""
    free_nrg = BSS.FreeEnergy.Relative(perturbable_system, engine="somd")


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_setup_gromacs(perturbable_system):
    """Test setup for a relative alchemical free energy leg using GROMACS."""
    free_nrg = BSS.FreeEnergy.Relative(perturbable_system, engine="gromacs")

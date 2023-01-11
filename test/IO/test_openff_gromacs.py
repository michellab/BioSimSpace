import BioSimSpace as BSS

import os
import pytest

# Check whether AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

# Check whether GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False,
    reason="Requires that both AMBER and GROMACS are installed.",
)
def test_molecule_combine():
    """Single point energy comparison to make sure that GROMACS
    can correctly combine Open Force Field molecules where
    duplicate atom types have different molecular parameters.
    """

    # Create two test molecules, parameterised using OpenFF.
    # These molecules will contain duplicate atom types, e.g. H1,
    # which have different parameters.

    # Toluene.
    m0 = BSS.Parameters.openff_unconstrained_2_0_0("Cc1ccccc1").getMolecule()

    # Benzene.
    m1 = BSS.Parameters.openff_unconstrained_2_0_0("c1ccccc1").getMolecule()

    # Translate m1 such that the molecules don't interact.
    m1.translate(3 * [50 * BSS.Units.Length.angstrom])

    # Create a single-step minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=1)

    # Create processes to single point calculations with GROMACS.
    p0 = BSS.Process.Gromacs(m0.toSystem(), protocol)
    p1 = BSS.Process.Gromacs(m1.toSystem(), protocol)
    p01 = BSS.Process.Gromacs((m0 + m1).toSystem(), protocol)

    # Modify the GROMACS configuration to run zero steps.
    config = p0.getConfig()
    config[1] = "nsteps = 0"
    p0.setConfig(config)
    p1.setConfig(config)
    p01.setConfig(config)

    # Run the processes and wait for them to finish.
    p0.start()
    p0.wait()
    p1.start()
    p1.wait()
    p01.start()
    p01.wait()

    # Get the potential energies.
    e0 = p0.getPotentialEnergy()
    e1 = p1.getPotentialEnergy()
    e01 = p01.getPotentialEnergy()

    # Make sure the energies are consistent.
    assert (e0 + e1).value() == pytest.approx(e01.value(), rel=1e-3)

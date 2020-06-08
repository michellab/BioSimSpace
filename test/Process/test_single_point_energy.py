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

@pytest.mark.skipif(has_amber is False or has_gromacs is False,
    reason="Requires that both AMBER and GROMACS are installed.")
def test_amber_gromacs():
    """Single point energy comparison between AMBER and GROMACS."""

    # Load the vacuum ubiquitin system.
    files = BSS.IO.glob("test/io/amber/ubiquitin/*")
    system = BSS.IO.readMolecules(files)

    # Create a single-step minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=1)

    # Create a process to run with AMBER.
    process_amb = BSS.Process.Amber(system, protocol)

    # Create a process to run with GROMACS.
    process_gmx = BSS.Process.Gromacs(system, protocol)

    # Modify the GROMACS configuration to run zero steps.
    config = process_gmx.getConfig()
    config[1] = "nsteps = 0"
    process_gmx.setConfig(config)

    # Run the AMBER process and wait for it to finish.
    process_amb.start()
    process_amb.wait()

    # Run the GROMACS process and wait for it to finish.
    process_gmx.start()
    process_gmx.wait()

    # Compare bond energies. (In kJ / mol)
    nrg_amb = process_amb.getBondEnergy().kj_per_mol().magnitude()
    nrg_gmx = process_gmx.getBondEnergy().kj_per_mol().magnitude()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

    # Compare angle energies. (In kJ / mol)
    nrg_amb = process_amb.getAngleEnergy().kj_per_mol().magnitude()
    nrg_gmx = process_gmx.getAngleEnergy().kj_per_mol().magnitude()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

    # Compare dihedral energies. (In kJ / mol)
    nrg_amb = process_amb.getDihedralEnergy().kj_per_mol().magnitude()
    nrg_gmx = process_gmx.getDihedralEnergy().kj_per_mol().magnitude()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

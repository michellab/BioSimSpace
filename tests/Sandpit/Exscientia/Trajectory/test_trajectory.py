import pytest

import BioSimSpace.Sandpit.Exscientia as BSS

try:
    from Sire.Base import wrap
except Exception:

    def wrap(arg):
        return arg


from tests.Sandpit.Exscientia.conftest import has_mdanalysis, has_mdtraj
from tests.conftest import root_fp


@pytest.fixture(scope="session")
def system():
    """A system object with the same topology as the trajectories."""
    return BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
    )


@pytest.fixture(scope="session")
def traj_sire(system):
    """A trajectory object using the Sire backend."""
    return BSS.Trajectory.Trajectory(
        trajectory=f"{root_fp}/input/ala.trr",
        topology=f"{root_fp}/input/ala.gro",
        system=system,
        backend="SIRE",
    )


@pytest.fixture(scope="session")
def traj_mdtraj(system):
    """A trajectory object using the MDTraj backend."""
    return BSS.Trajectory.Trajectory(
        trajectory=f"{root_fp}/input/ala.trr",
        topology=f"{root_fp}/input/ala.gro",
        system=system,
        backend="MDTRAJ",
    )


@pytest.fixture(scope="session")
def traj_mdanalysis(system):
    """A trajectory object using the MDAnalysis backend."""
    return BSS.Trajectory.Trajectory(
        trajectory=f"{root_fp}/input/ala.trr",
        topology=f"{root_fp}/input/ala.tpr",
        system=system,
        backend="MDANALYSIS",
    )


@pytest.fixture(scope="session")
def traj_mdanalysis_pdb(system):
    """A trajectory object using the MDAnalysis backend."""
    new_system = system.copy()
    new_system._sire_object.setProperty("fileformat", wrap("PDB"))
    return BSS.Trajectory.Trajectory(
        trajectory=f"{root_fp}/input/ala.trr",
        topology=f"{root_fp}/input/ala.tpr",
        system=new_system,
        backend="MDANALYSIS",
    )


@pytest.mark.skipif(
    has_mdanalysis is False or has_mdtraj is False,
    reason="Requires MDAnalysis and mdtraj to be installed.",
)
def test_frames(traj_sire, traj_mdtraj, traj_mdanalysis):
    """Make sure that the number of frames loaded by each backend agree."""
    assert traj_sire.nFrames() == traj_mdtraj.nFrames() == traj_mdanalysis.nFrames()


@pytest.mark.skipif(
    has_mdanalysis is False or has_mdtraj is False,
    reason="Requires MDAnalysis and mdtraj to be installed.",
)
def test_coords(traj_sire, traj_mdtraj, traj_mdanalysis):
    """Make sure that frames from both backends have comparable coordinates."""

    # Extract the first and last frame from each trajectory.
    frames0 = traj_sire.getFrames([0, -1])
    frames1 = traj_mdtraj.getFrames([0, -1])
    frames2 = traj_mdanalysis.getFrames([0, -1])

    # Make sure that all coordinates are approximately the same.
    for system0, system1, system2 in zip(frames0, frames1, frames2):
        for mol0, mol1, mol2 in zip(system0, system1, system2):
            for c0, c1, c2 in zip(
                mol0.coordinates(), mol1.coordinates(), mol2.coordinates()
            ):
                assert (
                    c0.x().value()
                    == pytest.approx(c1.x().value(), abs=1e-2)
                    == pytest.approx(c2.x().value(), abs=1e-2)
                )
                assert (
                    c0.y().value()
                    == pytest.approx(c1.y().value(), abs=1e-2)
                    == pytest.approx(c2.y().value(), abs=1e-2)
                )
                assert (
                    c0.z().value()
                    == pytest.approx(c1.z().value(), abs=1e-2)
                    == pytest.approx(c2.z().value(), abs=1e-2)
                )


@pytest.mark.skipif(
    has_mdanalysis is False or has_mdtraj is False,
    reason="Requires MDAnalysis and mdtraj to be installed.",
)
def test_coords_pdb(traj_mdtraj, traj_mdanalysis_pdb):
    """Make sure that frames from both backends have comparable coordinates
    when a PDB intermediate topology is used for reconstruction. This isn't
    needed for Sire.
    """

    # Extract the first and last frame from each trajectory.
    frames0 = traj_mdtraj.getFrames([0, -1])
    frames1 = traj_mdanalysis_pdb.getFrames([0, -1])

    # Make sure that all coordinates are approximately the same.
    for system0, system1 in zip(frames0, frames1):
        for mol0, mol1 in zip(system0, system1):
            for c0, c1 in zip(mol0.coordinates(), mol1.coordinates()):
                assert c0.x().value() == pytest.approx(c1.x().value(), abs=1e-2)
                assert c0.y().value() == pytest.approx(c1.y().value(), abs=1e-2)
                assert c0.z().value() == pytest.approx(c1.z().value(), abs=1e-2)


@pytest.mark.skipif(
    has_mdanalysis is False, reason="Requires MDAnalysis to be installed."
)
def test_velocities(traj_mdanalysis):
    """Make sure that the MDAnalysis format trajectory contains velocities."""

    # Extract the first and last frame from the trajectory.
    frames = traj_mdanalysis.getFrames([0, -1])

    # Make sure each molecule in each frame has a "velocity" property.
    for frame in frames:
        for mol in frame:
            assert mol._sire_object.hasProperty("velocity")


@pytest.mark.skipif(
    has_mdanalysis is False or has_mdtraj is False,
    reason="Requires MDAnalysis and mdtraj to be installed.",
)
def test_rmsd(traj_mdtraj, traj_mdanalysis):
    """Make sure that the RMSD computed by both backends is comparable."""

    # Compute the RMSD for a subset of atoms using the third frame
    # as a reference.
    rmsd0 = traj_mdtraj.rmsd(frame=3, atoms=[0, 10, 20, 30, 40])
    rmsd1 = traj_mdanalysis.rmsd(frame=3, atoms=[0, 10, 20, 30, 40])

    # Make sure the values are approximately the same.
    for v0, v1 in zip(rmsd0, rmsd1):
        assert v0.value() == pytest.approx(v1.value(), abs=1e-2)

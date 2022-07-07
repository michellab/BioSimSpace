import BioSimSpace.Sandpit.Exscientia as BSS

import pytest

@pytest.fixture
def traj_mdtraj(scope="session"):
    """A trajectory object using the MDTraj backend."""
    return BSS.Trajectory.Trajectory(
            trajectory="test/Sandpit/Exscientia/input/trajectories/ala.trr",
            topology="test/Sandpit/Exscientia/input/trajectories/ala.gro")

@pytest.fixture
def traj_mdanalysis(scope="session"):
    """A trajectory object using the MDAnalysis backend."""
    return BSS.Trajectory.Trajectory(
            trajectory="test/input/trajectories/ala.trr",
            topology="test/input/trajectories/ala.tpr")

def test_frames(traj_mdtraj, traj_mdanalysis):
    """Make sure that the number of frames loaded by each backend agree."""
    assert traj_mdtraj.nFrames() == traj_mdanalysis.nFrames()

def test_coords(traj_mdtraj, traj_mdanalysis):
    """Make sure that frames from both backends have comparable coordinates."""

    # Extract the first and last frame from each trajectory.
    frames0 = traj_mdtraj.getFrames([0, -1])
    frames1 = traj_mdanalysis.getFrames([0, -1])

    # Make sure that all coordinates are approximately the same.
    for system0, system1 in zip(frames0, frames1):
        for mol0, mol1 in zip(system0, system1):
            for c0, c1 in zip(mol0.coordinates(), mol1.coordinates()):
                assert c0.x().value() == pytest.approx(c1.x().value(), abs=1e-3)
                assert c0.y().value() == pytest.approx(c1.y().value(), abs=1e-3)
                assert c0.z().value() == pytest.approx(c1.z().value(), abs=1e-3)

def test_rmsd(traj_mdtraj, traj_mdanalysis):
    """Make sure that the RMSD computed by both backends is comparable."""

    # Compute the RMSD for a subset of atoms using the third frame
    # as a reference.
    rmsd0 = traj_mdtraj.rmsd(frame=3, atoms=[0, 10, 20, 30, 40])
    rmsd1 = traj_mdanalysis.rmsd(frame=3, atoms=[0, 10, 20, 30, 40])

    # Make sure the values are approximately the same.
    for v0, v1 in zip(rmsd0, rmsd1):
        assert v0.value() == pytest.approx(v1.value(), abs=1e-3)

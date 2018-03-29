# BioSimSpace.Trajectory

This sub-package provides functionality for reading and manipulating molecular
trajectories. A trajectory can be created by attaching to a running
[Process](../Process), or by passing an appropriate trajectory and topology
file. Internally, the `Trajectory` class makes use of [MDTraj](http://mdtraj.org/1.9.0)
and [MDAnalysis](https://www.mdanalysis.org) and can return trajectory
objects in either of these formats.

Some examples:

```python
import BioSimSpace as BSS

# Create a trajectory object from file.
trajectory = BSS.Trajectory(trajectory="ala.dcd", topology="ala.pdb")

# Create a trajectory object by attaching to a running process.
# (Here process is a BioSimSpace.Process object.)
trajectory = BSS.Trajectory(process=process)

# Get the latest trajectory in MDTraj format.
traj = trajectory.getTrajectory()

# Get the latest trajectory in MDAnalysis format.
traj = trajectory.getTrajectory(format="mdanalysis")

# Print the number of frames in the trajectory.
print(trajectory.nFrames())

# Get the first, second, and third frames.
frames = trajectory.getFrames([0, 1 ,2])

# Get every other frame.
frames = trajectory.getFrames([x for x in range(0, trajectory.nFrames(), 2)])
``

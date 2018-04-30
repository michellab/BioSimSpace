# BioSimSpace.Trajectory

This package provides functionality for reading and manipulating molecular
trajectories. A trajectory can be created by attaching to a running
[Process](../Process), or by passing an appropriate trajectory and topology
file. Internally, the `Trajectory` class makes use of [MDTraj](http://mdtraj.org/1.9.0)
and [MDAnalysis](https://www.mdanalysis.org) and can return trajectory
objects in either of these formats.

For particular trajectory formats, MDTraj trajectory relies on the
[NetCDF](https://en.wikipedia.org/wiki/NetCDF) format and the reader/writer
that comes as part of the [SciPy](https://www.scipy.org) package. Following
the update to NetCDF-5, the SciPy writer no longer conforms to the padding
requirements of the NetCDF standard. (See [here](https://github.com/Unidata/netcdf-c/issues/657)
for details.) This has been fixed in the development version of SciPy and
will be part of the next release. In the mean time, it is possible to patch
the installed version of SciPy using the following command:

```bash
curl https://raw.githubusercontent.com/scipy/scipy/master/scipy/io/netcdf.py \
    > $HOME/sire.app/lib/python3.5/site-packages/scipy/io/netcdf.py
```

Some examples of how to use the Trajectory package:

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
```

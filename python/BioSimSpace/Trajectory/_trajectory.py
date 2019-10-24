######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for reading and analysing molecular trajectories.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["getFrame", "Trajectory"]

import MDAnalysis as _mdanalysis
import mdtraj as _mdtraj
import os as _os
import shutil as _shutil
import warnings as _warnings

from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace.Process._process import Process as _Process
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace.Types import Time as _Time

from BioSimSpace import IO as _IO
from BioSimSpace import _SireWrappers as _SireWrappers

# A dictionary mapping the Sire file format extension to those expected by MDTraj.
_extensions = { "Gro87" : "gro",
                "PRM7"   : "parm7" }

def getFrame(trajectory, topology, index):
    """Extract a single frame from a trajectory file.

       Parameters
       ----------

       trajectory : str
           A trajectory file.

       topology : str
           A topology file.

       index : int
          The index of the frame.

       Returns
       -------

       frame : :class:`System <BioSimSpace._SireWrappers.System>`
           The System object of the corresponding frame.
    """

    if type(trajectory) is not str:
        raise TypeError("'trajectory' must be of type 'str'")

    if type(topology) is not str:
        raise TypeError("'topology' must be of type 'str'")

    if type(index) is not int:
        raise TypeError("'index' must be of type 'int'")

    # Try to load the frame.
    try:
        frame = _mdtraj.load_frame(trajectory, index, top=topology)
    except:
        # Get the file format of the topology file.
        try:
            # Load the topology file to determine the file format.
            file_format = _IO.readMolecules(topology).fileFormat()

            # Set the extension.
            extension = _extensions.get(file_format, file_format.lower())

            # Set the path to the temporary topology file.
            top_file = _os.getcwd() + "/.topology." + extension

            # Copy the topology to a file with the correct extension.
            _shutil.copyfile(topology, top_file)

            frame = _mdtraj.load_frame(trajectory, index, top=top_file)
        except:
            _os.remove(top_file)
            raise IOError("MDTraj failed to read frame %d from: traj=%s, top=%s" % (index, trajectory, topology))

        # Remove the temporary topology file.
        _os.remove(top_file)

    # The name of the frame coordinate file.
    frame_file = ".frame.nc"

    # Save the coordinates to file.
    frame.save(frame_file)

    # Load the frame into a System object.
    try:
        system = _System(_SireIO.MoleculeParser.read([topology, frame_file]))
    except Exception as e:
        _os.remove(frame_file)
        msg = "Failed to read trajectory frame: '%s'" % frame_file
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None

    # Remove the temporary frame coordinate file.
    _os.remove(frame_file)

    # Return the system.
    return system

class Trajectory():
    """A class for reading a manipulating biomolecular trajectories."""

    def __init__(self, process=None, trajectory=None, topology=None):
        """Constructor.

           Parameters
           ----------

           process : :class:`Process <BioSimSpace.Process>`
               A BioSimSpace process object.

           trajectory : str
               A trajectory file.

           topology : str
               A topology file.
        """

        # Set default member variables.
        self._process = None
        self._process_name = None
        self._traj_file = None
        self._top_file = None

        # Nothing to create a trajectory from.
        if process is None and trajectory is None:
            raise ValueError("Both 'process' and 'trajectory' keyword arguments are 'None'")

        # Both use cases active. Default to process.
        if not process is None and not trajectory is None:
            _warnings.warn("Both a process and trajectory file are specified! Defaulting to 'process'.")
            self._traj_file = None

        # BioSimSpace process.
        if process is not None:
            if isinstance(process, _Process):
                self._process = process
                self._process_name = process.__class__.__name__
                self._top_file = self._process._top_file

                # Check that the process can generate a trajectory.
                if not self._process._has_trajectory:
                    raise ValueError("BioSimSpace.Process.%s cannot generate a trajectory!" % self._process_name)

        # Trajectory and topology files.
        elif type(trajectory) is str and type(topology) is str:

            # Make sure the trajectory file exists.
            if not _os.path.isfile(trajectory):
                raise IOError("Trajectory file doesn't exist: '%s'" % trajectory)

            # Make sure the topology file exists.
            if not _os.path.isfile(topology):
                raise IOError("Topology file doesn't exist: '%s'" % topology)

            self._traj_file = trajectory
            self._top_file = topology

        # Invalid arguments.
        else:
            raise ValueError("BioSimSpace.Trajectory requires a BioSimSpace.Process object, "
                             "or a trajectory and topology file.")

        # Get the current trajectory.
        self._trajectory = self.getTrajectory()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Trajectory: nFrames=%d>" % self.nFrames()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Trajectory: nFrames=%d>" % self.nFrames()

    def getTrajectory(self, format="mdtraj"):
        """Get the current trajectory object.

           Parameters
           ----------

           format : str
               Whether to return an 'MDTraj' or 'MDAnalysis' object.

           Returns
           -------

           trajectory : mdtraj.core.trajectory.Trajectory, MDAnalysis.core.universe.Universe
               The trajectory in MDTraj or MDAnalysis format.
        """

        if format.upper() not in ["MDTRAJ", "MDANALYSIS"]:
            _warnings.warn("Invalid trajectory format. Using default (mdtraj).")
            format = "mdtraj"

        # Set the location of the trajectory and topology files.
        if self._process is not None:
            traj_file = self._process._traj_file

            # Weirdly, the GRO file is used as the topology.
            if self._process_name.upper() == "GROMACS":
                top_file = self._process._gro_file
            else:
                top_file = self._process._top_file
        else:
            traj_file = self._traj_file
            top_file = self._top_file

        # Check that the trajectory and topology files exist.
        if not _os.path.isfile(traj_file):
            raise IOError("Trajectory file doesn't exist: '%s'" % traj_file)

        if not _os.path.isfile(top_file):
            raise IOError("Topology file doesn't exist: '%s'" % top_file)

        # Load the topology file to determine the file format.
        file_format = _IO.readMolecules(top_file).fileFormat()

        # Set the extension.
        extension = _extensions.get(file_format, file_format.lower())

        # Set the path to the temporary topology file.
        new_top_file = _os.getcwd() + "/.topology." + extension

        # Copy the topology to a file with the correct extension.
        _shutil.copyfile(top_file, new_top_file)

        # Return an MDTraj object.
        if format == "mdtraj":

            try:
                traj = _mdtraj.load(traj_file, top=new_top_file)
            except:
                _warnings.warn("MDTraj failed to read: traj=%s, top=%s" % (traj_file, top_file))
                traj = None

            # Remove the temporary topology file.
            _os.remove(new_top_file)

            return traj

        # Return an MDAnalysis Universe.
        else:
            try:
                universe = _mdanalysis.Universe(new_top_file, traj_file)
            except:
                _warnings.warn("MDAnalysis failed to read: traj=%s, top=%s" % (traj_file, top_file))
                universe = None

            # Remove the temporary topology file.
            _os.remove(new_top_file)

            return universe

    def getFrames(self, indices=None):
        """Get trajectory frames as a list of System objects.

           Parameters
           ----------

           indices : [int], [:class:`Time <BioSimSpace.Types.Time>`]
               A list of trajectory frame indices, or time stamps (in ns).

           Returns
           -------

           frames : [:class:`System <BioSimSpace._SireWrappers.System>`]
               The list of System objects.
        """

        # The process is running. Grab the latest trajectory.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

            # There is no trajectory.
            if self._trajectory is None:
                return None

        # Store the number of frames.
        n_frames = self._trajectory.n_frames

        # Work out the frame spacing in nanoseconds.
        # TODO:
        # How can we do this in a robust way if the trajectory is loaded from file?
        # Some formats do not store time information as part of the trajectory.
        if n_frames > 1:
            if self._process is not None:
                time_interval = self._process._protocol.getRunTime() / self._process._protocol.getFrames()
            else:
                time_interval = self._trajectory.timestep / 1000

        # Create the indices array.

        # Default to all frames.
        if indices is None:
            indices = [x for x in range(0, self._trajectory.n_frames)]

        # A single frame index.
        elif type(indices) is int:
            indices = [indices]

        # A single time stamp.
        elif type(indices) is _Time:
            if n_frames > 1:
                # Round time stamp to nearest frame index.
                indices = [round(indices.nanoseconds().magnitude() / time_interval) - 1]
            else:
                raise _IncompatibleError("Cannot determine time stamps for a trajectory "
                                         "with only one frame!")

        # A list of frame indices.
        elif all(isinstance(x, int) for x in indices):
            pass

        # A list of time stamps.
        elif all(isinstance(x, _Time) for x in indices):
            if n_frames <= 1:
                raise _IncompatibleError("Cannot determine time stamps for a trajectory "
                                         "with only one frame!")

            # Round time stamps to nearest frame indices.
            indices = [round(x.nanoseconds().magnitude() / time_interval) - 1 for x in indices]

        # Unsupported argument.
        else:
            raise ValueError("Unsupported argument. Indices or time stamps "
                             "must be an 'int' or 'BioSimSpace.Types.Time', or list of 'int' or "
                             "'BioSimSpace.Types.Time' types.")

        # Intialise the list of frames.
        frames = []

        # Loop over all indices.
        for x in indices:

            # Make sure the frame index is within range.
            if x > 0 and x >= n_frames:
                raise ValueError("Frame index (%d) of of range (0 to %d)." % (x, n_frames - 1))
            elif x < -n_frames:
                raise ValueError("Frame index (%d) of of range (-1 to -%d)." % (x, n_frames))

            # The name of the frame coordinate file.
            frame_file = ".frame.nc"

            # Write the current frame as a NetCDF file.
            self._trajectory[x].save(frame_file)

            # Load the frame and create a System object.
            try:
                system = _System(_SireIO.MoleculeParser.read([self._top_file, frame_file]))
            except Exception as e:
                _os.remove(frame_file)
                msg = "Failed to read trajectory frame: '%s'" % frame_file
                if _isVerbose():
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None

            # Append the system to the list of frames.
            frames.append(system)

        # Remove the temporary frame coordinate file.
        _os.remove(frame_file)

        # Return the frames.
        return frames

    def nFrames(self):
        """Return the current number of trajectory frames.

           Returns
           -------

           nFrames : int
               The number of trajectory frames.
        """

        # First get the current MDTraj object.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

        # There is no trajectory.
        if self._trajectory is None:
            return 0
        else:
            return self._trajectory.n_frames

    def rmsd(self, frame=None, atoms=None):
        """Compute the root mean squared displacement.

           Parameters
           ----------

           frame : int
               The index of the reference frame.

           atoms : [int]
               A list of reference atom indices.

           Returns
           -------

           rmsd : [float]
               A list containing the RMSD value at each time point.
        """

        # Default to the first frame.
        if frame is None:
            frame = 0
        else:
            if type(frame) is not int:
                raise TypeError("'frame' must be of type 'int'")
            else:
                # Store the number of frames.
                n_frames = self._trajectory.n_frames

                # Make sure the frame index is within range.
                if frame > 0 and frame >= n_frames:
                    raise ValueError("Frame index (%d) of of range (0 to %d)." %s (frame, n_frames - 1))
                elif frame < -n_frames:
                    raise ValueError("Frame index (%d) of of range (-1 to -%d)." %s (frame, n_frames))

        if atoms is not None:
            # Check that all of the atom indices are integers.
            if not all(isinstance(x, int) for x in atoms):
                raise TypeError("'atom' indices must be of type 'int'")

        # Use MDTraj to compute the RMSD.
        try:
            rmsd = _mdtraj.rmsd(self._trajectory, self._trajectory, frame, atoms)
        except Exception as e:
            msg = "Atom indices not found in the system."
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        # Convert to a list and return.
        return list(rmsd)

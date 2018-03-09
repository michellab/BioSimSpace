"""
@package biosimspace
@author  Lester Hedges
@brief   A class for reading and manipulating biomolecular trajectories.
"""

from ..Process.process import Process

from Sire import try_import

from os import path
from shutil import copyfile
from warnings import warn

try:
    mdtraj = try_import("mdtraj")
except ImportError:
    raise ImportError("MDTraj is not installed. Please install mdtraj in order to use BioSimSpace.")

try:
    mdanalysis = try_import("MDAnalysis")
except ImportError:
    raise ImportError("MDAnalysis is not installed. Please install mdanalysis in order to use BioSimSpace.")

class Trajectory():
    """A class for reading a manipulating biomolecular trajectories."""

    def __init__(self, handle):
        """Constructor.

           Positional arguments:

           handle -- A handle to a BioSimSpace.Process object or a list of trajectory files.
        """

        # Set default member variables.
        self._process = None
        self._process_name = None
        self._traj_files = None

        # BioSimSpace process.
        if handle.__class__.__base__ == Process:
            self._process = handle
            self._process_name = handle.__class__.__name__

            # Check that the process can generate a trajectory.
            if not self._process._has_trajectory:
                raise ValueError("BioSimSpace.Process.%s cannot generate a trajectory!" % self._process_name)

        # List of trajectory files.
        elif all(isinstance(x, str) for x in handle):
            self._traj_files = handle

        # Single trajectory file.
        elif type(handle) is str:
            self._traj_files = [handle]

        # Invalid arguments.
        else:
            raise ValueError("The handle must be a BioSimSpace.Process object, "
                "a trajectory file, or a list of trajectory files.")

        # Make sure the trajectory files exist.
        if self._traj_files is not None:
            for file in self._traj_files:
                if not path.isfile(file):
                    raise IOError(('Trajectory file doesn\'t exist: "{x}"').format(x=file))

        # Get the current trajectory.
        self._trajectory = self.getTrajectory()

    def getTrajectory(self, format='mdtraj'):
        """Get the current trajectory object.

           Keyword arguments:

           format -- Whether to return a 'MDTraj' or 'MDAnalysis' object.
        """

        if format.upper() not in ['MDTRAJ', 'MDANALYSIS']:
            warn("Invalid trajectory format. Using default (mdtraj).")
            format = mdtraj

        # Get the trajectory from the process.
        if self._process is not None:

            # AMBER.
            if self._process_name.upper() == 'AMBER':

                # Path to the trajectory file.
                traj_file = "%s/%s.nc" % (self._process._work_dir, self._process._name)

                # MDTraj currently doesn't support the .prm7 extension, so we
                # need to copy the topology file to a temporary .parm7 file.
                # I've submitted a pull request to add support for the alternative
                # .prm7 extension.
                top_file = "%s/tmp.parm7" % self._process._work_dir
                copyfile(self._process._prm_file, top_file)

                # Set the topology format.
                top_format = 'PARM7'

            # NAMD.
            elif self._process_name.upper() == 'NAMD':

                # Path to the trajectory and topology files.
                traj_file = "%s/%s_out.dcd" % (self._process._work_dir, self._process._name)
                top_file = "%s/%s.pdb" % (self._process._work_dir, self._process._name)

                # Set the topology format.
                top_format = 'PDB'

            # Unsupported process.
            else:
                raise ValueError("BioSimSpace.Process.%s is unsupported!" % self._process_name)

            # Check that the trajectory and topology files exist.
            if not path.isfile(traj_file):
                raise IOError(('Trajectory file doesn\'t exist: "{x}"').format(x=traj_file))

            if not path.isfile(top_file):
                raise IOError(('Topology file doesn\'t exist: "{x}"').format(x=top_file))

            # Return an MDTraj object.
            if format is 'mdtraj':

                # Create the MDTraj object.
                traj = mdtraj.load(traj_file, top=top_file)

                # Delete the temporary .parm7 file.
                if self._process_name.upper() is 'AMBER':
                    remove(top_file)

                return traj

            # Return an MDAnalysis Universe.
            else:
                universe = mdanalysis.Universe(top_file, traj_file, topology_format=top_format)

                # Delete the temporary .parm7 file.
                if self._process_name.upper() is 'AMBER':
                    remove(top_file)

                return universe


        # TODO:
        # Reading trajectory from file isn't current supported.
        else:
            return None

    def getFrames(self, indices=None):
        """Get trajectory frames as a list of Sire systems.

           Keyword arguments:

           indices -- A list of trajectory frame indices, or time stamps (in ns).
        """

        # First get the current MDTraj object.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

        # TODO:
        # Reading trajectory from file isn't current supported.
        else:
            return None

        # There is no trajectory.
        if traj is None:
            return None

        # Work out the frame spacing in nanoseconds.
        time_interval = self._process._protocol.runtime / self._process._protocol.frames

        # Create the indices array.

        # Default to all frames.
        if indices is None:
            indices = [x for x in range(0, traj.n_frames)]

        # A single frame index.
        elif type(indices) is int:
            indices = [indices]

        # A list of frame indices.
        elif all(isinstance(x, int) for x in indices):
            pass

        # A single time stamp.
        elif type(indices) is float:
            if indices < 0:
                raise ValueError("Time stamp cannot be negative.")

            # Round time stamp to nearest frame index.
            indices = [round(indices / time_interval) - 1]

        # A list of time stamps.
        elif all(isinstance(x, float) for x in indices):
            # Make sure no time stamps are negative.
            for x in indices:
                if x < 0:
                    raise ValueError("Time stamp cannot be negative.")

            # Round time stamps to nearest frame indices.
            indices = [round(x / time_interval) - 1 for x in indices]

        # Unsupported argument.
        else:
            raise ValueError("Unsupported argument. Indices or time stamps "
                "must be an 'int' or 'float', or list of 'int' or 'float' types.")

        # Intialise the list of frames.
        frames = []

        # Store the maximum frame number.
        max_frame = traj.n_frames - 1

        # Loop over all indices.
        for x in indices:

            # Make sure the frame index is within range.
            if abs(x) > max_frame:
                raise ValueError("Frame index (%d) of of range (0-%d)." %s (x, max_frame))

            # The name of the frame coordinate file.
            frame_file = "%s/frame.nc" % self._work_dir

            # Write the current frame as a NetCDF file.
            # MDTraj is unable to write a CRD file that is compatible with the original topology!
            traj[x].save(frame_file)

            # Create a Sire system.
            system = MoleculeParser.read([self._prm_file, frame_file])

            # Append the system to the list of frames.
            frames.append(system)

        # Remove the temporary frame coordinate file.
        remove(frame_file)

        # Return the frames.
        return frames

    def nFrames(self):
        """Return the current number of trajectory frames."""

        # First get the current MDTraj object.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

        # TODO:
        # Reading trajectory from file isn't current supported.
        else:
            return None

        # There is no trajectory.
        if traj is None:
            return 0
        else:
            return traj.n_frames

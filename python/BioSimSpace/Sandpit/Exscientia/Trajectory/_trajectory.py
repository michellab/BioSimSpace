######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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
__email__ = "lester.hedges@gmail.com"

__all__ = ["getFrame", "Trajectory"]

from .._Utils import _try_import

_mdanalysis = _try_import("MDAnalysis")
_mdtraj = _try_import("mdtraj")
import os as _os
import shutil as _shutil
import warnings as _warnings

from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from .. import _isVerbose
from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Process._process import Process as _Process
from .._SireWrappers import System as _System
from ..Types import Time as _Time

from .. import IO as _IO
from .. import Units as _Units

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

    if not isinstance(trajectory, str):
        raise TypeError("'trajectory' must be of type 'str'")

    if not isinstance(topology, str):
        raise TypeError("'topology' must be of type 'str'")

    if not type(index) is int:
        raise TypeError("'index' must be of type 'int'")

    # Try to load the frame.
    try:
        frame = _mdtraj.load_frame(trajectory, index, top=topology)
    except:
        raise IOError("MDTraj failed to read frame %d from: traj=%s, top=%s" % (index, trajectory, topology))

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
        self._backend = None

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
                if self._process_name == "Gromacs":
                    self._top_file = self._process._gro_file
                else:
                    self._top_file = self._process._top_file

                # Check that the process can generate a trajectory.
                if not self._process._has_trajectory:
                    raise ValueError("BioSimSpace.Process.%s cannot generate a trajectory!" % self._process_name)

        # Trajectory and topology files.
        elif isinstance(trajectory, str) and isinstance(topology, str):

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
        self._trajectory = self.getTrajectory(format="AUTO")

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Trajectory: nFrames=%d, backend=%r>" \
            % (self.nFrames(), str(self._backend).lower())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Trajectory: nFrames=%d, backend=%r>" \
            % (self.nFrames(), str(self._backend).lower())

    def getTrajectory(self, format="mdtraj"):
        """Get the current trajectory object.

           Parameters
           ----------

           format : str
               Whether to return an 'MDTraj' or 'MDAnalysis' object.
               Use "auto" if you are happy with either format, i.e. it
               will try each backend in sequence and return an object
               from the first one that works.

           Returns
           -------

           trajectory : mdtraj.core.trajectory.Trajectory, MDAnalysis.core.universe.Universe
               The trajectory in MDTraj or MDAnalysis format.
        """

        # Strip whitespace and convert to upper case.
        format = format.replace(" ", "").upper()

        if format.replace(" ", "").upper() not in ["MDTRAJ", "MDANALYSIS", "AUTO"]:
            _warnings.warn("Invalid trajectory format. Using default (mdtraj).")
            format = "MDTRAJ"

        # If this object isn't bound to a Process and the format matches the
        # existing backend, then return the trajectory directly.
        if self._process is None:
            if format in ["MDTRAJ", "AUTO"] and self._backend == "MDTRAJ":
                return _copy.deepcopy(self._trajectory)
            elif format in ["MDANALYSIS", "AUTO"] and self._backend == "MDANALYSIS":
                return self._trajectory.copy()

        if format == "MDTRAJ" and self._backend == "MDANALYSIS":
            raise _IncompatibleError("This Trajectory object can only be used "
                                     "with the MDAnalysis backend.")

        # Set the location of the trajectory and topology files.
        if self._process is not None:
            traj_file = self._process._traj_file

            if self._process_name.upper() == "GROMACS":
                if format == "mdtraj":
                    top_file = self._process._gro_file
                else:
                    top_file = self._process._tpr_file
                    self._backend = "MDANALYSIS"
                    format = "MDANALYSIS"
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

        # Return an MDTraj object.
        if format in ["MDTRAJ", "AUTO"]:
            try:
                traj = _mdtraj.load(traj_file, top=top_file)
                if self._backend is None:
                    self._backend = "MDTRAJ"
                return traj
            except:
                if format == "MDTRAJ":
                    _warnings.warn("MDTraj failed to read: traj=%s, top=%s" % (traj_file, top_file))
                    return None

        # Return an MDAnalysis Universe.
        if format in ["MDANALYSIS", "AUTO"]:
            # Check the file extension.
            _, extension = _os.path.splitext(top_file)

            # Whether we've created a temporary duplicate topology file.
            is_temp_file = False

            # If this is a PRM7 file, copy to PARM7.
            if extension == ".prm7":
                # Set the path to the temporary topology file.
                new_top_file = _os.getcwd() + "/.topology.parm7"

                # Copy the topology to a file with the correct extension.
                _shutil.copyfile(top_file, new_top_file)
                top_file = new_top_file
                is_temp_file = True

            try:
                universe = _mdanalysis.Universe(top_file, traj_file)
                if self._backend is None:
                    self._backend = "MDANALYSIS"

                if is_temp_file:
                    _os.remove(top_file)
            except:
                if format == "MDANALYSIS":
                    _warnings.warn("MDAnalysis failed to read: traj=%s, top=%s" % (traj_file, top_file))
                else:
                    _warnings.warn("MDTraj and MDAnalysis failed to read: traj=%s, top=%s" % (traj_file, top_file))
                universe = None

                if is_temp_file:
                    _os.remove(top_file)

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
        n_frames = self.nFrames()

        # Work out the frame spacing in nanoseconds.
        # TODO:
        # How can we do this in a robust way if the trajectory is loaded from file?
        # Some formats do not store time information as part of the trajectory.
        if n_frames > 1:
            if self._process is not None and self._process._package_name != "OPENMM":
                time_interval = (self._process._protocol.getRunTime() / self._process._protocol.getRestartInterval())
                time_interval = time_interval.nanoseconds().value()
            else:
                if self._backend == "MDTRAJ":
                    time_interval = self._trajectory.timestep / 1000
                elif self._backend == "MDANALYSIS":
                    time_interval = self._trajectory.trajectory.totaltime / 1000

        # Create the indices array.

        # Default to all frames.
        if indices is None:
            indices = [x for x in range(0, n_frames)]

        # A single frame index.
        elif type(indices) is int:
            indices = [indices]

        # A single time stamp.
        elif isinstance(indices, _Time):
            if n_frames > 1:
                # Round time stamp to nearest frame index.
                indices = [round(indices.nanoseconds().value() / time_interval) - 1]
            else:
                raise _IncompatibleError("Cannot determine time stamps for a trajectory "
                                        "with only one frame!")

        # A list of frame indices.
        elif all(type(x) is int for x in indices):
            pass

        # A list of time stamps.
        elif all(isinstance(x, _Time) for x in indices):
            if n_frames <= 1:
                raise _IncompatibleError("Cannot determine time stamps for a trajectory "
                                        "with only one frame!")

            # Round time stamps to nearest frame indices.
            indices = [round(x.nanoseconds().value() / time_interval) - 1 for x in indices]

        # Unsupported argument.
        else:
            raise ValueError("Unsupported argument. Indices or time stamps "
                            "must be an 'int' or 'BioSimSpace.Types.Time', or list of 'int' or "
                            "'BioSimSpace.Types.Time' types.")

        # Initialise the list of frames.
        frames = []

        # Sort the indices.
        indices.sort()

        # Loop over all indices.
        for x in indices:

            # Make sure the frame index is within range.
            if x > 0 and x >= n_frames:
                raise ValueError("Frame index (%d) of of range (0 to %d)." % (x, n_frames - 1))
            elif x < -n_frames:
                raise ValueError("Frame index (%d) of of range (-1 to -%d)." % (x, n_frames))

            # Write the current frame as a NetCDF file.

            if self._backend == "MDTRAJ":
                frame_file = ".frame.nc"
                self._trajectory[x].save(frame_file)
            elif self._backend == "MDANALYSIS":
                frame_file = ".frame.pdb"
                self._trajectory.trajectory[x]
                with _warnings.catch_warnings():
                    _warnings.simplefilter("ignore")
                    self._trajectory.select_atoms("all").write(frame_file)

            # Load the frame and create a System object.
            try:
                if self._backend == "MDTRAJ":
                    system = _System(_SireIO.MoleculeParser.read([self._top_file, frame_file]))
                elif self._backend == "MDANALYSIS":
                    system = _System(_SireIO.MoleculeParser.read(frame_file))
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

        # Rewind the trajectory if the backend is MDAnalysis.
        if self._backend == "MDANALYSIS":
            self._trajectory.trajectory.rewind

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
            if self._backend == "MDTRAJ":
                return self._trajectory.n_frames
            elif self._backend == "MDANALYSIS":
                return self._trajectory.trajectory.n_frames

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

           rmsd : [:class:`Length <BioSimSpace.Types.Length>`]
               A list containing the RMSD value at each time point.
        """

        # Default to the first frame.
        if frame is None:
            frame = 0
        else:
            if not type(frame) is int:
                raise TypeError("'frame' must be of type 'int'")
            else:
                # Store the number of frames.
                n_frames = self.nFrames()

                # Make sure the frame index is within range.
                if frame > 0 and frame >= n_frames:
                    raise ValueError("Frame index (%d) of of range (0 to %d)." % (frame, n_frames - 1))
                elif frame < -n_frames:
                    raise ValueError("Frame index (%d) of of range (-1 to -%d)." % (frame, n_frames))

        if atoms is not None:
            # Check that all of the atom indices are integers.
            if not all(type(x) is int for x in atoms):
                raise TypeError("'atom' indices must be of type 'int'")

        if self._backend == "MDTRAJ":
            try:
                rmsd = _mdtraj.rmsd(self._trajectory, self._trajectory, frame, atoms)
            except Exception as e:
                msg = "Atom indices not found in the system."
                if _isVerbose():
                    raise ValueError(msg) from e
                else:
                    raise ValueError(msg) from None

            # Convert to a list and add units.
            rmsd = [_Units.Length.nanometer * float(x) for x in rmsd]

        else:
            _rms = _try_import("MDAnalysis.analysis.rms")

            # Create the atom selection.
            if atoms is None:
                select = "all"
            else:
                select = "bynum " + " ".join([str(x+1) for x in atoms])

            # Instantiate the RMSD object.
            R = _rms.RMSD(self._trajectory, self._trajectory, select, ref_frame=frame)

            # Run the analysis.
            R.run()

            # Get the RMSD results.
            results = R.results.rmsd

            rmsd = [_Units.Length.nanometer * float(x[-1]/10) for x in results]

        # Return the RMSD result.
        return rmsd

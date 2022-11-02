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
import copy as _copy
import os as _os
import shutil as _shutil
import tempfile as _tempfile
import uuid as _uuid
import warnings as _warnings

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol

from .. import _isVerbose
from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Process._process import Process as _Process
from .._SireWrappers import System as _System
from ..Types import Time as _Time

from .. import IO as _IO
from .. import Units as _Units

def getFrame(trajectory, topology, index, system=None, property_map={}):
    """Extract a single frame from a trajectory file.

       Parameters
       ----------

       trajectory : str
           A trajectory file.

       topology : str
           A topology file.

       index : int
          The index of the frame.

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           A BioSimSpace System object. If passed, then trajectory coordinates
           and velocities will be copied into this system and returned. This is
           useful if you wish to preserve all molecular properties, which might
           be lost during reading of the topology format required by the
           trajectory backend. It is assumed that he system contains the same
           number of molecules as each frame in the trajectory, and that the
           molecules are in the same order.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

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

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    if system is not None:
        if not isinstance(system, _System):
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        system = system.copy()

        # We assume the molecules are in the same order, so create a one-to-one
        # molecule mapping between the frame and reference.
        mapping = {_SireMol.MolIdx(x) : _SireMol.MolIdx(x) for x in range(0, system.nMolecules())}

        # Update the water topology to match topology/trajectory.
        system = _update_water_topology(system, topology, trajectory)

    # Create a temporary working directory.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Try to load the frame with MDTraj.
    errors = []
    is_mdanalysis = False
    pdb_file = work_dir + f"/{str(_uuid.uuid4())}.pdb"
    try:
        frame_file = work_dir + f"/{str(_uuid.uuid4())}.rst7"
        frame = _mdtraj.load_frame(trajectory, index, top=topology)
        frame.save(frame_file, force_overwrite=True)
        frame.save(pdb_file, force_overwrite=True)
    except Exception as e:
        is_mdanalysis = True
        errors.append(f"MDTraj: {str(e)}")
        # Try to load the frame with MDAnalysis.
        try:
            frame_file = work_dir + f"/{str(_uuid.uuid4())}.gro"
            universe = _mdanalysis.Universe(topology, trajectory)
            universe.trajectory.trajectory[index]
            with _warnings.catch_warnings():
                _warnings.simplefilter("ignore")
                universe.select_atoms("all").write(frame_file)
                universe.select_atoms("all").write(pdb_file)
        except Exception as e:
            errors.append(f"MDAnalysis: {str(e)}")
            msg = "MDTraj/MDAnalysis failed to read frame %d from: traj=%s, top=%s" % (index, trajectory, topology)
            if _isVerbose():
                raise IOError(msg + "\n" + "\n".join(errors))
            else:
                raise IOError(msg) from None

    # Try to update the coordinates/velocities in the reference system.
    if system is not None:
        # Parse the coordinates/velocites frame.
        try:
            if is_mdanalysis:
                frame = _SireIO.Gro87(frame_file)
            else:
                frame = _SireIO.AmberRst7(frame_file)
        except Exception as e:
            msg = "Failed to read trajectory frame: '%s'" % frame_file
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Parse the PDB frame.
        try:
            pdb = _SireIO.PDB2(pdb_file)
        except Exception as e:
            msg = "Failed to read PDB trajectory frame: '%s'" % pdb_file
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # The new_system object will contain a single molecule with the
        # coordinates of all of the atoms in the reference. As such, we
        # will need to split the system into molecules.
        new_system = _split_molecules(frame, pdb, system, work_dir, property_map)
        return _System(new_system)
        try:
            sire_system, _ = _SireIO.updateCoordinatesAndVelocities(
                    system._sire_object,
                    new_system,
                    mapping,
                    False,
                    property_map,
                    {})

            new_system = _System(sire_system)
        except Exception as e:
            msg = "Failed to copy trajectory coordinates/velocities into BioSimSpace system!"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        else:
            raise IOError("The trajectory frame is incompatible with the passed system!")

    # Load the frame directly to create a new System object.
    else:
        try:
            if is_mdanalysis:
                new_system = _System(_SireIO.MoleculeParser.read(frame_file))
            else:
                new_system = _System(_SireIO.MoleculeParser.read([frame_file, topology]))
        except Exception as e:
            msg = "Failed to read trajectory frame: '%s'" % frame_file
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    # Return the system.
    return new_system

class Trajectory():
    """A class for reading a manipulating biomolecular trajectories."""

    def __init__(self, process=None, trajectory=None,
                 topology=None, system=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           process : :class:`Process <BioSimSpace.Process>`
               A BioSimSpace process object.

           trajectory : str
               A trajectory file.

           topology : str
               A topology file.

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               A BioSimSpace System object. If passed, then trajectory coordinates
               and velocities will be copied into this system and returned. This is
               useful if you wish to preserve all molecular properties, which might
               be lost during reading of the topology format required by the
               trajectory backend.

           property_map : dict
              A dictionary that maps system "properties" to their user defined
              values. This allows the user to refer to properties with their
              own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Set default member variables.
        self._process = None
        self._traj_file = None
        self._top_file = None
        self._backend = None
        self._system = None
        self._property_map = {}

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
                process_name = process.__class__.__name__
                # Check that the process can generate a trajectory.
                if not self._process._has_trajectory:
                    raise ValueError("BioSimSpace.Process.%s cannot generate a trajectory!" % process_name)
            else:
                raise TypeError("'process' must be of type 'BioSimSpace.Process'.")

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

        if system is not None:
            if not isinstance(system, _System):
                raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
            self._system = system.copy()

            # We assume the molecules are in the same order, so create a one-to-one
            # molecule mapping between the frame and reference.
            self._mapping = {_SireMol.MolIdx(x) : _SireMol.MolIdx(x) for x in range(0, system.nMolecules())}

            if process is not None:
                # If this is a GROMACS process, convert the water topology to
                # GROMACS format. Use AMBER for everything else.
                if process._package_name == "GROMACS":
                    self._system._set_water_topology("GROMACS")
                else:
                    self._system._set_water_topology("AMBER")
            else:
                # Update the water topology to match topology/trajectory.
                self._system = _update_water_topology(self._system, self._top_file, self._traj_file)

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        self._property_map = property_map

        # Create a temporary working directory.
        self._tmp_dir = _tempfile.TemporaryDirectory()
        self._work_dir = self._tmp_dir.name

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

    def getTrajectory(self, format="auto"):
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

            if self._process._package_name == "GROMACS":
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

            # If this is a PRM7 file, copy to PARM7.
            if extension == ".prm7":
                # Set the path to the temporary topology file.
                new_top_file = self._work_dir + f"/{str(_uuid.uuid4())}.parm7"

                # Copy the topology to a file with the correct extension.
                _shutil.copyfile(top_file, new_top_file)
                top_file = new_top_file

            try:
                universe = _mdanalysis.Universe(top_file, traj_file)
                if self._backend is None:
                    self._backend = "MDANALYSIS"

            except:
                if format == "MDANALYSIS":
                    _warnings.warn("MDAnalysis failed to read: traj=%s, top=%s" % (traj_file, top_file))
                else:
                    _warnings.warn("MDTraj and MDAnalysis failed to read: traj=%s, top=%s" % (traj_file, top_file))
                universe = None

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

            # Write the current frame to file.

            pdb_file = self._work_dir + f"/{str(_uuid.uuid4())}.pdb"

            if self._backend == "MDTRAJ":
                frame_file = self._work_dir + f"/{str(_uuid.uuid4())}.rst7"
                self._trajectory[x].save(frame_file, force_overwrite=True)
                self._trajectory[x].save(pdb_file, force_overwrite=True)
            elif self._backend == "MDANALYSIS":
                frame_file = self._work_dir + f"/{str(_uuid.uuid4())}.gro"
                self._trajectory.trajectory[x]
                with _warnings.catch_warnings():
                    _warnings.simplefilter("ignore")
                    self._trajectory.select_atoms("all").write(frame_file)
                    self._trajectory.select_atoms("all").write(pdb_file)

            # Try to update the coordinates/velocities in the reference system.
            if self._system is not None:
                # Parse the coordinates/velocites frame.
                try:
                    if self._backend == "MDANALYSIS":
                        frame = _SireIO.Gro87(frame_file)
                    else:
                        frame = _SireIO.AmberRst7(frame_file)
                except Exception as e:
                    msg = "Failed to read trajectory frame: '%s'" % frame_file
                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg) from None

                # Parse the PDB file.
                try:
                    pdb = _SireIO.PDB2(pdb_file)
                except Exception as e:
                    msg = "Failed to read PDB trajectory frame: '%s'" % pdb_file
                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg) from None

                # The new_system object will contain a single molecule with the
                # coordinates of all of the atoms in the reference. As such, we
                # will need to split the system into molecules.
                new_system = _split_molecules(frame, pdb, self._system, self._work_dir, self._property_map)
                try:
                    sire_system, _ = _SireIO.updateCoordinatesAndVelocities(
                            self._system._sire_object,
                            new_system,
                            self._mapping,
                            False,
                            self._property_map,
                            {})

                    new_system = _System(sire_system)
                except Exception as e:
                    msg = "Failed to copy trajectory coordinates/velocities into BioSimSpace system!"
                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg) from None

            # Load the frame directly to create a new System object.
            else:
                try:
                    if self._backend == "MDANALYSIS":
                        new_system = _System(_SireIO.MoleculeParser.read(frame_file))
                    else:
                        new_system = _System(_SireIO.MoleculeParser.read([frame_file, self._top_file]))
                except Exception as e:
                    msg = "Failed to read trajectory frame: '%s'" % frame_file
                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg) from None

            # Append the system to the list of frames.
            frames.append(new_system)

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

def _split_molecules(frame, pdb, reference, work_dir, property_map={}):
    """Internal helper function to split molecules in a "squashed" system based
       on a reference, i.e. we use the number of atoms per molecule in the
       reference to break the system apart.

       Parameters
       ----------

       frame : Sire.IO.AmberRst7, SireIO.Gro87
           The trajectory frame, parsed to AmberRst7 or Gro87 format.

       pdb : Sire.IO.PDB2
           A PDB representation of the frame.

       reference : :class:`System <BioSimSpace._SireWrappers.System>`
           A BioSimSpace System object for reference topology.

       work_dir : str
           The working directory for the trajectory.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : Sire.System.System
           The passed system, split into the appropriate number of molecules,
           as a Sire System object.
    """

    if not isinstance(frame, (_SireIO.AmberRst7, _SireIO.Gro87)):
        raise TypeError("'frame' must be of type 'Sire.IO.AmberRst7' or 'Sire.IO.Gro87'")

    if not isinstance(pdb, _SireIO.PDB2):
        raise TypeError("'pdb' must be of type 'Sire.IO.PDB2'")

    if not isinstance(reference, _System):
        raise TypeError("'reference' must be of type 'BioSimSpace._SireWrappers.System'")

    if not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    if not _os.path.exists(work_dir):
        raise ValueError(f"'work_dir' doesn't exist: {work_dir}")

    # Get the coordinates and velocity property names.
    coords = property_map.get("coordinates", "coordinates")
    vel_prop = property_map.get("velocity", "velocity")

    # Whether the frame contains velocity information.
    has_vels = frame.hasVelocities()

    # Store the formats associated with the reference system.
    formats = reference.fileFormat()

    # Write the frame coordinates/velocities to file.
    coord_file = work_dir + f"/{str(_uuid.uuid4())}.coords"
    top_file = work_dir + f"/{str(_uuid.uuid4())}.top"
    frame.writeToFile(coord_file)

    # Whether we've parsed as a PDB file.
    is_pdb = False

    if "PRM7" in formats:
        try:
            top =  _SireIO.AmberPrm(reference._sire_object)
            top.writeToFile(top_file)
        except Exception as e:
            msg = "Unable to write reference system to AmberPrm7 format!"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    elif "GroTop" in formats:
        try:
            top =  _SireIO.GroTop(reference._sire_object)
            top.writeToFile(top_file)
        except Exception as e:
            msg = "Unable to write reference system to GroTop format!"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    elif "PSF" in formats:
        try:
            top =  _SireIO.CharmmPSF(reference._sire_object)
            top.writeToFile(top_file)
        except Exception as e:
            msg = "Unable to write reference system to CharmmPSF format!"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    # Unknown, try using PDB.
    else:
        is_pdb = True

        # Get the PDB records.
        pdb_lines = pdb.toLines()

        # Create a list to hold the new lines.
        new_lines = []

        # Find the first atom record.
        for idx, line in enumerate(pdb_lines):
            if line.startswith("ATOM"):
                break

        # Loop over all of the molecules in the reference, adding the records fo
        # each atom. We append a standalone TER record between each molecule,
        # which is used by Sire.IO.PDB2 to split molecules.
        for mol in reference:
            new_lines.extend(pdb_lines[idx:idx+mol.nAtoms()])
            new_lines.append("TER")
            idx += mol.nAtoms()

        # Recreate the PDB object using the updated records.
        pdb = _SireIO.PDB2(new_lines)

        # Convert to a system.
        split_system = pdb.toSystem()

    if not is_pdb:
        # Try to read the system back in, making sure that the numbering is unique.
        try:
            split_system = _SireIO.MoleculeParser.read([coord_file, top_file])

        except Exception as e:
            msg = "Unable to read trajectory frame!"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    return split_system

def _update_water_topology(system, topology, trajectory):
    """Internal helper function to update the water topology of the system
       so that it is consistent with the passed topology or trajectory.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           A BioSimSpace System object for reference topology.

       trajectory : str
           A trajectory file.

       topology : str
           A topology file.

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The passed system with updated water topology.
    """

    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    if not isinstance(trajectory, str):
        raise TypeError("'trajectory' must be of type 'str'")

    if not isinstance(topology, str):
        raise TypeError("'topology' must be of type 'str'")

    # GROMACS topology file.
    try:
        top = _SireIO.GroTop(topology)
        system._set_water_topology("GROMACS")
        matched_topology = True
    except:
        # GROMACS coordinate file.
        try:
            top = _SireIO.Gro87(topology)
            system._set_water_topology("GROMACS")
            matched_topology = True
        # Amber topology file.
        except:
            try:
                top = _SireIO.AmberPrm(topology)
                system._set_water_topology("AMBER")
                if top.toString() != "AmberPrm::null":
                    matched_topology = True
                else:
                    matched_topology = False
            except:
                matched_topology = False

    # If we couldn't determine the topology format, try to guess
    # from the extension.
    if not matched_topology:
        ext = _os.path.splitext(topology)
        if len(ext) == 2:
            ext = ext[1]
            if ext.upper() == ".TPR":
                system._set_water_topology("GROMACS")
                matched_topology = True

    # If nothing matched, default to AMBER water format.
    if not matched_topology:
        system._set_water_topology("AMBER")

    return system

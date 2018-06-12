######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from ..Process._process import Process as _Process
from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule

import MDAnalysis as _mdanalysis
import mdtraj as _mdtraj
import os as _os
import shutil as _shutil
import warnings as _warnings

__all__ = ["Trajectory"]

class Trajectory():
    """A class for reading a manipulating biomolecular trajectories."""

    def __init__(self, process=None, trajectory=None, topology=None):
        """Constructor.

           Keyword arguments:

           process    -- A BioSimSpace process object.
           trajectory -- A trajectory file.
           topology   -- A topology file.
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

    def getTrajectory(self, format="mdtraj"):
        """Get the current trajectory object.

           Keyword arguments:

           format -- Whether to return a 'MDTraj' or 'MDAnalysis' object.
        """

        if format.upper() not in ["MDTRAJ", "MDANALYSIS"]:
            _warnings.warn("Invalid trajectory format. Using default (mdtraj).")
            format = "mdtraj"

        # Get the trajectory from the process.
        if self._process is not None:

            # AMBER.
            if self._process_name.upper() == "AMBER":

                # Path to the trajectory file.
                traj_file = "%s/%s.nc" % (self._process._work_dir, self._process._name)

                # MDTraj currently doesn't support the .prm7 extension, so we
                # need to copy the topology file to a temporary .parm7 file.
                # I've submitted a pull request to add support for the alternative
                # .prm7 extension.
                top_file = "%s/tmp.parm7" % self._process._work_dir
                _shutil.copyfile(self._top_file, top_file)

                # Set the topology format.
                top_format = "PARM7"

            # NAMD.
            elif self._process_name.upper() == "NAMD":

                # Path to the trajectory and topology files.
                traj_file = "%s/%s_out.dcd" % (self._process._work_dir, self._process._name)
                top_file = self._top_file

                # Set the topology format.
                top_format = "PDB"

            # Unsupported process.
            else:
                raise ValueError("BioSimSpace.Process.%s is unsupported!" % self._process_name)

            # Check that the trajectory and topology files exist.
            if not _os.path.isfile(traj_file):
                raise IOError("Trajectory file doesn't exist: '%s'" % traj_file)

            if not _os.path.isfile(top_file):
                raise IOError("Topology file doesn't exist: '%s'" % top_file)

            # Return an MDTraj object.
            if format == "mdtraj":

                # Create the MDTraj object.
                try:
                    traj = _mdtraj.load(traj_file, top=top_file)
                except:
                    _warnings.warn("MDTraj failed to read: traj=%s, top=%s" % (traj_file, top_file))
                    traj = None

                # Delete the temporary .parm7 file.
                if self._process_name.upper() is "AMBER":
                    _os.remove(top_file)

                return traj

            # Return an MDAnalysis Universe.
            else:
                try:
                    universe = _mdanalysis.Universe(top_file, traj_file, topology_format=top_format)
                except:
                    _warnings.warn("MDAnalysis failed to read: traj=%s, top=%s" % (traj_file, top_file))
                    universe = None

                # Delete the temporary .parm7 file.
                if self._process_name.upper() is "AMBER":
                    _os.remove(top_file)

                return universe

        # Read trajectory from file.
        # TODO:
        # Make this more robust, i.e. determine topology format from the file
        # extension and rename files if needed.
        else:
            # Return an MDTraj object.
            if format == "mdtraj":

                try:
                    traj = _mdtraj.load(self._traj_file, top=self._top_file)
                except:
                    _warnings.warn("MDTraj failed to read: traj=%s, top=%s" % (self._traj_file, self._top_file))
                    traj = None

                return traj

            # Return an MDAnalysis Universe.
            else:
                try:
                    universe = _mdanalysis.Universe(self._top_file, self._traj_file, topology_format=top_format)
                except:
                    _warnings.warn("MDAnalysis failed to read: traj=%s, top=%s" % (self._traj_file, self._top_file))
                    universe = None

                return universe

    def getFrames(self, indices=None):
        """Get trajectory frames as a list of System objects.

           Keyword arguments:

           indices -- A list of trajectory frame indices, or time stamps (in ns).
        """

        # The process is running. Grab the latest trajectory.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

            # There is no trajectory.
            if self._trajectory is None:
                return None

        # Work out the frame spacing in nanoseconds.
        # TODO:
        # How can we do this in a robust way if the trajectory is loaded from file?
        # Some formats do not store time information as part of the trajectory.
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

        # Store the number of frames.
        n_frames = self._trajectory.n_frames

        # Loop over all indices.
        for x in indices:

            # Make sure the frame index is within range.
            if x > 0 and x >= n_frames:
                raise ValueError("Frame index (%d) of of range (0 to %d)." %s (x, n_frames - 1))
            elif x < -n_frames:
                raise ValueError("Frame index (%d) of of range (-1 to -%d)." %s (x, n_frames))

            # The name of the frame coordinate file.
            frame_file = ".frame.nc"

            # Write the current frame as a NetCDF file.
            self._trajectory[x].save(frame_file)

            # Create a Sire system.
            system = _System(_Sire.IO.MoleculeParser.read([self._top_file, frame_file]))

            # Append the system to the list of frames.
            frames.append(system)

        # Remove the temporary frame coordinate file.
        _os.remove(frame_file)

        # Return the frames.
        return frames

    def nFrames(self):
        """Return the current number of trajectory frames."""

        # First get the current MDTraj object.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

        # There is no trajectory.
        if self._trajectory is None:
            return 0
        else:
            return self._trajectory.n_frames

    def RMSD(self, frame=None, atoms=None, molecule=None):
        """Compute the root mean squared displacement.

           Keyword arguments:

           frame    -- The index of the reference frame.
           atoms    -- A list of reference atom indices.
           molecule -- The index of the reference molecule.
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

        if molecule is not None and atoms is not None:
            _warnings.warn("Cannot have a reference molecule and list of atoms. Defaulting to atoms.")
            molecule = None

        # Extract the molecule from the reference trajectory frame and generate a
        # list of atom indices.
        if molecule is not None:
            # Integer molecule index.
            if type(molecule) is int:
                try:
                    molecule = self.getFrames(frame)[0]._getSireSystem()[_Sire.Mol.MolIdx(molecule)]
                except:
                    raise ValueError("Missing molecule index '%d' in System" % molecule)
            # Sire.Mol.MolIdx index.
            elif type(molecule) is _Sire.Mol.MolIdx:
                try:
                    molecule = self.getFrames(frame)[0]._getSireSystem()[molecule]
                except:
                    raise ValueError("Missing '%s' in System" % molecule.toString())

            # A BioSimSpace Molecule object.
            elif type(molecule) is _Molecule:
                molecule = molecule._getSireMolecule()
            # A Sire.Mol.Molecule object.
            elif type(molecule) is _Sire.Mol.Molecule:
                pass
            else:
                raise TypeError("'molecule' must be of type 'int', 'BioSimSpace.Molecue', 'Sire.Mol.MolIdx', or 'Sire.Mol.Molecule'")

            # Initialise the list of atom indices.
            atoms = []

            # Loop over all of the atoms and add them to the list.
            # Atoms are numbered from one in Sire, so subtract one to get the python index.
            for atom in molecule.atoms():
                atoms.append(atom.number().value() - 1)

        if atoms is not None:
            # Check that all of the atom indices are integers.
            if not all(isinstance(x, int) for x in atoms):
                raise TypeError("'atom' indices must be of type 'int'")

        # Use MDTraj to compute the RMSD.
        try:
            rmsd = _mdtraj.rmsd(self._trajectory, self._trajectory, frame, atoms)
        except:
            raise ValueError("Atom indices not found in the system.")

        # Convert to a list and return.
        return list(rmsd)

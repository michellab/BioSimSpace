"""
@package biosimspace
@author  Lester Hedges
@brief   A class for reading and manipulating biomolecular trajectories.
"""

from ..Process.process import Process

from Sire import try_import
from Sire.IO import MoleculeParser
from Sire.Mol import AtomNum, Molecule, MolIdx

from os import path, remove
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
            warn("Both a process and trajectory file are specified! Defaulting to 'process'.")
            self._traj_file = None

        # BioSimSpace process.
        if process is not None:
            if process.__class__.__base__ == Process:
                self._process = process
                self._process_name = process.__class__.__name__
                self._top_file = self._process._top_file

                # Check that the process can generate a trajectory.
                if not self._process._has_trajectory:
                    raise ValueError("BioSimSpace.Process.%s cannot generate a trajectory!" % self._process_name)

        # Trajectory and topology files.
        elif type(trajectory) is str and type(topology) is str:

            # Make sure the trajectory file exists.
            if not path.isfile(trajectory):
                raise IOError("Trajectory file doesn't exist: '%s'" % trajectory)

            # Make sure the topology file exists.
            if not path.isfile(topology):
                raise IOError("Topology file doesn't exist: '%s'" % topology)

            self._traj_file = trajectory
            self._top_file = topology

        # Invalid arguments.
        else:
            raise ValueError("BioSimSpace.Trajectory requires a BioSimSpace.Process object, "
                "or a trajectory and topology file.")

        # Get the current trajectory.
        self._trajectory = self.getTrajectory()

    def getTrajectory(self, format='mdtraj'):
        """Get the current trajectory object.

           Keyword arguments:

           format -- Whether to return a 'MDTraj' or 'MDAnalysis' object.
        """

        if format.upper() not in ["MDTRAJ", "MDANALYSIS"]:
            warn("Invalid trajectory format. Using default (mdtraj).")
            format = mdtraj

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
                copyfile(self._top_file, top_file)

                # Set the topology format.
                top_format = 'PARM7'

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
            if not path.isfile(traj_file):
                raise IOError("Trajectory file doesn't exist: '%s'" % traj_file)

            if not path.isfile(top_file):
                raise IOError("Topology file doesn't exist: '%s'" % top_file)

            # Return an MDTraj object.
            if format is 'mdtraj':

                # Create the MDTraj object.
                traj = mdtraj.load(traj_file, top=top_file)

                # Delete the temporary .parm7 file.
                if self._process_name.upper() is "AMBER":
                    remove(top_file)

                return traj

            # Return an MDAnalysis Universe.
            else:
                universe = mdanalysis.Universe(top_file, traj_file, topology_format=top_format)

                # Delete the temporary .parm7 file.
                if self._process_name.upper() is "AMBER":
                    remove(top_file)

                return universe

        # Read trajectory from file.
        # TODO:
        # Make this more robust, i.e. determine topology format from the file
        # extension and rename files if needed.
        else:
            # Return an MDTraj object.
            if format is "mdtraj":

                # Create the MDTraj object.
                return mdtraj.load(self._traj_file, top=self._top_file)

            # Return an MDAnalysis Universe.
            else:
                return mdanalysis.Universe(self._top_file, self._traj_file)

    def getFrames(self, indices=None):
        """Get trajectory frames as a list of Sire systems.

           Keyword arguments:

           indices -- A list of trajectory frame indices, or time stamps (in ns).
        """

        # The process is running. Grab the latest trajectory.
        if self._process is not None and self._process.isRunning():
            self._trajectory = self.getTrajectory()

            # There is no trajectory.
            if traj is None:
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
        max_frame = self._trajectory.n_frames - 1

        # Loop over all indices.
        for x in indices:

            # Make sure the frame index is within range.
            if abs(x) > max_frame:
                raise ValueError("Frame index (%d) of of range (0-%d)." %s (x, max_frame))

            # The name of the frame coordinate file.
            frame_file = ".frame.nc"

            # Write the current frame as a NetCDF file.
            self._trajectory[x].save(frame_file)

            # Create a Sire system.
            system = MoleculeParser.read([self._top_file, frame_file])

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
                # Store the maximum frame number.
                max_frame = self._trajectory.n_frames - 1
                if abs(frame) > max_frame:
                    raise ValueError("Frame index (%d) of of range (0-%d)." %s (frame, max_frame))

        if molecule is not None and atoms is not None:
            warn("Cannot have a reference molecule and list of atoms. Defaulting to atoms.")
            molecule = None

        # Extract the molecule from the reference trajectory frame and generate a
        # list of atom indices.
        if molecule is not None:
            # Integer molecule index.
            if type(molecule) is int:
                try:
                    molecule = self.getFrames(frame)[0][MolIdx(molecule)]
                except:
                    raise ValueError("Missing 'MolIdx(%d)' in Sire.System" % molecule)
            # Sire.Mol.MolIdx index.
            elif type(molecule) is MolIdx:
                try:
                    molecule = self.getFrames(frame)[0][molecule]
                except:
                    raise ValueError("Missing '%s' in Sire.System" % molecule.toString())
            # A Sire.Mol.Molecule object.
            elif type(molecule) is Molecule:
                pass
            else:
                raise TypeError("'molecule' must be of type 'int', 'Sire.Mol.MolIdx', or 'Sire.Mol.Molecule'")

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
            rmsd = mdtraj.rmsd(self._trajectory, self._trajectory, frame, atoms)
        except:
            raise ValueError("Atom indices not found in the system.")

        # Convert to a list and return.
        return list(rmsd)

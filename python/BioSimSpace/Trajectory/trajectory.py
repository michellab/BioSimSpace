"""
@package biosimspace
@author  Lester Hedges
@brief   A class for reading and manipulating biomolecular trajectories.
"""

from ..Process.process import Process

from Sire import try_import

from os import path
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

    def __init__(process=None, trajectory=None, topology=None):
        """Constructor.

           Keyword arguments:

           process    -- A BioSimSpace process object.
           trajectory -- A trajectory file.
           topology   -- A topology file.
        """

        # Nothing to create a trajectory from.
        if process is None and trajectory is None:
            raise ValueError("Both 'process' and 'trajectory' keyword arguments are 'None'")

        # Both use cases active. Default to process.
        if not process is None and not trajectory is None:
            warn("Both a process and trajectory file are specified! Defaulting to 'process'.")
            self._trajectory_file = None
            self._topology_file = None

        if not process is None:
            # BioSimSpace process.
            if process.__class__.__base__ == Process:
                self._process = process
            else:
                raise ValueError("Process must be of type 'BioSimSpace.Process!")
        else:
            self._process = None

            # Make sure the trajectory and topology files exist.

            if not path.isfile(trajectory):
                raise IOError(('Trajectory file doesn\'t exist: "{x}"').format(x=trajectory))
            else:
                self._trajectory_file = trajectory

            if not path.isfile(topology):
                raise IOError(('Topology file doesn\'t exist: "{x}"').format(x=topology))
            else:
                self._topology_file = topology

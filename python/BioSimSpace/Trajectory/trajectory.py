"""
@package biosimspace
@author  Lester Hedges
@brief   A class for reading and manipulating biomolecular trajectories.
"""

try:
    from Sire import try_import
    pygtail = try_import("mdtraj")
except ImportError:
    raise ImportError("MDTraj is not installed. Please install mdtraj in order to use BioSimSpace.")

class Trajectory():
    """A class for reading a manipulating biomolecular trajectories."""

    def __init__():
        return

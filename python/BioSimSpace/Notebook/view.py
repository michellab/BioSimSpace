"""
@package biosimspace
@author  Lester Hedges
@brief   A class for handling interactive molecular visualisations.
"""

from Sire.IO import PDB2
from Sire.System import System

from ..Process.process import Process

import tempfile

try:
    from Sire import try_import
    nv = try_import("nglview")
except ImportError:
    raise ImportError("NGLview is not installed. Please install nglview in order to use BioSimSpace.Notebook.")

class View():
    """A class for handling interactive molecular visualisations."""

    def __init__(self, handle):
        """Constructor.

           Keyword arguments:

           handle -- A handle to a Sire.System or BioSimSpace.Process object.
        """

        # Check the handle.

        # BioSimSpace process.
        if handle.__class__.__base__ == Process:
            self._handle = handle
            self._is_process = True

        # Sire system.
        elif handle.__class__ == System:
            self._handle = handle
            self._is_process = False

        else:
            raise ValueError("The handle must be a Sire.System or BioSimSpace.Process object.")

        # Create a temporary workspace for the view object.
        self._tmp_dir = tempfile.TemporaryDirectory()
        self._work_dir = self._tmp_dir.name

        # Zero the number of views.
        self._num_views = 0

    def system(self):
        """View the entire molecular system."""

        # Get the latest system from the process.
        if self._is_process:
            system = self._handle.getSystem()

            # No system.
            if system is None:
                return

        else:
            system = handle

        # Increment the number of views.
        self._num_views += 1

        # Create the file name.
        filename = "%s/view_%04d.pdb" % (self._work_dir, self._num_views)

        # Create a PDB object.
        pdb = PDB2(system)

        # Write the PDB to file.
        pdb.writeToFile(filename)

        # Create the NGLview.
        view = nv.show_file(filename, gui=True)
        view

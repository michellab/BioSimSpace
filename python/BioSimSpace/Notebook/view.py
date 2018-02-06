"""
@package biosimspace
@author  Lester Hedges
@brief   A class for handling interactive molecular visualisations.
"""

import Sire.Mol
import Sire.System

from Sire.IO import PDB2

from ..Process.process import Process

from os import chdir, getcwd
from shutil import copyfile

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
        elif handle.__class__ == Sire.System.System:
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

        # Create and return the view.
        return self._create_view(system)

    def molecule(self, index=0):
        """View a specific molecule.

           Keyword arguments:

           index -- The molecule index.
        """

        # Get the latest system from the process.
        if self._is_process:
            system = self._handle.getSystem()

            # No system.
            if system is None:
                return

        else:
            system = handle

        # Extract the molecule numbers.
        molnums = system.molNums()

        if index < 0 or index > len(molnums):
            raise ValueError("Molecule index is out of range!")

        # Create a new system and add a single molecule.
        s = Sire.System.System("BioSimSpace molecule")
        m = Sire.Mol.MoleculeGroup("all")
        m.add(system[molnums[index]])
        s._old_add(m)

        # Create and return the view.
        return self._create_view(s)

    def reload(self, index=None):
        """Reload a particular view.

           Keyword arguments:

           index -- The view index.
        """

        if index < 1 or index > self._num_views:
            raise ValueError("View index (%d) is out of range: [1--%d]" % (index, self._num_views))

        # Default to the most recent view.
        if index is None:
            index = self._num_views

        # Create and return the view.
        return self._create_view(view=index)

    def nViews(self):
        """Return the number of views."""
        return self._num_views

    def _create_view(self, system=None, view=None):
        """Helper function to create the NGLview object.

           Keyword arguments:

           system -- A Sire molecular system.
           view   -- The index of an existing view.
        """

        if system is None and view is None:
            raise ValueError("Both 'system' and 'view' cannot be 'None'.")

        elif system is not None and view is not None:
            raise ValueError("One of 'system' or 'view' must be 'None'.")

        # Increment the number of views.
        if view is None:
            self._num_views += 1
            view = self._num_views

        # Create the file name.
        filename = "%s/view_%04d.pdb" % (self._work_dir, view)

        # Create a PDB object and write to file.
        if not system is None:
            pdb = PDB2(system)
            pdb.writeToFile(filename)

        # Create the NGLview object.
        view = nv.show_file(filename)

        # Return the view and display it.
        return view.display()

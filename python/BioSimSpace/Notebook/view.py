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

import nglview
import tempfile

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

    def system(self, gui=True):
        """View the entire molecular system.

           Keyword arguments:

           gui -- Whether to display the gui.
        """

        # Get the latest system from the process.
        if self._is_process:
            system = self._handle.getSystem()

            # No system.
            if system is None:
                return

        else:
            system = handle

        # Create and return the view.
        return self._create_view(system, gui=gui)

    def molecules(self, indices=None, gui=True):
        """View specific molecules.

           Keyword arguments:

           indices -- A list of molecule indices.
           gui     -- Whether to display the gui.
        """

        # Return a view of the entire system.
        if indices is None:
            return self.system(gui=gui)

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

        # Create a new system.
        s = Sire.System.System("BioSimSpace molecule")
        m = Sire.Mol.MoleculeGroup("all")

        # Loop over all of the indices.
        for index in indices:
            if index < 0 or index > len(molnums):
                raise ValueError("Molecule index is out of range!")

            # Add the molecule.
            m.add(system[molnums[index]])

        # Add all of the molecules to the system.
        s._old_add(m)

        # Create and return the view.
        return self._create_view(s, gui=gui)

    def molecule(self, index=0, gui=True):
        """View a specific molecule.

           Keyword arguments:

           index -- The molecule index.
           gui   -- Whether to display the gui.
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

        # Make sure the index is valid.
        if index < 0 or index > len(molnums):
            raise ValueError("Molecule index is out of range!")

        # Create a new system and add a single molecule.
        s = Sire.System.System("BioSimSpace molecule")
        m = Sire.Mol.MoleculeGroup("all")
        m.add(system[molnums[index]])
        s._old_add(m)

        # Create and return the view.
        return self._create_view(s, gui=gui)

    def reload(self, index=None, gui=True):
        """Reload a particular view.

           Keyword arguments:

           index -- The view index.
           gui   -- Whether to display the gui.
        """

        if index < 1 or index > self._num_views:
            raise ValueError("View index (%d) is out of range: [1--%d]" % (index, self._num_views))

        # Default to the most recent view.
        if index is None:
            index = self._num_views

        # Create and return the view.
        return self._create_view(view=index, gui=gui)

    def nViews(self):
        """Return the number of views."""
        return self._num_views

    def _create_view(self, system=None, view=None, gui=True):
        """Helper function to create the NGLview object.

           Keyword arguments:

           system -- A Sire molecular system.
           view   -- The index of an existing view.
           gui    -- Whether to display the gui.
        """

        if system is None and view is None:
            raise ValueError("Both 'system' and 'view' cannot be 'None'.")

        elif system is not None and view is not None:
            raise ValueError("One of 'system' or 'view' must be 'None'.")

        # Make sure gui flag is valid.
        if gui not in [True, False]:
            gui = True

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
        view = nglview.show_file(filename)

        # Return the view and display it.
        return view.display(gui=gui)

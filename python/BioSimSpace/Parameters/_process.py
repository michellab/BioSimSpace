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
Functionality running parameterisation protocols as a background process.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from BioSimSpace import _is_notebook

from . import Protocol as _Protocol
from .._SireWrappers import Molecule as _Molecule

import glob as _glob
import os as _os
import queue as _queue
import sys as _sys
import tempfile as _tempfile
import threading as _threading
import zipfile as _zipfile

if _is_notebook():
    from IPython.display import FileLink as _FileLink

__all__ = ["Process"]

class Process():
    """A class for running parameterisation protocols as a background process."""

    def __init__(self, molecule, protocol, work_dir=None, autostart=False):
        """Constructor

           Positional arguments:

           molecule  -- The molecule to parameterise.
           protocol  -- The parameterisation protocol.

           Keyword arguments:

           work_dir  -- The working directory for the process.
           autostart -- Whether to automatically start the process.
        """

        # Validate arguments.

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

        if not isinstance(protocol, _Protocol._Protocol):
            raise TypeError("'protocol' must be of type 'BioSimSpace.Parameters.Protocol'")

        if work_dir is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        if type(autostart) is not bool:
            raise TypeError("'autostart' must be of type 'bool'")

        # Set attributes.
        self._molecule = molecule
        self._protocol = protocol
        self._new_molecule = None

        # Create a hash for the object.
        self._hash = hash((molecule, protocol)) % ((_sys.maxsize + 1) * 2)

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir)

        # Flag that the process hasn't been started.
        self._is_started = False

        # Initialise the queue and thread.
        self._queue = None
        self._thread = None

        # Start the process.
        if autostart:
            self.start()

    def start(self):
        """Start the process."""

        # Flag that the process has been started.
        if self._is_started:
            return None
        else:
            self._is_started = True

        # Create the queue.
        self._queue = _queue.Queue()

        # Run the protocol as a background thread.
        self._thread = _threading.Thread(target=self._protocol.run,
                                         args=[self._molecule, self._work_dir, self._queue])
        self._thread.start()

    def getMolecule(self):
        """Get the parameterised molecule. This method blocks until
           parameterisation is complete.
        """

        # Start the process, if it's not already started.
        if not self._is_started:
            self._start()

        # Block the thread until it finishes.
        if self._thread.is_alive():
            self._thread.join()

            # Get the parameterise molecule from the thread function.
            self._new_molecule = self._queue.get()

        # No molecule was return, parameterisation failed.
        if self._new_molecule is None:
            zipname = "%s.zip" % self._hash

            # Append the files to the archive.
            with _zipfile.ZipFile(zipname, "w") as zip:
                # Loop over all of the output files.
                for file in _glob.glob("%s/*" % self._work_dir):
                    zip.write(file, arcname=_os.path.basename(file))

            # Return a link to the archive.
            if _is_notebook():
                print("Parameterisation failed! Check output:")
                return _FileLink(zipname)
            # Return the name of the zip archive.
            else:
                print("Parameterisation failed! Check output: '%s.zip'" % self._hash)
                return None
        else:
            return self._new_molecule

    def getHash(self):
        """Get the object hash."""
        return self._hash

    def workDir(self):
        """Return the working directory."""
        return self._work_dir

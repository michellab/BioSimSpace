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

# TODO:
# Work out a way to safely kill running processes.
#
# This is hard because the process launches a thread which itself calls a
# Protocol.run method, in which multiple subprocesses can be launched. The
# thread would need access to the PID of the subprocess in order to kill them,
# and this would affect the logic of the run method (once a subprocess has been
# killed the method should exit).

# Alternatively, one could use a multiprocessing.Process instead of a thread,
# which has a terminate method. However, communication between the Process and
# the run method requires the return type of the method to be  picklable, which
# isn't the case for our Molecule object.

from BioSimSpace import _is_notebook

from . import Protocol as _Protocol
from .._Exceptions import ParameterisationError as _ParameterisationError
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

           Positional arguments
           --------------------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule to parameterise.

           protocol : BioSimSpace.Parameters.Protocol
               The parameterisation protocol.


           Keyword arguments
           -----------------

           work_dir : str
               The working directory for the process.

           autostart : bool
               Whether to automatically start the process.
        """

        # Validate arguments.

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

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
        self._is_error = False
        self._zipfile = None

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
                _os.makedirs(work_dir, exist_ok=True)

        # Flag that the process hasn't started/finished.
        self._is_started = False
        self._is_finished = False

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

        # Create the thread.
        self._thread = _threading.Thread(target=self._protocol.run,
                                         args=[self._molecule, self._work_dir, self._queue])

        # Start the thread.
        self._thread.start()

    def getMolecule(self):
        """Get the parameterised molecule. This method blocks until
           parameterisation is complete.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The parameterised molecule.
        """

        # Start the process, if it's not already started.
        if not self._is_started:
            self._start()

        # Block the thread until it finishes.
        if not self._is_finished:
            self._thread.join()

            # Get the parameterise molecule from the thread function.
            self._new_molecule = self._queue.get()

            # Flag that the thread has finished.
            self._is_finished = True

            # No molecule was return, parameterisation failed.
            if self._new_molecule is None:
                self._is_error = True
                raise _ParameterisationError("Parameterisation failed! "
                    + "Run 'getOutput() to see intermediate files.")
            else:
                # Fix the charges so that the total is integer values.
                self._new_molecule._fixCharge()

        return self._new_molecule

    def isError(self):
        """Return whether there was a parameterisation error.

           Returns
           -------

           is_error : bool
               Whether there was an error during parameterisation.
        """

        # Try to get the parameterised molecule.
        self.getMolecule()

        return self._is_error

    def getOutput(self, filename=None):
        """Return the output of the parameterisation process.

           Keyword arguments
           -----------------

           filename : str
               The name to write the output to.


           Returns
           -------

           file_link : str, IPython.lib.display.FileLink
               The name of, or link to, a zipfile containing the output.
        """

        # Try to get the molecule.
        self.getMolecule()

        if self._zipfile is None or filename is not None:
            if filename is not None:
                zipname = "%s.zip" % filename

                # If the user has passed a directory, make sure that is exists.
                if _os.path.basename(filename) != filename:
                    dirname = _os.path.dirname(filename)
                    # Create the directory if it doesn't already exist.
                    if not _os.path.isdir(dirname):
                        _os.makedirs(dirname, exist_ok=True)
            else:
                zipname = "%s.zip" % self._hash

            # Append the files to the archive.
            with _zipfile.ZipFile(zipname, "w") as zip:
                # Loop over all of the output files.
                for file in _glob.glob("%s/*" % self._work_dir):
                    zip.write(file, arcname=_os.path.basename(file))

            # Store the location of the zip file.
            self._zipfile = zipname

        # Return a link to the archive.
        if _is_notebook():
            return _FileLink(self._zipfile)
        else:
            return self._zipfile

    def getHash(self):
        """Get the object hash."""
        return self._hash

    def workDir(self):
        """Return the working directory."""
        return self._work_dir

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
Functionality running parameterisation protocols as a background process.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Process"]

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
# the run method requires the return type of the method to be picklable, which
# isn't the case for our Molecule object.

import glob as _glob
import os as _os
import queue as _queue
import sys as _sys
import tempfile as _tempfile
import threading as _threading
import warnings as _warnings
import zipfile as _zipfile

from BioSimSpace import _is_notebook
from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import ParameterisationError as _ParameterisationError
from BioSimSpace._SireWrappers import Molecule as _Molecule

from . import Protocol as _Protocol

if _is_notebook:
    from IPython.display import FileLink as _FileLink

def _wrap_protocol(protocol_function, process):
    """A simple decorator function to wrap the running of parameterisation
       protocols and catch exceptions.

       Parameters
       ----------

       protocol_function : function
           The protocol function.

       process : BioSimSpace.Parameters.Process
           A handle to the parent process.
    """
    try:
        protocol_function(process._molecule,
                          process._work_dir,
                          process._queue)
    except Exception as e:
        # Record that an error has been thrown.
        process._is_error = True
        process._last_error = e

        # Add None to the queue (no molecule).
        if process._queue is not None:
            process._queue.put(None)

        # Return to the user directory.
        _os.chdir(process._dir)

class Process():
    """A class for running parameterisation protocols as a background process."""

    def __init__(self, molecule, protocol, work_dir=None, auto_start=False):
        """Constructor

           Parameters
           ----------

           molecule : BioSimSpace._SireWrappers.Molecule, str
               The molecule to parameterise, either as a Molecule object, or
               SMILES string.

           protocol : BioSimSpace.Parameters.Protocol
               The parameterisation protocol.

           work_dir : str
               The working directory for the process.

           auto_start : bool
               Whether to automatically start the process.
        """

        # Validate arguments.

        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'")

        if not isinstance(protocol, _Protocol._Protocol):
            raise TypeError("'protocol' must be of type 'BioSimSpace.Parameters.Protocol'")

        if work_dir is not None and not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

        if not isinstance(auto_start, bool):
            raise TypeError("'auto_start' must be of type 'bool'")

        # Warn the user if the molecule contains no hydrogens.
        if isinstance(molecule, _Molecule):
            if molecule.search("element H").nResults() == 0:
                _warnings.warn("Attempting to parameterise a molecule without hydrogen atoms!")

        # Set attributes.
        self._molecule = molecule
        self._protocol = protocol
        self._new_molecule = None
        self._is_error = False
        self._last_error = None
        self._zipfile = None

        # Store the directory from which the process was launched.
        self._dir = _os.getcwd()

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
        if auto_start:
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
        self._thread = _threading.Thread(target=_wrap_protocol, args=[self._protocol.run, self])

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

            # Fix the charges so that the total is integer valued.
            if self._new_molecule is not None:
                self._new_molecule._fixCharge(property_map=self._protocol._property_map)

        # If there was an problem, return the last error.
        if self._is_error:
            if _isVerbose():
                raise _ParameterisationError("Parameterisation failed! Last error: '%s'" % str(self._last_error))
            else:
                raise _ParameterisationError("Parameterisation failed! Last error: '%s'" % str(self._last_error)) from None

        return self._new_molecule

    def isRunning(self):
        """Return whether the parameterising protocol is still running."""
        return not self._is_finished

    def isError(self):
        """Return whether there was a parameterisation error.

           Returns
           -------

           is_error : bool
               Whether there was an error during parameterisation.
        """
        return self._is_error

    def getOutput(self, filename=None, file_link=False):
        """Return a zip file containing the output from the parameterisation process.

           Parameters
           ----------

           filename : str
               The name to write the output to.

           file_link : bool
               Whether to return a FileLink when working in Jupyter.

           Returns
           -------

           file_link : str, IPython.lib.display.FileLink
               The name of, or link to, a zipfile containing the output.
        """

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

            # Store the location of the zip file. Only do so if the
            # parameterisation has finished.
            if not self._is_finished:
                self._zipfile = zipname

        # Return a link to the archive.
        if _is_notebook:
            if file_link and self._zipfile is not None:
                # Create a FileLink to the archive.
                f_link = _FileLink(self._zipfile)

                # Set the download attribute so that JupyterLab doesn't try to open the file.
                f_link.html_link_str = f"<a href='%s' target='_blank' download='{self._zipfile}'>%s</a>"

                # Return a link to the archive.
                return f_link
            else:
                return self._zipfile
        else:
            return self._zipfile

    def getHash(self):
        """Get the object hash."""
        return self._hash

    def workDir(self):
        """Return the working directory."""
        return self._work_dir

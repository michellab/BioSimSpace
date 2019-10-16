######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
Functionality for running background tasks.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Task"]

import glob as _glob
import os as _os
import tempfile as _tempfile
import threading as _threading
import zipfile as _zipfile

from BioSimSpace import _is_notebook

def _wrap_task(task):
    """A simple wrapper function to run a background tasks and catch exceptions.

       Parameters
       ----------

       task : :class:`Task <BioSimSpace.Process.Task>`
           A handle to the task object.
    """

    # Try to run the task and grab the result.
    try:
        task._result = task._run()

    # Catch any exception raised in the _run method.
    except Exception as e:
        task._result = e

class Task():
    """Base class for running a background task."""

    def __init__(self, name=None, work_dir=None, auto_start=False):
        """Constructor.

           Parameters
           ----------

           name : str
               The name of the task.

           work_dir : str
               The working directory for the task.

           auto_start : bool
               Whether to immediately start the task.
        """

        # Don't allow user to create an instance of this base class.
        if type(self) is Task:
            raise Exception("<Task> must be subclassed.")

        # Validate inputs.

        if name is not None and type(name) is not str:
            raise TypeError("'name' must be of type 'str'")

        if work_dir is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        if type(auto_start) is not bool:
            raise TypeError("'auto_start' must be of type 'bool'")

        # Initialise status flags.
        self._is_started = False
        self._is_finished = False
        self._is_error = False

        # Initialise the result of the task.
        self._result = None

        # Initialise the error message.
        self._error_message = None

        # Set the task name.
        self._name = name

        # Initialise the zip file name.
        self._zipfile = None

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._tmp_dir = None

            # Use full path.
            if work_dir[0] != "/":
                work_dir = _os.getcwd() + "/" + work_dir
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir, exist_ok=True)

        # Start the task.
        if auto_start:
            self.start()

    def start(self):
        """Start the task."""

        # The task is already running.
        if self._is_started:
            if not self._is_finished:
                return None

        # Update status flags.
        self._is_started = True
        self._is_finished = False
        self._is_error = False

        # Reset the error message.
        self._error_message = None

        # Create the thread.
        self._thread = _threading.Thread(target=_wrap_task, args=[self])

        # Start the task.
        self._thread.start()

    def workDir(self):
        """Return the working directory for the task.

           Returns
           -------

           work_dir : str
               The path of the working directory.
        """
        return self._work_dir

    def getResult(self):
        """Get the result of the task. This will block until the task finishes."""

        if not self._is_started:
            return None

        # Block until the task finishes.
        if not self._is_finished:
            self._thread.join()
            self._is_finished = True

        # If there was a problem, return the error message.
        if isinstance(self._result, Exception):
            self._is_error = True
            self._error_message = str(self._result)
            raise self._result

        return self._result

    def isStarted(self):
        """Return whether the task has been started.

           Returns
           -------

           is_started : bool
               Whether the task has been started.
        """
        return self._is_started

    def isError(self):
        """Return whether the task errored.

           Returns
           -------

           is_error : bool
               Whether the task errored.
        """
        return self._is_error

    def isFinished(self):
        """Return whether the task has finished.

           Returns
           -------

           is_finished : bool
               Whether the task has finished.
        """
        return self._is_finished

    def getException(self):
        """Return the exception.

           Returns
           -------

           exception : Exception
               The exception.
        """
        if self._is_error:
            return self._result
        else:
            return None

    def getErrorMessage(self):
        """Return the error message.

           Returns
           -------

           error : str
               The error message.
        """
        if self._is_error:
            return self._error_message
        else:
            return None

    def getOutput(self, filename=None, file_link=False):
        """Return a zip file containing the output of the task.

           Parameters
           ----------

           filename : str
               The name of the output archive.

           file_link : bool
               Whether to return a FileLink when working in Jupyter.

           Returns
           -------

           file_link : str, IPython.lib.display.FileLink
               The name of, or link to, a zipfile containing the output.
        """

        # The user has specified a working directory, return its name.
        if self._tmp_dir is None:
            return self._work_dir

        # Don't recreate an existing zip file.
        if self._zipfile is None:

            # Create the file name.
            if filename is None:
                if self._name is None:
                    filename = "output"
                else:
                    filename = self._name
            else:
                if type(filename) is not str:
                    raise TypeError("'filename' must be of type 'str'.")

            # Create the name of the zip file.
            zipname = "%s.zip" % filename

            # If the user has passed a directory, make sure that is exists.
            if _os.path.basename(filename) != filename:
                dirname = _os.path.dirname(filename)
                # Create the directory if it doesn't already exist.
                if not _os.path.isdir(dirname):
                    _os.makedirs(dirname, exist_ok=True)

            # Append the files to the archive.
            with _zipfile.ZipFile(zipname, "w") as zip:
                # Loop over all of the output files.
                for file in _glob.glob("%s/*" % self._work_dir):
                    zip.write(file, arcname=_os.path.basename(file))

            # Store the location of the zip file. Only do so if the
            # task has finished.
            if not self._is_finished:
                self._zipfile = zipname

        # Return a link to the archive.
        if _is_notebook:
            if file_link:
                return _FileLink(zipname)
            else:
                return zipname
        else:
            return zipname

    def _run(self):
        """User-defined method to run the specific background task."""

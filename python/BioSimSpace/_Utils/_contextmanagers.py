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
Custom context managers.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["cd", "stdout_redirected", "stderr_redirected"]

from contextlib import contextmanager as _contextmanager

import os as _os
import sys as _sys

# Adapted from: http://ralsina.me/weblog/posts/BB963.html
@_contextmanager
def cd(work_dir):
    """Execute the context in the directory "work_dir"

       Parameters
       ----------

       work_dir : str
           The working directory for the context.
    """

    # Validate the input.
    if type(work_dir) is not str:
        raise TypeError("'work_dir' must be of type 'str'")

    # Store the current directory.
    old_dir = _os.getcwd()

    # Create the working directory if it doesn't exist.
    if not _os.path.isdir(work_dir):
        _os.makedirs(work_dir)

    # Change to the new directory.
    _os.chdir(work_dir)

    # Execute the context.
    try:
        yield

    # Return to original directory.
    finally:
        _os.chdir(old_dir)

# Adapted from: https://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python/17954769#17954769
@_contextmanager
def stdout_redirected(to=_os.devnull):
    """ Redirect stdout to a file.

        Parameters
        ----------

        to : str
            The file to which stdout will be redirected.
    """

    # Validate the input.
    if type(to) is not str:
        raise TypeError("'to' must be of type 'str'")

    fd = sys.stdout.fileno()

    # Assert that Python and C stdio write using the same file descriptor
    # Assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        _sys.stdout.close()		    # + implicit flush()
        _os.dup2(to.fileno(), fd)	    # fd writes to 'to' file
        _sys.stdout = _os.fdopen(fd, "w")   # Python writes to fd

    with _os.fdopen(_os.dup(fd), "w") as old_stdout:
        with open(to, "w") as file:
            _redirect_stdout(to=file)
        try:
            yield   # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

# Adapted from: https://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python/17954769#17954769
@_contextmanager
def stderr_redirected(to=_os.devnull):
    """ Redirect stderr to a file.

        Parameters
        ----------

        to : str
            The file to which stderr will be redirected.
    """

    # Validate the input.
    if type(to) is not str:
        raise TypeError("'to' must be of type 'str'")

    fd = _sys.stderr.fileno()

    ##### Assert that Python and C stdio write using the same file descriptor
    ##### Assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stderr")) == fd == 1

    def _redirect_stderr(to):
        _sys.stderr.close()                 # + implicit flush()
        _os.dup2(to.fileno(), fd)	    # fd writes to 'to' file
        _sys.stderr = _os.fdopen(fd, "w")   # Python writes to fd

    with _os.fdopen(_os.dup(fd), "w") as old_stderr:
        with open(to, "w") as file:
            _redirect_stderr(to=file)
        try:
            yield   # allow code to be run with the redirected stderr
        finally:
            _redirect_stderr(to=old_stderr) # restore stderr.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

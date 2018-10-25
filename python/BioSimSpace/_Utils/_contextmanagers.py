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
Custom context managers.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from contextlib import contextmanager as _contextmanager
import os as _os

__all__ = ["cd"]

@_contextmanager
def cd(path):
    # Store the current directory.
    old_dir = _os.getcwd()

    # Create the new directory if it doesn't exist.
    if not _os.path.isdir(path):
        _os.makedirs(path)

    # Change to the new directory.
    _os.chdir(path)

    # Execute context.
    try:
        yield

    # Return to original directory.
    finally:
        _os.chdir(old_dir)

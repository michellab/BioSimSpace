######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
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

"""Custom context managers."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["create_workdir"]


import os as _os
import tempfile as _tempfile


def create_workdir(work_dir=None):
    """
    Execute the context in the directory "work_dir".

    Parameters
    ----------

    work_dir : str
        The working directory for the context. If None, then a temporary
        working directory will be created.


    Returns
    -------

    work_dir : str
        The absolute path to the working directory.

    tmp_dir : tempfile.TemporaryDirectory()
        The temporary directory object, if created. This will be None
        unless 'work_dir' passed to the function is None.
    """

    # Validate the input.
    if work_dir and not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    # Temporary directory.
    if work_dir is None:
        tmp_dir = _tempfile.TemporaryDirectory()
        work_dir = tmp_dir.name

    # User specified directory.
    else:
        tmp_dir = None

        # Use absolute path.
        if not _os.path.isabs(work_dir):
            work_dir = _os.path.abspath(work_dir)

        # Create the directory if it doesn't already exist.
        if not _os.path.isdir(work_dir):
            _os.makedirs(work_dir, exist_ok=True)

    return work_dir, tmp_dir

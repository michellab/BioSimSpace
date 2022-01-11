#####################################################################
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
Functionality for streaming wrapped Sire objects.
"""

__author__ = "Lester Hedges"
__email__  = "lester.hedges@gmail.com"

__all__ = ["save", "load"]

import os as _os
import uuid as _uuid

from Sire import Mol as _SireMol
from Sire import Qt as _Qt
from Sire import Stream as _SireStream
from Sire import System as _SireSystem

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import StreamError as _StreamError
from BioSimSpace._SireWrappers._sire_wrapper import SireWrapper as _SireWrapper
from BioSimSpace import _SireWrappers

def save(sire_object, filebase=None):
    """Save a wrapped Sire object to a binary data stream. Objects can be
       streamed to file, or to a bytes object, which will be returned.

       Parameters
       ----------

       sire_object : :class:`System <BioSimSpace._SireWrappers.SireWrapper>`
           A wrapped Sire object.

       filebase : str
           The base name of the binary output file. If none, then the object
           will be streamed to a bytes object, which will be returned.

       Returns
       -------

       stream : bytes
           The streamed object. None will be returned when streaming to file.
    """

    # Validate input.

    if not isinstance(sire_object, (_SireWrapper, _SireWrappers.SearchResult)):
        raise TypeError("'sire_object' must be of type 'BioSimSpace._SireWrappers.SireWrapper'.")

    if filebase is not None:
        if not isinstance(filebase, str):
            raise TypeError("'filebase' must be of type 'str'.")

    try:
        # Stream to a bytes object.
        if filebase is None:
            # Stream the object to an s3 file.
            filebase = "." + _uuid.uuid4().hex
            _SireStream.save(sire_object._sire_object, filebase + ".s3")

            # Read the file as a binary data stream.
            contents = open(filebase + ".s3", "rb").read()

            # Remove the intermediate file.
            _os.remove(filebase + ".s3")

            return contents

        # Stream to file.
        else:
            _SireStream.save(sire_object._sire_object, f"{filebase}.s3")

            return None

    except Exception as e:
        msg = f"Failed to stream {sire_object}."
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None

def load(stream):
    """Load a wrapped Sire object from a binary data stream. Objects can be
       streamed from file, or from a bytes object.

       Parameters
       ----------

       stream : str, bytes
           The path to the binary file / stream containing the object.
    """

    # Validate input.

    if isinstance(stream, str):
        is_file = True
    elif isinstance(stream, bytes):
        is_file = False
    else:
        raise TypeError("'stream' must be of type 'str' or 'bytes'.")

    if is_file:
        if not _os.path.isfile(stream):
            raise ValueError(f"The binary file '{file}' doesn't exist!")

    # Try to load the object.

    try:
        # Load a bytes object via an intermediate s3 file.
        if not is_file:
            # Write the binary data stream to an s3 file.
            filename = "." + _uuid.uuid4().hex + ".s3"
            with open(filename, "wb") as file:
                file.write(stream)

            # Load the binary stream and remove the intermediate file.
            sire_object = _SireStream.load(filename)
            _os.remove(filename)

        # Stream from file directly.
        else:
            sire_object = _SireStream.load(stream)

        # Construct the wrapped object.
        if isinstance(sire_object, _SireSystem.System):
            return _SireWrappers.System(sire_object)
        elif isinstance(sire_object, _SireMol.Molecule):
            return _SireWrappers.Molecule(sire_object)
        elif isinstance(sire_object, _SireMol.MoleculeGroup):
            return _SireWrappers.Molecules(sire_object)
        elif isinstance(sire_object, _SireMol.Residue):
            return _SireWrappers.Residue(sire_object)
        elif isinstance(sire_object, _SireMol.Atom):
            return _SireWrappers.Atom(sire_object)
        elif isinstance(sire_object, _SireMol.SelectResult):
            return _SireWrappers.SearchResult(sire_object)
        else:
            raise _StreamError(f"Unable to stream object of type {type(sire_object)}.")

    except Exception as e:
        msg = f"Failed to stream {object} from '{stream}'."
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None

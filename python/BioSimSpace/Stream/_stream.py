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

"""
Functionality for streaming wrapped Sire objects.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["save", "load"]

import os as _os

from sire.legacy import Mol as _SireMol
from sire.legacy import Stream as _SireStream
from sire.legacy import System as _SireSystem

from .. import _isVerbose
from .._Exceptions import StreamError as _StreamError
from .._SireWrappers._sire_wrapper import SireWrapper as _SireWrapper
from .. import _SireWrappers


def save(sire_object, filebase):
    """
    Stream a wrapped Sire object to file.

    Parameters
    ----------

    sire_object : :class:`System <BioSimSpace._SireWrappers.SireWrapper>`
        A wrapped Sire object.

    filebase : str
        The base name of the binary output file.
    """

    # Validate input.

    if not isinstance(sire_object, _SireWrapper) and not isinstance(
        sire_object, _SireWrappers.SearchResult
    ):
        raise TypeError(
            "'sire_object' must be of type 'BioSimSpace._SireWrappers.SireWrapper'."
        )

    if type(filebase) is not str:
        raise TypeError("'filebase' must be of type 'str'.")

    try:
        sire_object = _add_metadata(sire_object)
    except:
        raise _StreamError("Unable to add metadata to streamed object!")

    try:
        _SireStream.save(sire_object._sire_object, f"{filebase}.bss")
    except Exception as e:
        msg = f"Failed to stream {sire_object} to file '{filebase}.bss'."
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None


def load(file):
    """
    Stream a wrapped Sire object from file.

    Parameters
    ----------

    file : str
        The path to the binary file containing the streamed object.
    """

    # Validate input.

    if not _os.path.isfile(file):
        raise TypeError(f"The binary file '{file}' doesn't exist!")

    # Try to load the object.

    try:
        # Stream from file.
        sire_object = _SireStream.load(file)

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
        elif "Selector" in sire_object.__class__.__name__:
            return _SireWrappers.SearchResult(sire_object)
        else:
            raise _StreamError(f"Unable to stream object of type {type(sire_object)}.")

    except Exception as e:
        msg = f"Failed to stream {object} from file '{file}'."
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None


def _add_metadata(sire_object):
    """
    Internal function to tag a Sire object with metadata.

    Parameters
    ----------

    sire_object : :class:`System <BioSimSpace._SireWrappers.SireWrapper>`
        The wrapped Sire object to stream.

    Returns
    -------

    sire_object : :class:`System <BioSimSpace._SireWrappers.SireWrapper>`
        The tagged wrapped Sire object.
    """

    if not isinstance(sire_object, _SireWrapper) and not isinstance(
        sire_object, _SireWrappers.SearchResult
    ):
        raise TypeError(
            "'sire_object' must be of type 'BioSimSpace._SireWrappers.SireWrapper'."
        )

    from sire import __version__ as _sire_version
    from sire import __revisionid__ as _sire_revisionid
    from .. import __version__ as _bss_version

    # Create a copy of the object.
    try:
        _sire_object = sire_object.copy()
    except:
        _sire_object = sire_object

    # Work out the name of the Sandpit.
    try:
        sandpit = sire_object.__module__.split("Sandpit")[1].split(".")[1]
    except:
        sandpit = "None"

    # Generate the metadata.
    metadata = {
        "bss_version": _bss_version,
        "sire_version": _sire_version,
        "sire_revisionid": _sire_revisionid,
        "sandpit": sandpit,
    }

    # Tag the object.
    try:
        # Set a system level property.
        _sire_object._sire_object.setProperty("metadata", metadata)
    except:
        # Set an object (molecule, residue, or atom) level property.
        # Make sure to convert back to the correct type.
        try:
            _sire_object = _sire_object.__class__(
                sire_object._sire_object.edit()
                .setProperty("metadata", metadata)
                .molecule()
                .commit()
            )
        except:
            pass

    return _sire_object

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

"""Functionality for converting between molecular representations."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["supportedFormats", "to"]

from rdkit.Chem.rdchem import Mol as _RDMol

from sire import convert as _sire_convert

import sire.legacy.Mol as _SireMol
import sire.legacy.System as _SireSystem

from .._Exceptions import ConversionError as _ConversionError
from .. import _SireWrappers


def supportedFormats():
    """
    Return a list of the supported formats.

    Returns
    -------

    formats : [str]
        The supported formats.
    """
    return _sire_convert.supported_formats()


def to(obj, format="biosimspace", property_map={}):
    """
    Convert an object to a specified format.

    Parameters
    ----------

    obj :
        The input object to convert.

    format : str
        The format to convert to.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    converted_obj :
       The object in the converted format.
    """

    # Validate the input.

    if not isinstance(format, str):
        raise TypeError("'format' must be of type 'str'.")

    # Convert to lower case and strip whitespace.
    format = format.lower().replace(" ", "")

    if not format in supportedFormats():
        raise ValueError(
            f"Unsupported format '{format}', options are: {', '.join(supportedFormats())}."
        )

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    # Now do some work to handle different types of input. Where possible, we
    # try to convert to the "expected" format, i.e. atom to atom, residue
    # to residue, etc. This means that we need to use sire.legacy.Mol.PartialMolecule
    # to extract the sub-components we require.

    # BioSimSpace objects, i.e. wrapped Sire objects.
    if isinstance(obj, _SireWrappers._sire_wrapper.SireWrapper):
        # Return the object itself.
        if format == "biosimspace":
            return obj

        # Return the underlying object.
        elif format == "sire":
            if isinstance(obj, _SireWrappers.Molecules):
                return _SireMol.SelectorMol(obj._sire_object)
            else:
                return obj._sire_object

        elif format == "rdkit":
            # Handle System and Molecules objects.
            if isinstance(obj, (_SireWrappers.System, _SireWrappers.Molecules)):
                try:
                    return _sire_convert.to(
                        obj._sire_object.molecules(), "rdkit", map=property_map
                    )
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")
            # Handle Residues and Atoms by extracting the PartialMolecule.
            else:
                try:
                    return _sire_convert.to(
                        _SireMol.PartialMolecule(obj._sire_object).extract(),
                        "rdkit",
                        map=property_map,
                    )
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

    # Sire objects.
    elif isinstance(
        obj,
        (
            _SireSystem.System,
            _SireMol.SelectorMol,
            _SireMol.MoleculeGroup,
            _SireMol.Molecules,
            _SireMol.Molecule,
            _SireMol.Residue,
            _SireMol.Atom,
        ),
    ):
        if format == "biosimspace":
            if isinstance(obj, _SireSystem.System):
                try:
                    return _SireWrappers.System(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.SelectorMol):
                try:
                    return _SireWrappers.Molecules(obj.toMoleculeGroup())
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Molecules):
                try:
                    molgrp = _SireMol.MoleculeGroup("all")
                    molgrp.add(obj)
                    return _SireWrappers.Molecules(molgrp)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.MoleculeGroup):
                try:
                    return _SireWrappers.Molecules(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Molecule):
                try:
                    return _SireWrappers.Molecule(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Residue):
                try:
                    return _SireWrappers.Residue(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Atom):
                try:
                    return _SireWrappers.Atom(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

        elif format == "rdkit":
            if isinstance(obj, (_SireSystem.System, _SireMol.MoleculeGroup)):
                try:
                    return _sire_convert.to(obj.molecules(), "rdkit", map=property_map)
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

            elif isinstance(obj, (_SireMol.Molecules, _SireMol.Molecule)):
                try:
                    return _sire_convert.to(obj, "rdkit", map=property_map)
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

            elif isinstance(obj, (_SireMol.Residue, _SireMol.Atom)):
                try:
                    return _sire_convert.to(
                        _SireMol.PartialMolecule(obj).extract(),
                        "rdkit",
                        map=property_map,
                    )
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

    # RDKit objects.
    elif isinstance(obj, _RDMol):
        try:
            # Convert the object, then return as a single residue, atom,
            # or entire molecule.
            new_obj = _sire_convert.to(obj, format, map=property_map)

            if format == "biosimspace":
                if new_obj.nAtoms() == 1:
                    return new_obj.getAtoms()[0]
                elif new_obj.nResidues() == 1:
                    return new_obj.getResidues()[0]
                else:
                    return new_obj

            elif format == "sire":
                if new_obj.nAtoms() == 1:
                    return new_obj.atoms()[0]
                if new_obj.nResidues() == 1:
                    return new_obj.residues()[0]
                else:
                    return new_obj
        except:
            raise _ConversionError(f"Unable to convert object to {format} format!")

    # Lists or tuples of objects.
    elif isinstance(obj, (list, tuple)):
        initial_type = type(obj[0])
        if not all(isinstance(x, initial_type) for x in obj):
            raise ValueError("Containers must contain objects of the same type!")

        if isinstance(obj[0], (_RDMol, _SireWrappers.Molecule)):
            try:
                new_obj = _sire_convert.to(obj, format, map=property_map)
                if format == "biosimspace":
                    new_obj = _SireWrappers.Molecules(new_obj)
                return new_obj
            except:
                raise _ConversionError(f"Unable to convert object to {format} format!")

        else:
            raise TypeError(f"Cannot convert object of type '{type(obj)}'")

    else:
        raise TypeError(f"Cannot convert object of type '{type(obj)}'")

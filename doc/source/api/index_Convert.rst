.. _ref-Box:

BioSimSpace.Convert
===================

The *Convert* package contains tools for converting between different molecular
representations. Due to the different way that molecules are represented, full
round trip conversion isn't yet possible. For example, starting with a BioSimSpace
system object, you perform multiple conversions and end up back with a list of
BioSimSpace molecules, i.e. those in the original system. Where possible, we try
to return converted objects as the smallest possible object, i.e. a single residue
molecule will be returned as a residue, a single atom residue will be returned as
an atom. As always, you can use the native functionality of the returned object
to convert to a different object if desired.

We currently only support conversion to an `OpenMM <https://openmm.org>`_ context
from specific BioSimSpace and `Sire <https://github.com/openbiosim/sire>`_ objects.
At present the conversion is one-way only.

.. automodule:: BioSimSpace.Convert

.. toctree::
   :maxdepth: 1

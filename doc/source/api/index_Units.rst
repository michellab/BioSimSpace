.. _ref-Units:

BioSimSpace.Units
=================

The *Units* package provides a set of common physical units required by BioSimSpace.
This provides convenient, short-hand access to :ref:`ref-Types`.

Some examples:

.. code-block:: python

   import BioSimSpace as BSS

   # Create a length of 14.3 Angstrom.
   length = BSS.Types.Length(14.3, "Angstrom") # Long-winded way.
   length = 14.3*BSS.Units.Length.angstrom     # Simplified way.

   # Create an area by multiplying two lengths.
   area = 3.6*BSS.Units.Length.nanometer * 12*BSS.Units.Length.angstrom

   # Create a time of 100 milliseconds.
   time = 100*BSS.Units.Time.millisecond

.. automodule:: BioSimSpace.Units

.. toctree::
   :maxdepth: 1

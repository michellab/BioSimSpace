
BioSimSpace.Types
=================

This package provides additional data types required by BioSimSpace. These are
types with units, such as temperature and time. This allows for type validation,
e.g. when passing a ``Temperature`` to a ``Protocol`` object, as well as conversion
between units.

The package is intrinsically linked with the `BioSimSpace.Units <../Units>`_
package, which provides a convenient shorthand for generating objects of a
given physical unit.

As an example:

.. code-block:: python

   import BioSimSpace as BSS.

   # Directly create an object of type length.
   length = BSS.Types.Length(10, "ANGSTROM")

   # Alternatively, we could have done.
   length = BSS.Types.Length("10 A")

   # or
   length = 10 * BSS.Units.Length.angstrom

   # Convert the length to centimeters.
   cms = length.centimeters()

   # Multiply two lengths to create an area.
   area = length * BSS.Units.Length.picometer

   # Multiply by an area to create a volume.
   volume = length * area

   # Get the magnitude and unit of the volume.
   magnitude = volume.magnitude()
   unit = volume.unit()

   # Convert the volume to meters cubed.
   m3 = volume.meters3()

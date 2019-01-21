
BioSimSpace.Units
=================

This package provides a set of common physical units required by BioSimSpace.

Units are intrinsically linked to `BioSimSpace.Types <../Types>`_\ , i.e. units can
be used to generate an object of the type associated with the physical unit.

For example:

.. code-block:: python

   import BioSimSpace as BSS

   # This will return a temperature object.
   temp = 300 * BSS.Units.Temperature.kelvin

   # This will return a length object.
   length = 10 * BSS.Units.Length.angstrom

``Length``
--------------


* ``meter``
* ``centimeter``
* ``millimeter``
* ``nanometer``
* ``angstrom``
* ``picometer``

``Area``
------------


* ``meter2``
* ``nanometer2``
* ``angstrom2``
* ``picometer2``

``Volume``
--------------


* ``meter3``
* ``nanometer3``
* ``angstrom3``
* ``picometer3``

``Energy``
--------------


* ``kcal_per_mol``
* ``kj_per_mol``
* ``kt``

``Pressure``
----------------


* ``atm``
* ``bar``

``Time``
------------


* ``day``
* ``hour``
* ``minute``
* ``second``
* ``millisecond``
* ``nanosecond``
* ``picosecond``
* ``femtosecond``

``Temperature``
-------------------


* ``celsius``
* ``fahrenheit``
* ``kelvin``

Note that temperature units are non-multiplicative. This is because the relations
between the different units include both a scale factor and an offset. Consequently
there is ambiguity when performing operations involving temperature types.

For example, ``10C + 20C`` could be interpreted as ``30C`` or ``576.3K``\ , i.e.
``283.15K + 293.15K``. There are similar problems for multiplication and division,
e.g. is ``10 * 1C`` equal to ``10C``\ , or ``283.15K``\ ? Because of this, BioSimSpace will
report an error when certain operations are performed:

.. code-block:: python

   ValueError: Ambiguous operation with offset unit: 'CELSIUS'

(Note that this is inspired by the excellent units package `Pint <http://pint.readthedocs.io/en/latest/nonmult.html>`_.)

To perform operations, one can explicitly use the Type interface:

.. code-block:: python

   # Add two temperatures together in Celsius to create a new temperature object.
   temp = BSS.Types.Temperature(BSS.Types.Temperature(10, "C").magnitude() +
                                BSS.Types.Temperature(80, "F").celsius().magnitude(), "C")

Alternatively, one can set the ``allow_offset`` variable to ``True``
to allow operations with offset units, e.g.

.. code-block:: python

   # Allow operations with offset units.
   BSS.Units.allow_offset = True

   # Can now perform regular multiplication and addition with offset units (result = 30C).
   temp = 10*BSS.Units.Temperature.celsius + 20*BSS.Units.Temperature.celsius

   # When mixing units, the left-hand operand takes precendence, i.e. the right-hand
   # is converted to the same type as the left before the operation.

   # Result = 3.3333C
   temp = 10*BSS.Units.Temperature.celsius + 20*BSS.Units.Temperature.fahrenheit

   # Result = 70F
   temp = 20*BSS.Units.Temperature.fahrenheit + 10*BSS.Units.Temperature.celsius

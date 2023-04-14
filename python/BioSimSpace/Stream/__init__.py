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
.. currentmodule:: BioSimSpace.Stream

Functions
=========

.. autosummary::
    :toctree: generated/

    save
    load

Examples
========

Stream a :class:`System <BioSimSpace._SireWrappers.System>` object to and from
file.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecular system.
   system0 = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"])

   # Stream to file.
   BSS.Stream.save(system0, "system")

   # Alternatively, stream the system object directly.
   system0.save("system")

   # Stream from file.
   system1 = BSS.Stream.load("system.s3")
"""

from ._stream import *

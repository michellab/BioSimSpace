######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
.. currentmodule:: BioSimSpace.Units

.. autosummary::
    :toctree: generated/

    Area.meter2
    Area.angstrom2
    Area.nanometer2
    Area.picometer2

    Charge.electron_charge
    Charge.coulomb

    Energy.kcal_per_mol
    Energy.kj_per_mol
    Energy.kt

    Length.meter
    Length.angstrom
    Length.nanometer
    Length.picometer

    Pressure.atm
    Pressure.bar

    Temperature.kelvin
    Temperature.celsius
    Temperature.fahrenheit

    Time.day
    Time.hour
    Time.minute
    Time.second
    Time.millisecond
    Time.nanosecond
    Time.picosecond
    Time.femtosecond

    Volume.meter3
    Volume.angstrom3
    Volume.nanometer3
    Volume.picometer3
"""

from . import Area
from . import Charge
from . import Energy
from . import Length
from . import Pressure
from . import Temperature
from . import Time
from . import Volume

# Whether to allow operations between offset units, see here for details:
# http://pint.readthedocs.io/en/latest/nonmult.html
allow_offset = False

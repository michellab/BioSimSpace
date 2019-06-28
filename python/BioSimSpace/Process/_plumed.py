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
Utility functions for interfacing with PLUMED.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["createPlumedConfig"]

from ..Metadynamics import CollectiveVariable as _CollectiveVariable
from ..Protocol import Metadynamics as _Metadynamics
from .._SireWrappers import System as _System

def createPlumedConfig(system, protocol):
    """Create PLUMED configuration file parameters.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           A BioSimSpace system object.

       protocol : :class:`Protocol.Metadynamics <BioSimSpace.Protocol.Metadynamics>`
           The metadynamics protocol.

       Returns
       -------

       config : [str]
           The list of PLUMED configuration strings.
    """

    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    if type(protocol) is not _Metadynamics:
        raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.Metadynamics'")

    config = []

    return config

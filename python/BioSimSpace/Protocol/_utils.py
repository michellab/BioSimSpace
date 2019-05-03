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
Utility functions.
"""

from ._custom import *
from ._equilibration import *
from ._free_energy import *
from ._minimisation import *
from ._production import *

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["protocols", "createProtocol"]

# Create a list of the supported protocols.
_protocols = []
_protocol_dict = {}
import sys as _sys
_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var != "Custom":
        _protocols.append(_var)
        _protocol_dict[_var] = getattr(_namespace, _var)
del _namespace
del _sys
del _var

def protocols():
    """Return a list of the supported Molecular Dynamics protocols.

       Returns
       -------

       protocols : [str]
          A list of the supported Molecular Dynamics protocols.
    """
    return _protocols

def createProtocol(protocol):
    """Create a default simulation protocol.

       Parameters
       ----------

       protocol : str
           The name of the simulation protocol.

       Returns
       --------

       protocol : :class:`Protocol <BioSimSpace.Protocol>`
           The chosen simulation protocol.
    """

    if protocol not in _protocols:
        raise KeyError("Unsupported protocol '%s', supported protocols are %s" % (protocol, _protocols))

    return _protocol_dict[protocol]()

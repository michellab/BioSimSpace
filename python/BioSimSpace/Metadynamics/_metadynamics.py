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
Functionality for initialising metadynamics simulation processes.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["run"]

from BioSimSpace._SireWrappers import System as _System
from BioSimSpace import Process as _Process
from BioSimSpace import Protocol as _Protocol

def run(system, protocol, auto_start=True, name="metamd", work_dir=None,
     seed=None, property_map={}):
    """Auto-configure and run a metadynamics process.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The molecular system.

       protocol : :class:`Protocol <BioSimSpace.Protocol.Metadynamics>`
           The metadynamics protocol.

       auto_start : bool
           Whether to start the process automatically.

       name : str
           The name of the process.

       work_dir : str
           The working directory for the process.

       seed : int
           A random number seed.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       process : :class:`Process <BioSimSpace.Process>`
           A process to run the molecular dynamics protocol.
    """

    # Check that the system is valid.
    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    # Check that the protocol is valid.
    if type(protocol) is not _Protocol.Metadynamics:
        raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.Metadynamics'")

    # Validate optional arguments.

    if type(auto_start) is not bool:
        raise TypeError("'auto_start' must be of type 'bool'")

    if type(name) is not str:
        raise TypeError("'name' must be of type 'str'")

    if work_dir is not None:
        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

    if seed is not None:
        if type(seed) is not int:
            raise TypeError("'seed' must be of type 'int'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create the process object.

    process = _Process.Gromacs(system, protocol, name=name,
        work_dir=work_dir, seed=seed, property_map=property_map)

    # Start the process.
    if auto_start:
        return process.start()
    else:
        return process

######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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
__email__ = "lester.hedges@gmail.com"

__all__ = ["run"]

from .._SireWrappers import System as _System
from .. import Process as _Process
from .. import Protocol as _Protocol

# Import common objects from BioSimSpace.MD._md
from ..MD._md import _file_extensions, _md_engines, _find_md_engines

def run(system, protocol, engine="AUTO", gpu_support=False, auto_start=True,
    name="metamd", work_dir=None, seed=None, property_map={},
    ignore_warnings=False, show_errors=True):
    """Auto-configure and run a metadynamics process.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The molecular system.

       protocol : :class:`Protocol <BioSimSpace.Protocol.Metadynamics>`
           The metadynamics protocol.

       engine : str
           The molecular dynamics engine to use. If "AUTO", then a matching
           engine will automatically be chosen. Supported engines can be
           found using 'BioSimSpace.Metadynamics.engines()'.

       gpu_support : bool
           Whether to choose an engine with GPU support.

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

       ignore_warnings : bool
           Whether to ignore warnings when generating the binary run file.
           This option is specific to GROMACS and will be ignored when a
           different molecular dynamics engine is chosen.

       show_errors : bool
           Whether to show warning/error messages when generating the binary
           run file. This option is specific to GROMACS and will be ignored
           when a different molecular dynamics engine is chosen.

       Returns
       -------

       process : :class:`Process <BioSimSpace.Process>`
           A process to run the molecular dynamics protocol.
    """

    # Check that the system is valid.
    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    # Check that the protocol is valid.
    if not isinstance(protocol, _Protocol.Metadynamics):
        raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.Metadynamics'")

    # Validate optional arguments.

    if not isinstance(engine, str):
        raise TypeError("'engine' must be of type 'str'.")
    md_engine = engine.upper().replace(" ", "")

    if not isinstance(gpu_support, bool):
        raise TypeError("'gpu_support' must be of type 'bool'")

    if not isinstance(auto_start, bool):
        raise TypeError("'auto_start' must be of type 'bool'")

    if not isinstance(name, str):
        raise TypeError("'name' must be of type 'str'")

    if work_dir is not None:
        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

    if seed is not None:
        if not type(seed) is int:
            raise TypeError("'seed' must be of type 'int'")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    if not isinstance(ignore_warnings, bool):
        raise ValueError("'ignore_warnings' must be of type 'bool.")

    if not isinstance(show_errors, bool):
        raise ValueError("'show_errors' must be of type 'bool.")

    # Find a molecular dynamics engine and executable.
    engines, exes = _find_md_engines(system, protocol, md_engine, gpu_support)

    # Create the process object, return the first supported engine that can
    # instantiate a process.

    for engine, exe in zip(engines, exes):
        try:
            # AMBER.
            if engine == "AMBER":
                process = _Process.Amber(system, protocol, exe=exe, name=name,
                    work_dir=work_dir, seed=seed, property_map=property_map)

            # GROMACS.
            elif engine == "GROMACS":
                process = _Process.Gromacs(system, protocol, exe=exe, name=name,
                    work_dir=work_dir, seed=seed, property_map=property_map,
                    ignore_warnings=ignore_warnings, show_errors=show_errors)

            # OPENMM.
            elif engine == "OPENMM":
                if gpu_support:
                    platform = "CUDA"
                else:
                    platform = "CPU"
                # Don't pass the executable name through so that this works on Windows too.
                process = _Process.OpenMM(system, protocol, exe=None, name=name,
                    work_dir=work_dir, seed=seed, property_map=property_map, platform=platform)

            # Start the process.
            if auto_start:
                return process.start()
            else:
                return process

        except:
            pass

    # If we got here, then we couldn't create a process.
    if md_engine == "AUTO":
        raise Exception(f"Unable to create a process using any supported engine: {engines}")
    else:
        raise Exception(f"Unable to create a process using the chosen engine: {md_engine}")

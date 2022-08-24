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
Functionality for absolute free-energy simulations. Extremely similar relative,
but adds functionality for dealing with receptor-ligand restraints.
"""

__all__ = ["Absolute", "getData"]

from collections import defaultdict as _defaultdict, OrderedDict as _OrderedDict
import copy as _copy
from glob import glob as _glob
import math as _math
import shlex as _shlex
import sys as _sys
import os as _os
import re as _re
import shutil as _shutil
import subprocess as _subprocess
import tempfile as _tempfile
import warnings as _warnings
import zipfile as _zipfile

try:
    import alchemlyb as _alchemlyb
    from alchemlyb.postprocessors.units import R_kJmol, kJ2kcal
    from alchemlyb.parsing.gmx import extract_u_nk as _gmx_extract_u_nk
    from alchemlyb.parsing.gmx import extract_dHdl as _gmx_extract_dHdl
    from alchemlyb.parsing.amber import extract_u_nk as _amber_extract_u_nk
    from alchemlyb.parsing.amber import extract_dHdl as _amber_extract_dHdl
    from alchemlyb.preprocessing.subsampling import statistical_inefficiency as _statistical_inefficiency
    from alchemlyb.estimators import AutoMBAR as _AutoMBAR
    from alchemlyb.estimators import TI as _TI
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol
    is_alchemlyb = True
except:
    print('Please install alchemlyb via pip for analysis using it.')
    is_alchemlyb = False

import numpy as _np
import pandas as _pd 
from Sire.Base import getBinDir as _getBinDir
from Sire.Base import getShareDir as _getShareDir

from .. import _gmx_exe
from .. import _is_notebook
from .._Exceptions import AnalysisError as _AnalysisError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from .._Utils import cd as _cd
from .. import Process as _Process
from .. import Protocol as _Protocol
from .. import Types as _Types
from .. import Units as _Units
from ._restraint import Restraint as _Restraint
from ._relative import Relative as _Relative

from ..MD._md import _find_md_engines

if _is_notebook:
    from IPython.display import FileLink as _FileLink

# Check that the analyse_freenrg script exists.
if _sys.platform != "win32":
    _analyse_freenrg = _os.path.join(_getBinDir(), "analyse_freenrg")
else:
    _analyse_freenrg = _os.path.join(_os.path.normpath(_getShareDir()), "scripts", "analyse_freenrg.py")
if not _os.path.isfile(_analyse_freenrg):
    raise _MissingSoftwareError("Cannot find free energy analysis script in expected location: '%s'" % _analyse_freenrg)
if _sys.platform == "win32":
    _analyse_freenrg = "%s %s" % (_os.path.join(_os.path.normpath(_getBinDir()), "sire_python.exe"), _analyse_freenrg)

class Absolute(_Relative):
    """Class for configuring and running absolute free-energy perturbation simulations."""

    # Create a list of supported molecular dynamics engines.
    _engines = ["GROMACS", "SOMD"]

    def __init__(self, system, protocol=None, work_dir=None, engine=None,
            gpu_support=False, setup_only=False, ignore_warnings=False,
            show_errors=True, extra_options=None, extra_lines=None,
            estimator='MBAR', restraint=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the perturbation leg. This must contain
               a single perturbable molecule and is assumed to have already
               been equilibrated.

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`, \
               The simulation protocol.

           work_dir : str
               The working directory for the free-energy perturbation
               simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "AMBER", "GROMACS", or "SOMD". If this argument is omitted
               then BioSimSpace will choose an appropriate engine for you.

           gpu_support : bool
               Whether the engine must have GPU support.

           setup_only: bool
               Whether to only support simulation setup. If True, then no
               simulation processes objects will be created, only the directory
               hierarchy and input files to run a simulation externally. This
               can be useful when you don't intend to use BioSimSpace to run
               the simulation. Note that a 'work_dir' must also be specified.

           ignore_warnings : bool
               Whether to ignore warnings when generating the binary run file.
               This option is specific to GROMACS and will be ignored when a
               different molecular dynamics engine is chosen.

           show_errors : bool
               Whether to show warning/error messages when generating the binary
               run file. This option is specific to GROMACS and will be ignored
               when a different molecular dynamics engine is chosen.

           extra_options : dict
               A dictionary containing extra options. Overrides the ones generated from the protocol.

           extra_lines : list
               A list of extra lines to be put at the end of the script.
           
           estimator : str
               Estimator used for the analysis - must be either 'MBAR' or 'TI'.

           restraint : :class:`Restraint <BioSimSpace.FreeEnergy.Restraint>`
               The Restraint object that contains information for the ABFE
               calculations.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """
        # Call the base class constructor.
        super().__init__(system, protocol, work_dir, engine,
                         gpu_support, setup_only, ignore_warnings,
                         show_errors, extra_options, extra_lines,
                         estimator, property_map)

        # Strip whitespace from engine and convert to upper case.
        engine = engine.replace(" ", "").upper()
        if engine not in ['GROMACS', 'SOMD']:
            raise NotImplementedError(f'ABFE functionality for MD Engine {engine} not implemented.')

        # Check that if a restraint is passed (bound leg simulation) it is valid.
        # For free leg simulations, the restraint will be None.
        if restraint is not None:
            if not isinstance(restraint, _Restraint):
                raise TypeError("'restraint' must be of type 'BioSimSpace.FreeEnergy.Restraint'.")
            else:
                # Ensure that the system is compatible with the restraint
                restraint.system = self._system

        self._restraint = restraint

        # Re-initialise the process runner to overwrite self._restraint = None (set in Relative)
        self._initialise_runner(self._system)

def getData(name="data", file_link=False, work_dir=None):
    """Return a link to a zip file containing the data files required for
       post-simulation analysis.

       Parameters
       ----------

       name : str
           The name of the zip file.

       file_link : bool
           Whether to return a FileLink when working in Jupyter.

       work_dir : str
           The working directory for the simulation.

       Returns
       -------

       output : str, IPython.display.FileLink
           A path, or file link, to an archive of the process input.
    """
    return Absolute.getData(name=name, file_link=file_link, work_dir=work_dir)

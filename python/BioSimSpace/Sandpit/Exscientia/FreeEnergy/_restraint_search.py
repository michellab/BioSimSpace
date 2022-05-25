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
Functionality for selecting receptor-ligand restraints from unrestrained
simulations.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["RestraintSearch"]

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
import MDAnalysis as mda
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

class RestraintSearch():
    """Class for running unrestrained simulations from which receptor-ligand 
    restraints are selected"""

    # Create a list of supported molecular dynamics engines.
    _engines = ["AMBER", "GROMACS", "SOMD"]

    def __init__(self, system, protocol=None, work_dir=None, engine=None,
            gpu_support=False, setup_only=False, ignore_warnings=False,
            show_errors=True, extra_options=None, extra_lines=None,
            property_map={}):
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
           
           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Validate the input.

        if not isinstance(system, _System):
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        if protocol is not None:
            if isinstance(protocol, _Protocol.Production):
                self._protocol = protocol
            else:
                raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.Production'")
        else:
            # Use a default protocol.
            self._protocol = _Protocol.Equilibration()

        self._extra_options = extra_options if extra_options is not None else {}
        self._extra_lines = extra_lines if extra_lines is not None else []

        if not isinstance(setup_only, bool):
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            if setup_only:
                raise ValueError("A 'work_dir' must be specified when 'setup_only' is True!")
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir, exist_ok=True)

        # There must be a single molecule to be decoupled (or annihilated).
        if system.nDecoupledMolecules() != 1:
            raise ValueError("The system must contain a single molecule to be decoupled! "
                                "Use the 'BioSimSpace.Align.Decouple' function to mark a molecule"
                                " to be decoupled.")

        # Validate the user specified molecular dynamics engine.
        self._exe = None
        if engine is not None:
            if not isinstance(engine, str):
                raise TypeError("'engine' must be of type 'str'.")

            # Strip whitespace from engine and convert to upper case.
            engine = engine.replace(" ", "").upper()

            # Check that the engine is supported.
            if engine not in self._engines:
                raise ValueError("Unsupported molecular dynamics engine '%s'. "
                                 "Supported engines are: %r." % ", ".join(self._engines))

            # Make sure GROMACS is installed if GROMACS engine is selected.
            if engine == "GROMACS":
                if _gmx_exe is None:
                    raise _MissingSoftwareError("Cannot use GROMACS engine as GROMACS is not installed!")
                self._exe = _gmx_exe

            elif engine == "AMBER":
                # Find a molecular dynamics engine and executable.
                engines, exes = _find_md_engines(system, protocol, engine, gpu_support)
                if not exes:
                    raise _MissingSoftwareError("Cannot use AMBER engine as AMBER is not installed!")
                elif len(exes) > 1:
                    _warnings.warn(f"Multiple AMBER engines were found. Proceeding with {exes[0]}...")
                self._exe = exes[0]

                if self._protocol.getPerturbationType() != "full":
                    raise NotImplementedError("AMBER currently only supports the 'full' perturbation "
                                              "type. Please use engine='SOMD' when running multistep "
                                              "perturbation types.")
        else:
            # Use SOMD as a default.
            engine = "SOMD"

        # Set the engine.
        self._engine = engine

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool'.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool'.")
        self._show_errors = show_errors

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")
        self._property_map = property_map

        # Create fake instance methods for 'analyse'. These
        # pass instance data through to the staticmethod versions.
        self.analyse = self._analyse
        # self.analyse_all_repeats = self._analyse_all_repeats

        # Initialise the process.
        self._initialise_process(self._system)

    def start(self):
        """Start the simulation."""
        if self._setup_only:
            _warnings.warn("No process exists! Object created in 'setup_only' mode.")
        else:
            self._process.start()

    def wait(self):
        """Wait for the simulation to finish."""
        if self._setup_only:
            _warnings.warn("No processes exist! Object created in 'setup_only' mode.")
        else:
            self._process.wait()

    def kill(self):
        """Kill the process."""
        self._process.kill()

    def workDir(self):
        """Return the working directory.

           Returns
           -------

           work_dir : str
               The path of the working directory.
        """
        return self._work_dir

    @staticmethod
    def analyse(work_dir, rest_type='Boresch', 
                append_to_lig_selection="resname LIG and not name H*",
                append_to_recept_selection="protein and not name H*",
                cutoff=10): # In Angstrom
        """Analyse existing trajectory from a simulation working directory and
        select restraints which best mimic the strongest receptor-ligand
        interactions.

           Parameters
           ----------

           work_dir : str
               The working directory for the simulation.

           rest_type: str
               The type of restraints to select (currently only Boresch is available).
               Default is Boresch.

           append_to_lig_selection: str
               Appends the supplied string to the default atom selection which chooses
               the atoms in the ligand to consider as potential anchor points. The default
               atom selection is f'resname {ligand_resname} and not name H*'. Uses the
               mdanalysis atom selection language. For example, 'not name O*' will result
               in an atom selection of f'resname {ligand_resname} and not name H* and not 
               name O*'.

           append_to_recept_selection: str
               Appends the supplied string to the default atom selection which chooses
               the atoms in the receptor to consider as potential anchor points. The default
               atom selection is 'protein and not name H*'. Uses the
               mdanalysis atom selection language. For example, 'not name N*' will result
               in an atom selection of 'protein and not name H* and not 
               name N*'.

           cutoff: float
               The greatest distance between ligand and receptor anchor atoms, in
               Angstrom. Receptor anchors further than cutoff Angstroms from the closest
               ligand anchors will not be included in the search for potential anchor points.

           Returns
           -------

           restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
               The restraints of `rest_type` which best mimic the strongest receptor-ligand
               interactions.

           normal_frame : :class:`System <BioSimSpace._SireWrappers.System>`
               The configuration of the system with the lowest restraint energy
               observed during the final percent_traj % of the trajectory.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        # TODO Implement checks on other inputs

        # TODO Implement analysis. Outline of functions required

        def getLowVariancePairs():
            pass

        def getBoreschAnchors():
            pass

        # TODO Implement restrainedDOF class which holds relevant DOF? Implement in different module?
        # Then do something like:
        # low_var_pairs = getLowVariancePairs(...)
        # anchors = getBoreschAnchors(low_var_pairs, ...)
        # restrainedDOF_obj = restrainedDOF(anchors=anchors, rest_type="Boresch", ...)
        # best_restraints = restrainedDOF_obj.getOptimalParams(...)
        # 
        # Could also have e.g. restrainedDOF_obj.plot()

    def _analyse(self, rest_type='Boresch', 
                lig_selection="resname LIG and not name H*",
                recept_selection="protein and not name H*",
                cutoff=10): # In Angstrom
        """Analyse trajectory and select restraints which best mimic strongest receptor-
        ligand interactions.

           Parameters
           ----------

           rest_type: str
               The type of restraints to select (currently only Boresch is available).
               Default is Boresch.

           append_to_lig_selection: str
               Appends the supplied string to the default atom selection which chooses
               the atoms in the ligand to consider as potential anchor points. The default
               atom selection is f'resname {ligand_resname} and not name H*'. Uses the
               mdanalysis atom selection language. For example, 'not name O*' will result
               in an atom selection of f'resname {ligand_resname} and not name H* and not 
               name O*'.

           append_to_recept_selection: str
               Appends the supplied string to the default atom selection which chooses
               the atoms in the receptor to consider as potential anchor points. The default
               atom selection is 'protein and not name H*'. Uses the
               mdanalysis atom selection language. For example, 'not name N*' will result
               in an atom selection of 'protein and not name H* and not 
               name N*'.

           cutoff: float
               The greatest distance between ligand and receptor anchor atoms, in
               Angstrom. Receptor anchors further than cutoff Angstroms from the closest
               ligand anchors will not be included in the search for potential anchor points.

           Returns
           -------

           restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
               The restraints of `rest_type` which best mimic the strongest receptor-ligand
               interactions.

           normal_frame : :class:`System <BioSimSpace._SireWrappers.System>`
               The configuration of the system with the lowest restraint energy
               observed during the trajectory.
        """

        # Return the result of calling the staticmethod, passing in the working
        # directory of this object.
        return RestraintSearch.analyse(self._work_dir, rest_type=rest_type, 
                lig_selection=lig_selection,
                recept_selection=recept_selection,
                cutoff=cutoff) 


    def _initialise_process(self, system):
        """Internal helper function to initialise the process.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

        """
        pass
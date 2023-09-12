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

"""Functionality for relative free-energy simulations."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Relative", "getData"]

import copy as _copy
import json as _json
import math as _math
import numpy as _np
import os as _os
import pandas as _pd
import pathlib as _pathlib
import pyarrow.parquet as _pq
import re as _re
import shutil as _shutil
import subprocess as _subprocess
import sys as _sys
import warnings as _warnings
import zipfile as _zipfile
import alchemlyb as _alchemlyb

from pytest import approx
from scipy.constants import proton_mass
from scipy.constants import physical_constants
hydrogen_amu = proton_mass/(physical_constants["atomic mass constant"][0])

from .._Utils import _assert_imported, _have_imported, _try_import

# alchemlyb isn't available for all variants of Python that we support, so we
# need to try_import it.
_alchemlyb = _try_import("alchemlyb")

if _have_imported(_alchemlyb):
    import logging as _logging

    # Silence pymbar warnings on startup.
    _logger = _logging.getLogger("pymbar")
    _logger.setLevel(_logging.ERROR)

    # Handle alchemlyb MBAR API changes.
    try:
        from alchemlyb.estimators import AutoMBAR as _AutoMBAR
    except ImportError:
        from alchemlyb.estimators import MBAR as _AutoMBAR
    from alchemlyb.estimators import TI as _TI
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol
    from alchemlyb.parsing.amber import extract_dHdl as _amber_extract_dHdl
    from alchemlyb.parsing.amber import extract_u_nk as _amber_extract_u_nk
    from alchemlyb.parsing.gmx import extract_dHdl as _gmx_extract_dHdl
    from alchemlyb.parsing.gmx import extract_u_nk as _gmx_extract_u_nk
    from alchemlyb.preprocessing.subsampling import (
        equilibrium_detection as _equilibrium_detection,
    )
    from alchemlyb.preprocessing.subsampling import (
        statistical_inefficiency as _statistical_inefficiency,
    )
    from alchemlyb.preprocessing.subsampling import slicing as _slicing
    from alchemlyb.preprocessing.subsampling import decorrelate_u_nk, decorrelate_dhdl
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol
    from alchemlyb.postprocessors.units import kJ2kcal as _kJ2kcal
    from alchemlyb.postprocessors.units import R_kJmol as _R_kJmol

from .._Utils import _assert_imported, _have_imported, _try_import

# alchemlyb isn't available for all variants of Python that we support, so we
# need to try_import it.
_alchemlyb = _try_import("alchemlyb")

if _have_imported(_alchemlyb):
    import logging as _logging

    # Silence pymbar warnings on startup.
    _logger = _logging.getLogger("pymbar")
    _logger.setLevel(_logging.ERROR)

    # Handle alchemlyb MBAR API changes.
    try:
        from alchemlyb.estimators import AutoMBAR as _AutoMBAR
    except ImportError:
        from alchemlyb.estimators import MBAR as _AutoMBAR
    from alchemlyb.estimators import TI as _TI
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol
    from alchemlyb.parsing.amber import extract_dHdl as _amber_extract_dHdl
    from alchemlyb.parsing.amber import extract_u_nk as _amber_extract_u_nk
    from alchemlyb.parsing.gmx import extract_dHdl as _gmx_extract_dHdl
    from alchemlyb.parsing.gmx import extract_u_nk as _gmx_extract_u_nk
    from alchemlyb.preprocessing.subsampling import (
        equilibrium_detection as _equilibrium_detection,
    )
    from alchemlyb.preprocessing.subsampling import (
        statistical_inefficiency as _statistical_inefficiency,
    )
    from alchemlyb.preprocessing.subsampling import slicing as _slicing
    from alchemlyb.preprocessing.subsampling import decorrelate_u_nk, decorrelate_dhdl
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol
    from alchemlyb.postprocessors.units import kJ2kcal as _kJ2kcal
    from alchemlyb.postprocessors.units import R_kJmol as _R_kJmol

from sire.legacy.Base import getBinDir as _getBinDir
from sire.legacy.Base import getShareDir as _getShareDir

from .. import _gmx_exe, _gmx_version
from .. import _is_notebook
from .. import _isVerbose
from .._Exceptions import AnalysisError as _AnalysisError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._Exceptions import IncompatibleError as _IncompatibleError
from .._SireWrappers import System as _System
from .._Utils import cd as _cd
from .. import Process as _Process
from .. import Protocol as _Protocol
from .. import Types as _Types
from .. import Units as _Units
from .. import _Utils

from BioSimSpace.MD._md import _find_md_engines

if _is_notebook:
    from IPython.display import FileLink as _FileLink

# Check that the analyse_freenrg script exists.
if _sys.platform != "win32":
    _analyse_freenrg = _os.path.join(_getBinDir(), "analyse_freenrg")
else:
    _analyse_freenrg = _os.path.join(
        _os.path.normpath(_getShareDir()), "scripts", "analyse_freenrg.py"
    )
if not _os.path.isfile(_analyse_freenrg):
    raise _MissingSoftwareError(
        "Cannot find free energy analysis script in expected location: '%s'"
        % _analyse_freenrg
    )
if _sys.platform == "win32":
    _analyse_freenrg = "%s %s" % (
        _os.path.join(_os.path.normpath(_getBinDir()), "sire_python.exe"),
        _analyse_freenrg,
    )


class Relative:
    """Class for configuring and running relative free-energy perturbation simulations."""

    # Create a list of supported molecular dynamics engines.
    _engines = ["AMBER", "GROMACS", "SOMD"]

    def __init__(
        self,
        system,
        protocol=None,
        work_dir=None,
        engine=None,
        gpu_support=False,
        setup_only=False,
        ignore_warnings=False,
        show_errors=True,
        extra_options=None,
        extra_lines=None,
        estimator='MBAR',
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the perturbation leg. This must contain
            a single perturbable molecule and is assumed to have already
            been equilibrated.

        protocol : :class:`Protocol.FreeEnergyMinimisation <BioSimSpace.Protocol.FreeEnergyMinimisation>`, \
                   :class:`Protocol.FreeEnergyEquilibration <BioSimSpace.Protocol.FreeEnergyEquilibration>`, \
                   :class:`Protocol.FreeEnergyProduction <BioSimSpace.Protocol.FreeEnergyProduction>`
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

        setup_only : bool
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

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Validate the input.

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        if protocol is not None:
            if isinstance(protocol, _Protocol._FreeEnergyMixin):
                self._protocol = protocol
            else:
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.FreeEnergy'"
                )
        else:
            # Use a default protocol.
            self._protocol = _Protocol.FreeEnergy()

        if not isinstance(setup_only, bool):
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            if setup_only:
                raise ValueError(
                    "A 'work_dir' must be specified when 'setup_only' is True!"
                )
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir, exist_ok=True)

        # Validate the user specified molecular dynamics engine.
        self._exe = None
        if engine is not None:
            if not isinstance(engine, str):
                raise TypeError("'engine' must be of type 'str'.")

            # Strip whitespace from engine and convert to upper case.
            engine = engine.replace(" ", "").upper()

            # Check that the engine is supported.
            if engine not in self._engines:
                raise ValueError(
                    f"Unsupported molecular dynamics engine {engine}. "
                    "Supported engines are: %r." % ", ".join(self._engines)
                )

            # Make sure GROMACS is installed if GROMACS engine is selected.
            if engine == "GROMACS":
                if _gmx_exe is None:
                    raise _MissingSoftwareError(
                        "Cannot use GROMACS engine as GROMACS is not installed!"
                    )

                # The system must have a perturbable molecule.
                if system.nPerturbableMolecules() == 0:
                    raise ValueError(
                        "The system must contain a perturbable molecule! "
                        "Use the 'BioSimSpace.Align' package to map and merge molecules."
                    )

                if self._protocol.getPerturbationType() != "full":
                    raise NotImplementedError(
                        "GROMACS currently only supports the 'full' perturbation "
                        "type. Please use engine='SOMD' when running multistep "
                        "perturbation types."
                        )
                self._exe = _gmx_exe
            elif engine == "AMBER":
                # Find a molecular dynamics engine and executable.
                engines, exes = _find_md_engines(
                    system, protocol, engine, gpu_support)
                if not exes:
                    raise _MissingSoftwareError(
                        "Cannot use AMBER engine as AMBER is not installed!")
                elif len(exes) > 1:
                    _warnings.warn(
                        f"Multiple AMBER engines were found. Proceeding with {exes[0]}...")
                self._exe = exes[0]

                # The system must have a perturbable molecule.
                if system.nPerturbableMolecules() == 0:
                    raise ValueError("The system must contain a perturbable molecule! "
                                     "Use the 'BioSimSpace.Align' package to map and merge molecules.")

                if self._protocol.getPerturbationType() != "full":
                    raise NotImplementedError("AMBER currently only supports the 'full' perturbation "
                                              "type. Please use engine='SOMD' when running multistep "
                                              "perturbation types.")
        else:
            # Use SOMD as a default.
            engine = "SOMD"

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError(
                    "The system must contain a single perturbable molecule! "
                    "Use the 'BioSimSpace.Align' package to map and merge molecules."
                )

        # Set the engine.
        self._engine = engine

        # Validate the protocol.
        if protocol is not None:
            from ..Protocol._free_energy_mixin import _FreeEnergyMixin

            if isinstance(protocol, _FreeEnergyMixin):
                if engine == "SOMD" and not isinstance(
                    protocol, _Protocol.FreeEnergyProduction
                ):
                    raise ValueError(
                        "Currently SOMD only supports protocols of type 'BioSimSpace.Protocol.FreeEnergyProduction'"
                    )

                self._protocol = protocol
            else:
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.FreeEnergy'"
                )
        else:
            # Use a default protocol.
            self._protocol = _Protocol.FreeEnergy()

        # Check that multi-step perturbation isn't specified if GROMACS is the chosen engine.
        if engine == "GROMACS":
            if self._protocol.getPerturbationType() != "full":
                raise NotImplementedError(
                    "GROMACS currently only supports the 'full' perturbation "
                    "type. Please use engine='SOMD' when running multistep "
                    "perturbation types."
                )

        if not isinstance(setup_only, bool):
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        if work_dir is None and setup_only:
            raise ValueError(
                "A 'work_dir' must be specified when 'setup_only' is True!"
            )

        # Create the working directory.
        self._work_dir = _Utils.WorkDir(work_dir)

        # HMR check
        # by default, if the timestep is greater than 4 fs this should be true.
        if self._protocol.getHmr():

            # Set the expected HMR factor.
            hmr_factor = self._protocol.getHmrFactor()
            if hmr_factor == "auto":
                if self._engine == "AMBER" or self._engine == "GROMACS":
                    hmr_factor = 3
                elif self._engine == "SOMD":
                    self._extra_options["hydrogen mass repartitioning factor"] = "1.5"
            else:
                if self._engine == "SOMD":
                    self._extra_options["hydrogen mass repartitioning factor"] = str(
                        hmr_factor)

            # Extract the molecules to which HMR applies.
            # set auto based on engine
            if self._protocol.getHmrWater() == "auto" or self._protocol.getHmrWater() == False:
                water = "no"
                # water_mols = system.getWaterMolecules()
                # molecules = system - water_mols
                molecules = self._system.search("not water", property_map)
            elif self._protocol.getHmrWater() == True:
                water = "yes"
                molecules = self._system.getMolecules()

            # TODO update search to whole system and mass search term once this is updated in the Sire API
            # use length of no of H to find average mass of H and if correctly repartitioned.
            # mass = (molecules.search("{element H}[0]"))[0]._sire_object.evaluate().mass().value()
            # Currently checking the mass (in g per mol) of the first H only.

            for mol in molecules:
                try:
                    h = mol.search("{element H}[0]")
                    # mass = h._sire_object.property(mass_prop).value()
                    mass = h[0]._sire_object.evaluate().mass(
                        property_map).value()
                    found_h = True
                except:
                    found_h = False
                if found_h:
                    break

            if not found_h:
                # in some cases, may not be able to find the mass based on element
                # if the only molecule is perturbable. Search for atom name with H.
                # need mass0 as this is one of the perturbable molecule properties.
                mass_prop = property_map.get("mass0", "mass0")
                molecule = self._system.getPerturbableMolecules()[0]
                for atom in molecule.getAtoms():
                    # check if the atom is a H atom
                    if atom.name().startswith("H"):
                        mass = atom._sire_object.property(mass_prop).value()
                        # check in case the mass at 0 is 0 as perturbable molecule.
                        if mass > 1:
                            found_h = True
                    if found_h:
                        break

            # error if cant find a H mass
            if not found_h:
                raise _IncompatibleError(
                    "Can't find the mass of a H in the system.")

            # Check that the mass matches what is expected.
            if self._engine == "SOMD":
                # values should be in amu
                if mass == approx(hydrogen_amu, rel=1e-2):
                    repartition = False
                else:
                    raise _IncompatibleError(
                        "Please do not pass an already repartitioned system in for use with SOMD.")

            # check for amber or gromacs repartitioning
            elif self._engine == "AMBER" or self._engine == "GROMACS":
                # check if the system has been repartitioned at all. If not, repartition.
                if mass == approx(hydrogen_amu, rel=1e-2):
                    repartition = True
                # check if system as been repartitioned with the decided factor.
                elif mass == approx(hmr_factor * hydrogen_amu, rel=1e-2):
                    repartition = False
                # finally, check if the system is repartitioned with the wrong factor.
                elif mass != approx(hmr_factor * hydrogen_amu, rel=1e-2):
                    raise _IncompatibleError("""
                    The system is repartitioned at a factor different from that specified in 'hmr_factor'
                    or at the auto default for this engine (3 for AMBER and GROMACS, None for SOMD (as this is specified in the cfg file)).
                    Please pass a correctly partitioned or entirely unpartitioned system.""")

            # Repartition if necessary.
            if repartition:
                _warnings.warn(
                    f"The passed system is being repartitioned according to a factor of '{hmr_factor}'.")
                self._system.repartitionHydrogenMass(
                    factor=hmr_factor, water=water, property_map=property_map)
            else:
                if self._engine != "SOMD":
                    _warnings.warn(
                        "The passed system is already repartitioned. Proceeding without additional repartitioning.")

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool'.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool'.")
        self._show_errors = show_errors

        # Check that the estimator is either MBAR or TI.
        if not isinstance(estimator, str):
            raise ValueError("'estimator' must be of type 'str'.")
        if estimator not in ['MBAR', 'TI']:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")
        self._estimator = estimator

        # Check the extra options.
        if not isinstance(extra_options, dict):
            raise TypeError("'extra_options' must be of type 'dict'.")
        else:
            keys = extra_options.keys()
            if not all(isinstance(k, str) for k in keys):
                raise TypeError("Keys of 'extra_options' must be of type 'str'.")
        self._extra_options = extra_options

        # Check the extra lines.
        if not isinstance(extra_lines, list):
            raise TypeError("'extra_lines' must be of type 'list'.")
        else:
            if not all(isinstance(line, str) for line in extra_lines):
                raise TypeError("Lines in 'extra_lines' must be of type 'str'.")
        self._extra_lines = extra_lines

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")
        self._property_map = property_map

        # Create fake instance methods for 'analyse', 'checkOverlap',
        # and 'difference'. These pass instance data through to the
        # staticmethod versions.
        self.analyse = self._analyse
        self.checkOverlap = self._checkOverlap
        self.difference = self._difference
        self.check_overlap = self._check_overlap
        self.plot = self._plot

        # Initialise the process runner.
        self._initialise_runner(self._system)

    def run(self, serial=True):
        """
        Run the simulation.

        Parameters
        ----------

        serial : bool
            Whether to run the individual processes for the lambda windows
            in serial.
        """
        if not isinstance(serial, bool):
            raise TypeError("'serial' must be of type 'bool'.")

        if self._setup_only:
            _warnings.warn(
                "No processes exist! Object created in 'setup_only' mode.")
        else:
            self._runner.startAll(serial=serial)

    def wait(self):
        """Wait for the simulation to finish."""
        if self._setup_only:
            _warnings.warn(
                "No processes exist! Object created in 'setup_only' mode.")
        else:
            self._runner.wait()

    def kill(self, index):
        """
        Kill a process for a specific lambda window.

        Parameters
        ----------

        index : int
            The index of the lambda window.
        """
        self._runner.kill(index)

    def killAll(self):
        """Kill any running processes for all lambda windows."""

        self._runner.killAll()

    def workDir(self):
        """
        Return the working directory.

        Returns
        -------

        work_dir : str
            The path of the working directory.
        """
        return str(self._work_dir)

    def getData(self, name="data", file_link=False, work_dir=None):
        """
        Return a link to a zip file containing the data files required for
        post-simulation analysis.

        Parameters
        ----------

        name : str
            The name of the zip file.

        file_link : bool
            Whether to return a FileLink when working in Jupyter.

        work_dir : str
            The working directory for the free-energy perturbation
            simulation.

        Returns
        -------

        output : str, IPython.display.FileLink
            A path, or file link, to an archive of the process input.
        """

        if self._work_dir is None:
            raise ValueError("'work_dir' must be set!")
        else:
            if not isinstance(work_dir, str):
                raise TypeError("'work_dir' must be of type 'str'.")
            if not _os.path.isdir(work_dir):
                raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(name, str):
            raise TypeError("'name' must be of type 'str'")

        # Generate the zip file name.
        zipname = "%s.zip" % name

        # Get the current working directory.
        cwd = _os.getcwd()

        # Change into the working directory.
        with _cd(work_dir):
            # Specify the path to glob.
            glob_path = _pathlib.Path(work_dir)

            # First try SOMD data.
            files = glob_path.glob("**/gradients.dat")

            if len(files) == 0:
                files = glob_path.glob("**/[!bar]*.xvg")

                if len(files) == 0:
                    files = glob_path.glob("**/*.out")

                    if len(files) == 0:
                        raise ValueError(
                            f"Couldn't find any analysis files in '{work_dir}'"
                        )

            # Write to the zip file.
            with _zipfile.ZipFile(cwd + f"/{zipname}", "w") as zip:
                for file in files:
                    zip.write(file)

        # Return a link to the archive.
        if _is_notebook:
            if file_link:
                # Create a FileLink to the archive.
                f_link = _FileLink(zipname)

                # Set the download attribute so that JupyterLab doesn't try to open the file.
                f_link.html_link_str = (
                    f"<a href='%s' target='_blank' download='{zipname}'>%s</a>"
                )

                # Return a link to the archive.
                return f_link
            else:
                return zipname
        # Return the path to the archive.
        else:
            return zipname

    @staticmethod
    def analyse(work_dir, estimator="MBAR", method="alchemlyb", **kwargs):
        """
        Analyse existing free-energy data from a simulation working directory.

        Parameters
        ----------

        work_dir : str
            The working directory for the simulation.

        estimator : str
            The estimator to use for the free-energy analysis. ("MBAR" or "TI")

        method : str
            The method to use for the free-energy analysis. ("alchemlyb" or "native")

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : np.matrix, None
            The overlap matrix. This gives the overlap between each lambda
            window. This parameter is only computed when available for the
            specified estimator and engine, otherwise None will be returned.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")

        if not isinstance(method, str):
            raise TypeError("'method' must be of type 'str'.")
        method = method.replace("-", "").upper()

        if method == "ALCHEMLYB":
            _assert_imported("alchemlyb")

        function_glob_dict = {
            "SOMD": (Relative._analyse_somd, "**/simfile.dat"),
            "SOMD2": (Relative._analyse_somd2, "**/*.parquet"),
            "GROMACS": (Relative._analyse_gromacs, "**/[!bar]*.xvg"),
            "AMBER": (Relative._analyse_amber, "**/*.out"),
        }

        for engine, (func, mask) in function_glob_dict.items():
            glob_path = _pathlib.Path(work_dir)
            data = sorted(glob_path.glob(mask))
            if data and engine == "AMBER":
                if method != "ALCHEMLYB":
                    raise _AnalysisError(
                        f"AMBER can only use the 'alchemlyb' analysis method."
                    )
            if data and engine == "SOMD" and estimator == "TI" and method == "native":
                raise _AnalysisError(
                    "SOMD cannot use 'TI' estimator with 'native' analysis method."
                )
            if data and engine == "SOMD2":
                if method != "ALCHEMLYB":
                    raise _AnalysisError(
                        f"SOMD2 can only use the 'alchemlyb' analysis method."
                    )
            if data and engine == "GROMACS" and method == "NATIVE":
                _warnings.warn(
                    "Native GROMACS analysis cannot do MBAR/TI. BAR will be used."
                )
            if data:
                return func(work_dir, estimator, method, **kwargs)

        raise ValueError(
            "Couldn't find any SOMD, SOMD2, GROMACS or AMBER free-energy output?"
        )

    @staticmethod
    def checkOverlap(overlap, threshold=0.03):
        """
        Check the overlap of an FEP leg.

        Parameters
        ----------

        overlap : numpy.matrix
            The overlap matrix. This gives the overlap between lambda windows.

        threshold : float
            The threshold value used to check the off-diagonals.

        Returns
        -------

        is_okay : boolean
             True if the overlap is okay, False if any off-diagonals are less
             than the threshold value.

        num_low : int
            The number of off-diagonals that are less than the threshold value.

        """
        if not isinstance(overlap, _np.matrix):
            raise TypeError("'overlap' must be of type 'numpy.matrix'.")

        if not isinstance(threshold, float):
            raise TypeError("'threshold' must be of type 'float'.")
        if threshold < 0.0 or threshold > 1.0:
            raise ValueError("'threshold' must be between 0 and 1.")

        # Get all off diagonal elements.
        off_diagonal = (_np.diagonal(overlap, 1)).tolist()
        for x in (_np.diagonal(overlap, -1)).tolist():
            off_diagonal.append(x)

        # Check if any off diagonals are less than the threshold value.
        num_low = 0
        is_okay = False
        for od in off_diagonal:
            if od < threshold:
                num_low += 1
        if num_low > 0:
            _warnings.warn(
                f"Overlap is poor: {num_low} off-diagonals are less than {threshold}."
            )
        else:
            is_okay = True

        return is_okay, num_low

    @staticmethod
    def difference(pmf, pmf_ref):
        """
        Compute the relative free-energy difference between two perturbation
        legs.

        Parameters
        ----------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF) of interest. The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        pmf_ref : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The reference potential of mean force (PMF). The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The relative free-energy difference and its associated error.
        """

        if not isinstance(pmf, list):
            raise TypeError("'pmf' must be of type 'list'.")
        if not isinstance(pmf_ref, list):
            raise TypeError("'pmf_ref' must be of type 'list'.")

        for rec in pmf:
            if not isinstance(rec, tuple):
                raise TypeError(
                    "'pmf1' must contain tuples containing lambda "
                    "values and the associated free-energy and error."
                )
            else:
                if len(rec) != 3:
                    raise ValueError(
                        "Each tuple in 'pmf1' must contain three items: "
                        "a lambda value and the associated free energy "
                        "and error."
                    )
                for val in rec[1:]:
                    if not isinstance(val, _Types.Energy):
                        raise TypeError(
                            "'pmf' must contain 'BioSimSpace.Types.Energy' types."
                        )

        for rec in pmf_ref:
            if not isinstance(rec, tuple):
                raise TypeError(
                    "'pmf_ref' must contain tuples containing lambda "
                    "values and the associated free-energy and error."
                )
            else:
                if len(rec) != 3:
                    raise ValueError(
                        "Each tuple in 'pmf_ref' must contain three items: "
                        "a lambda value and the associated free energy "
                        "and error."
                    )
                for val in rec[1:]:
                    if not isinstance(val, _Types.Energy):
                        raise TypeError(
                            "'pmf_ref' must contain 'BioSimSpace.Types.Energy' types."
                        )

        # Work out the difference in free energy.
        free_energy = (pmf[-1][1] - pmf[0][1]) - (pmf_ref[-1][1] - pmf_ref[0][1])

        # Propagate the errors. (These add in quadrature.)

        # Measure.
        error0 = _math.sqrt(
            (pmf[-1][2].value() * pmf[-1][2].value())
            + (pmf[0][2].value() * pmf[0][2].value())
        )

        # Reference.
        error1 = _math.sqrt(
            (pmf_ref[-1][2].value() * pmf_ref[-1][2].value())
            + (pmf_ref[0][2].value() * pmf_ref[0][2].value())
        )

        # Error for free-energy difference.
        error = (
            _math.sqrt((error0 * error0) + (error1 * error1))
            * _Units.Energy.kcal_per_mol
        )

        return (free_energy, error)

    @staticmethod
    def _get_data(files, temperatures, engine, estimator):
        """
        files : list(pathlib.Path)
            List of files for all lambda values to analyse. Should be sorted.

        temperatures : list(float)
            List of temperatures at which the simulation was carried out at for each lambda window.
            Index of the temperature value should match it's corresponding lambda window index in files.

        lambdas : list(float)
            Sorted list of lambda values used for the simulation.

        engine : str
            Engine with which the simulation was run.

        estimator : str
            The estimator to use for the analysis. Options are "MBAR" or "TI".

        Returns
        -------

        data : list(pandas.DataFrame)
            A list of dataframes containing the data for each lambda window.
        """

        if not isinstance(files, (tuple, list)):
            raise TypeError("'files' must be of type 'list' or 'tuple'.")
        if not all(isinstance(x, _pathlib.Path) for x in files):
            raise TypeError("'files' must be a list of 'pathlib.Path' types.")

        if not isinstance(temperatures, (tuple, list)):
            raise TypeError("'temperatures' must be of type 'list' or 'tuple'.")
        if not all(isinstance(x, float) for x in temperatures):
            raise TypeError("'temperatures' must be a list of 'float' types.")

        if not isinstance(engine, str):
            raise TypeError("'engine' must be of type 'str'.")
        engine = engine.replace(" ", "").upper()
        if not engine in Relative._engines_analysis:
            raise ValueError(
                f"Unsupported engine '{engine}'. Options are: {', '.join(Relative._engines_analysis)}"
            )

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        estimator = estimator.replace(" ", "").upper()
        if not estimator in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        if estimator == "MBAR":
            is_mbar = True
        else:
            is_mbar = False

        from functools import partial

        function_dict = {
            "SOMD": partial(Relative._somd_extract, estimator=estimator),
            "SOMD2": partial(Relative._somd2_extract, estimator=estimator),
            "GROMACS": _gmx_extract_u_nk if is_mbar else _gmx_extract_dHdl,
            "AMBER": _amber_extract_u_nk if is_mbar else _amber_extract_dHdl,
        }

        # Extract the data.
        func = function_dict[engine]
        try:
            data = [func(file, T=temp) for file, temp in zip(files, temperatures)]
        except Exception as e:
            msg = "Could not extract the data from the provided files!"
            if _isVerbose():
                raise _AnalysisError(msg) from e
            else:
                raise _AnalysisError(msg) from None

        return data

    @staticmethod
    def _get_u_nk(files, temperatures, engine):
        """
        Get the u_nk dataframes for MBAR analysis.

        Parameters
        ----------

        files : list(pathlib.Path)
            A list of data files.

        temperatures : list(float)
            A list of temperatures.

        engine : str
            The simulation engine used to generate the data.

        Returns
        -------

        u_nk : list(pandas.DataFrame)
            A list of dataframes containing the u_nk data for each lambda window.
        """
        return Relative._get_data(files, temperatures, engine, "MBAR")

    @staticmethod
    def _get_dh_dl(files, temperatures, engine):
        """
        Get the dh_dl dataframes for TI analysis.

        Parameters
        ----------

        files : list(pathlib.Path)
            A list of data files.

        temperatures : list(float)
            A list of temperatures.

        engine : str
            The simulation engine used to generate the data.

        Returns
        -------

        dh_dl : list(pandas.DataFrame)
            A list of dataframes containing the u_nk data for each lambda window.
        """
        return Relative._get_data(files, temperatures, engine, "TI")

    def _analyse(self, estimator="MBAR"):
        """
        Analyse free-energy data for this object.

        Parameters
        ----------

        estimator : str
            The estimator to use for the free-energy analysis. ("MBAR" or "TI")

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : np.matrix, None
            The overlap matrix. This gives the overlap between each lambda
            window. This parameter is only computed when available for the
            specified estimator and engine, otherwise None will be returned.
        """

        # Return the result of calling the staticmethod, passing in the working
        # directory of this object and the specified estimator.
        return Relative.analyse(str(self._work_dir), estimator=estimator)

    @staticmethod
    def _somd_extract(simfile, T=None, estimator="MBAR"):
        """
        Extract required data from a SOMD output file (simfile.dat).

        Parameters
        ----------

        simfile : str
            Path to the simfile.dat file.

        T : float
            Temperature in Kelvin at which the simulations were performed;
            needed to generated the reduced potential (in units of kT).

        estimator : str
            The estimator that the returned data will be used with. This can
            be either 'MBAR' or 'TI'.

        Returns
        -------

        data : pandas.DataFrame
            Either: Reduced potential for each alchemical state (k) for each
            frame (n) for MBAR, or dH/dl as a function of time for this lambda
            window for TI.
        """

        if not isinstance(simfile, _pathlib.Path):
            raise TypeError("'simfile' must be of type 'pathlib.Path'.")
        if not _os.path.isfile(simfile):
            raise ValueError("'simfile' doesn't exist!")

        if not isinstance(T, float):
            raise TypeError("'T' must be of type 'float'.")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        if estimator.replace(" ", "").upper() not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        # Flag that we're using MBAR.
        if estimator == "MBAR":
            is_mbar = True
        else:
            is_mbar = False

        # For dhdl need to consider the temperature, as the gradient is in
        # kcal/mol in the simfile.dat.
        if not is_mbar:
            k_b = _R_kJmol * _kJ2kcal
            beta = 1 / (k_b * T)

        # Flags to check if we have found the required records.
        found_lambda = False
        found_array = False
        found_time = False

        # Process the file.
        with open(simfile, "r") as f:
            # Terms to search for in the record lines.
            start_w = "#Generating lambda is"
            start_a = "#Alchemical array is"
            start_t = " and "
            end_t = " ps"
            # Read the file line-by-line.
            for line in f.readlines():
                if start_w in line:
                    lambda_win = float(line.replace(start_w, "").strip())
                    if lambda_win is not None:
                        found_lambda = True
                if start_a in line:
                    lambda_array = (
                        (line.replace(start_a, ""))
                        .strip()
                        .replace("(", "")
                        .replace(")", "")
                        .replace(" ", "")
                    ).split(",")
                    lambda_array = [float(lam) for lam in lambda_array]
                    if lambda_array is not None:
                        found_array = True
                if start_t and end_t in line:
                    sim_length = float(
                        ((line.split(start_t)[1]).split(end_t)[0]).strip()
                    )
                    if sim_length is not None:
                        found_time = True
                # All records found, break the loop.
                if found_lambda:
                    if found_array:
                        if found_time:
                            break

        if not found_lambda:
            raise ValueError(
                f"The lambda window was not detected in the SOMD output file: {file}"
            )

        if not found_array:
            raise ValueError(
                f"The lambda array was not detected in the SOMD output file: {file}"
            )

        if not found_time:
            raise ValueError(
                f"The simulation time was not detected in the SOMD output file: {file}"
            )

        # TODO: get header from the file instead of like this.
        header = [
            "step",
            "potential_kcal/mol",
            "gradient_kcal/mol",
            "forward_Metropolis",
            "backward_Metropolis",
        ]
        header.extend(lambda_array)

        file_df = _pd.read_fwf(
            simfile, skipinitialspace=True, skiprows=13, header=None, names=header
        )

        time_step = sim_length / len(file_df["step"])
        time_rows = _np.arange(0, len(file_df["step"]), 1)
        time = _np.arange(0, sim_length, time_step)

        # For MBAR, results in list of lists where each list is the 0 to 1
        # window values that lambda value. For TI, it is a list of gradients
        # at that lambda.
        if is_mbar:
            results = (
                file_df.iloc[:, 5:].subtract(file_df[lambda_win], axis=0).to_numpy()
            )
        else:
            gradient = file_df["gradient_kcal/mol"].to_numpy()
            results = gradient * beta

        # Turn into a dataframe that can be processed by alchemlyb.
        if is_mbar:
            df = _pd.DataFrame(
                results,
                columns=_np.array(lambda_array, dtype=_np.float64),
                index=_pd.MultiIndex.from_arrays(
                    [time, _np.repeat(lambda_win, len(time))], names=["time", "lambdas"]
                ),
            )
        else:
            df = _pd.DataFrame(
                results,
                columns=["fep"],
                index=_pd.MultiIndex.from_arrays(
                    [time, _np.repeat(lambda_win, len(time))],
                    names=["time", "fep-lambdas"],
                ),
            )
        df.attrs["temperature"] = T
        df.attrs["energy_unit"] = "kT"

        return df

    @staticmethod
    def _somd2_extract(parquet_file, T=None, estimator="MBAR"):
        """
        Extract required data from a SOMD2 output file (parquet file).

        Parameters
        ----------

        parquet_file : str
            Path to the parquet file.

        T : float
            Temperature in Kelvin at which the simulations were performed;
            needed to generated the reduced potential (in units of kT).

        estimator : str
            The estimator that the returned data will be used with. This can
            be either 'MBAR' or 'TI'.

        Returns
        -------

        data : pandas.DataFrame
            Either: Reduced potential for each alchemical state (k) for each
            frame (n) for MBAR, or dH/dl as a function of time for this lambda
            window for TI.
        """

        if not isinstance(parquet_file, _pathlib.Path):
            raise TypeError("'parquet_file' must be of type 'pathlib.Path'.")
        if not _os.path.isfile(parquet_file):
            raise ValueError("'parquet_file' doesn't exist!")

        if not isinstance(T, float):
            raise TypeError("'T' must be of type 'float'.")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        if estimator.replace(" ", "").upper() not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        # Flag that we're using MBAR.
        if estimator == "MBAR":
            is_mbar = True
        else:
            is_mbar = False

        # Beta factor for computing reduced potentials and gradients.
        k_b = _R_kJmol * _kJ2kcal
        beta = 1 / (k_b * T)

        # Try to read the file.
        try:
            table = _pq.read_table(parquet_file)
        except:
            raise IOError(f"Could not read the SOMD2 parquet file: {parquet_file}")

        # Try to extract the metadata.
        try:
            metadata = _json.loads(table.schema.metadata["somd2".encode()])
            temperature = float(metadata["temperature"])
        except:
            raise IOError(
                f"Could not read the SOMD2 metadta from parquet file: {parquet_file}"
            )

        # Validate metadata required by all analysis methods.
        try:
            lam = float(metadata["lambda"])
        except:
            raise ValueError("Parquet metadata does not contain lambda value.")
        try:
            lambda_array = metadata["lambda_array"]
        except:
            raise ValueError("Parquet metadata does not contain lambda array.")

        # Make sure that the temperature is correct.
        if not T == temperature:
            raise ValueError(
                f"The temperature in the parquet metadata '{temperature:%.3f}' "
                f"does not match the specified temperature '{T:%.3f}'."
            )

        # Convert to a pandas dataframe.
        df = table.to_pandas()

        if is_mbar:
            df = df[[str(x) for x in lambda_array]].copy()
            lambdas = len(df.index) * [lam]
            multiindex = _pd.MultiIndex.from_arrays(
                [df.index, lambdas], names=["time", "lambdas"]
            )
            df.index = multiindex
            df = beta * df.subtract(df[str(lam)], axis=0)

            # Set the temperature and energy unit.
            df.attrs["temperature"] = T
            df.attrs["energy_unit"] = "kT"

            # Convert column index to float.
            df.columns = df.columns.astype(float)

            return df.dropna()

        else:
            columns_lambdas = df.columns[
                _pd.to_numeric(df.columns, errors="coerce").to_series().notnull()
            ]

            if len(columns_lambdas) > 3 and lambda_array is None:
                raise ValueError(
                    "More than 3 lambda values in the dataframe but no lambda array provided?"
                )
            try:
                lam_below = max(
                    [
                        float(lambda_val)
                        for lambda_val in columns_lambdas
                        if float(lambda_val) < lam
                    ]
                )
            except ValueError:
                lam_below = None
            try:
                lam_above = min(
                    [
                        float(lambda_val)
                        for lambda_val in columns_lambdas
                        if float(lambda_val) > lam
                    ]
                )
            except ValueError:
                lam_above = None

            # Compute gradient using finite differences.
            if lam_below is None:
                double_incr = (lam_above - lam) * 2
                grad = (df[str(lam_above)] - df[str(lam)]) * 2 / double_incr
                back_m = _np.exp(beta * (df[str(lam_above)] - df[str(lam)]))
                forward_m = _np.exp(-1 * beta * (df[str(lam_above)] - df[str(lam)]))
            elif lam_above is None:
                double_incr = (lam - lam_below) * 2
                grad = (df[str(lam)] - df[str(lam_below)]) * 2 / double_incr
                back_m = _np.exp(-1 * beta * (df[str(lam_below)] - df[str(lam)]))
                forward_m = _np.exp(beta * (df[str(lam_below)] - df[str(lam)]))
            else:
                double_incr = lam_above - lam_below
                grad = (df[str(lam_above)] - df[str(lam_below)]) / double_incr
                back_m = _np.exp(beta * (df[str(lam)] - df[str(lam_below)]))
                forward_m = _np.exp(beta * (df[str(lam)] - df[str(lam_above)]))

            grad.name = "gradient"
            back_m.name = "backward_mc"
            forward_m.name = "forward_mc"

            if lambda_array is not None:
                df[[str(i) for i in lambda_array]] = df[
                    [str(i) for i in lambda_array]
                ].apply(lambda x: x * -1 * beta)

            df = _pd.concat(
                [
                    df,
                    _pd.DataFrame(grad),
                    _pd.DataFrame(back_m),
                    _pd.DataFrame(forward_m),
                ],
                axis=1,
            )

            try:
                time = list(df["time"])
            except KeyError:
                time = list(df.index)

            lambdas = len(time) * [lam]

            # Create a multi-index from the two lists
            multi_index = _pd.MultiIndex.from_tuples(
                zip(time, lambdas), names=["time", "fep-lambdas"]
            )

            # Get the gradients.
            grads = list(df["gradient"])

            # Create a DataFrame with the multi-index
            df = _pd.DataFrame({"fep": grads}, index=multi_index)

            # Set the temperature and energy unit.
            df.attrs["temperature"] = T
            df.attrs["energy_unit"] = "kT"

            return df.dropna()

    @staticmethod
    def _preprocess_data(data, estimator, **kwargs):
        """
        Preprocess FEP simulation data.

        Parameters
        ----------

        data : list
                List of extracted dHdl or u_nk data.

        estimator : str
               The estimator ('MBAR' or 'TI') used.

        Returns
        -------

        processed_data : pandas.DataFrame
            Dataframe of dHdl or u_nk data processed using automated equilibration
            detection followed by statistical inefficiency.
        """

        if not isinstance(data, (list, _pd.DataFrame)):
            raise TypeError("'data' must be of type 'list' or 'pandas.DataFrame'.")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        if not estimator.replace(" ", "").upper() in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        # Assign defaults in case not passed via kwargs.
        auto_eq = False
        stat_ineff = False
        truncate = False
        truncate_keep = "start"

        # Parse kwargs.
        for key, value in kwargs.items():
            key = key.replace(" ", "").replace("_", "").upper()
            if key == "AUTOEQUILIBRATION":
                auto_eq = value
            if key == "STATISTICALINEFFICIENCY":
                stat_ineff = value
            if key == "TRUNCATEPERCENTAGE":
                truncate = value
            if key == "TRUNCATEKEEP":
                truncate_keep = value

        # First truncate data.
        raw_data = data
        if truncate:
            # Get the upper and lower bounds for truncate.
            data_len = len(data[0])
            data_step = round((data[0].index[-1][0] - data[0].index[-2][0]), 1)
            data_kept = data_len * (truncate / 100)
            data_time = data_kept * data_step
            if truncate_keep == "start":
                truncate_lower = 0
                truncate_upper = data_time - data_step
            if truncate_keep == "end":
                truncate_lower = (data_len * data_step) - data_time
                truncate_upper = (data_len * data_step) - data_step

            try:
                data = [
                    _slicing(i, lower=truncate_lower, upper=truncate_upper)
                    for i in raw_data
                ]
            except:
                _warnings.warn("Could not truncate data.")
                data = raw_data
        else:
            data = raw_data

        # If using auto equilibration, do we want to remove burn in?
        # Perform a statistical inefficiency check afterwards.
        if stat_ineff:
            if estimator == "MBAR":
                decorrelated_data = [
                    decorrelate_u_nk(i, method="dE", remove_burnin=auto_eq)
                    for i in data
                ]

            elif estimator == "TI":
                decorrelated_data = [
                    decorrelate_dhdl(i, remove_burnin=auto_eq) for i in data
                ]

            sampled_data = decorrelated_data

            for i in decorrelated_data:
                if len(i.iloc[:, 0]) < 50:
                    _warnings.warn(
                        "Less than 50 samples as a result of preprocessing. No preprocessing will be performed."
                    )
                    sampled_data = data
        else:
            # Need stats ineff for auto eq to run as well.
            if auto_eq:
                _warnings.warn(
                    "Auto equilibration can only be detected if statistical inefficiency is run as well."
                )
            sampled_data = data

        # Concatanate in alchemlyb format.
        processed_data = _alchemlyb.concat(sampled_data)

        return processed_data

    @staticmethod
    def _analyse_internal(files, temperatures, lambdas, engine, estimator, **kwargs):
        """
        Analyse existing free-energy data using MBAR and the alchemlyb library.

        Parameters
        ----------

        files : list(pathlib.Path)
            List of files for all lambda values to analyse. Should be sorted.

        temperatures : list(float)
            List of temperatures at which the simulation was carried out at for each lambda window.
            Index of the temperature value should match it's corresponding lambda window index in files.

        lambdas : list(float)
            Sorted list of lambda values used for the simulation.

        engine : str
            Engine with which the simulation was run.

        estimator : str
            The estimator to use for the analysis. Options are "MBAR" or "TI".

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : numpy.matrix, None
            The overlap matrix. This gives the overlap between each lambda
            window. Returns None if overlap isn't supported for the chosen
            estimator or engine.
        """

        if not isinstance(files, (tuple, list)):
            raise TypeError("'files' must be of type 'list' or 'tuple'.")
        if not all(isinstance(x, _pathlib.Path) for x in files):
            raise TypeError("'files' must be a list of 'pathlib.Path' types.")

        if not isinstance(temperatures, (tuple, list)):
            raise TypeError("'temperatures' must be of type 'list' or 'tuple'.")
        if not all(isinstance(x, float) for x in temperatures):
            raise TypeError("'temperatures' must be a list of 'float' types.")

        if not isinstance(lambdas, (tuple, list)):
            raise TypeError("'lambdas' must be of type 'list' or 'tuple'.")
        if not all(isinstance(x, float) for x in lambdas):
            raise TypeError("'lambdas' must be a list of 'float' types.")

        if not isinstance(engine, str):
            raise TypeError("'engine' must be of type 'str'.")
        if not engine.replace(" ", "").upper() in Relative._engines_analysis:
            raise ValueError(
                f"Unsupported engine '{engine}'. Options are: {', '.join(Relative._engines_analysis)}"
            )

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        estimator = estimator.replace(" ", "").upper()
        if not estimator in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        if estimator == "MBAR":
            is_mbar = True
        else:
            is_mbar = False

        # Extract the data.
        try:
            data = Relative._get_data(files, temperatures, engine, estimator)
        except Exception as e:
            msg = "Could not extract the data from the provided files!"
            if _isVerbose():
                raise _AnalysisError(msg) from e
            else:
                raise _AnalysisError(msg) from None

        # Preprocess the data.
        try:
            processed_data = Relative._preprocess_data(data, estimator, **kwargs)
        except:
            _warnings.warn("Could not preprocess the data!")
            processed_data = _alchemlyb.concat(data)

        mbar_method = None
        if is_mbar:
            # Check kwargs in case there is an mbar_method and then use this
            for key, value in kwargs.items():
                key = key.replace(" ", "").replace("_", "").upper()
                if key == "MBARMETHOD":
                    mbar_method = value

        # MBAR analysis using a specified method.
        if mbar_method:
            try:
                alchem = _AutoMBAR(method=mbar_method)
                alchem.fit(processed_data)
            except Exception as e:
                msg = f"MBAR free-energy analysis failed with {mbar_method} as mbar_method!"
                if _isVerbose():
                    raise _AnalysisError(msg) from e
                else:
                    raise _AnalysisError(msg) from None
        # Standard analysis.
        else:
            try:
                if is_mbar:
                    alchem = _AutoMBAR().fit(processed_data)
                else:
                    alchem = _TI().fit(processed_data)
            except Exception as e:
                msg = f"{estimator} free-energy analysis failed!"
                if _isVerbose():
                    raise _AnalysisError(msg) from e
                else:
                    raise _AnalysisError(msg) from None

        # Extract the data from the results.
        data = []
        # Convert the data frames to kcal/mol.
        delta_f_ = _to_kcalmol(alchem.delta_f_)
        d_delta_f_ = _to_kcalmol(alchem.d_delta_f_)
        for lambda_, t in zip(lambdas, temperatures):
            x = lambdas.index(lambda_)
            mbar_value = delta_f_.iloc[0, x]
            mbar_error = d_delta_f_.iloc[0, x]

            # Append the data.
            data.append(
                (
                    lambda_,
                    (mbar_value) * _Units.Energy.kcal_per_mol,
                    (mbar_error) * _Units.Energy.kcal_per_mol,
                )
            )

        if is_mbar:
            return (data, _np.matrix(alchem.overlap_matrix))
        else:
            return (data, None)

    @staticmethod
    def _analyse_amber(work_dir=None, estimator="MBAR", method="alchemlyb", **kwargs):
        """
        Analyse the AMBER free energy data.

           Parameters
           ----------

           work_dir : str
               The path to the working directory.

        estimator : str
            The estimator ('MBAR' or 'TI') used.

        method : str
            The method used to analyse the data. Options are "alchemlyb" or "native".

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap or dHdl : numpy.matrix or alchemlyb.estimators.ti_.TI
            For MBAR, this returns the overlap matrix for the overlap between
            each lambda window. For TI, this returns None.
        """

        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        if estimator.replace(" ", "").upper() not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        if not isinstance(method, str):
            raise TypeError("'method' must be of type 'str'.")
        if method.replace(" ", "").upper() != "ALCHEMLYB":
            raise ValueError(
                "Only 'alchemyb' is supported as an analysis method for AMBER."
            )

        # Find the output files and work out the lambda windows from the directory names.
        glob_path = _pathlib.Path(work_dir)
        files = sorted(glob_path.glob("**/*.out"))
        lambdas = []
        for file in files:
            for part in file.parts:
                if "lambda" in part:
                    lambdas.append(float(part.split("_")[-1]))

        # Find the temperature for each lambda window.
        temperatures = []
        for file, lambda_ in zip(files, lambdas):
            found_temperature = False
            with open(file) as f:
                for line in f.readlines():
                    if not found_temperature:
                        match = _re.search(r"temp0=([\d.]+)", line)
                        if match is not None:
                            temperatures += [float(match.group(1))]
                            found_temperature = True
                            break

                if not found_temperature:
                    raise ValueError(
                        "The temperature was not detected in the AMBER output file."
                    )

        if temperatures[0] != temperatures[-1]:
            raise ValueError("The temperatures at the endstates don't match!")

        return Relative._analyse_internal(
            files, temperatures, lambdas, "AMBER", estimator, **kwargs
        )

    @staticmethod
    def _analyse_gromacs(work_dir=None, estimator="MBAR", method="alchemlyb", **kwargs):
        """
        Analyse GROMACS free energy data.

        Parameters
        ----------

        work_dir : str
            The path to the working directory.

        estimator : str
            The estimator ('MBAR' or 'TI') used.

        method : str
            The method used to analyse the data. Options are "alchemlyb" or "native".

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap or dHdl : numpy.matrix or alchemlyb.estimators.ti_.TI
            For MBAR, this returns the overlap matrix for the overlap between
            each lambda window. For TI, this returns None.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        if not estimator.replace(" ", "").upper() in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        if not isinstance(method, str):
            raise TypeError("'method' must be of type 'str'.")
        method = method.replace(" ", "").upper()
        if method not in ["ALCHEMLYB", "NATIVE"]:
            raise ValueError("'method' must be either 'alchemlyb' or 'native'.")

        if method == "ALCHEMLYB":
            # Find the output files and work out the lambda windows from the directory names.
            glob_path = _pathlib.Path(work_dir)
            files = sorted(glob_path.glob("**/[!bar]*.xvg"))
            lambdas = []
            for file in files:
                for part in file.parts:
                    if "lambda" in part:
                        lambdas.append(float(part.split("_")[-1]))

            # Find the temperature at each lambda window
            temperatures = []
            for file in files:
                found_temperature = False
                with open(file, "r") as f:
                    for line in f.readlines():
                        t = None
                        start = "T ="
                        end = "(K)"
                        if start and end in line:
                            t = int(((line.split(start)[1]).split(end)[0]).strip())
                            temperatures.append(float(t))
                            if t is not None:
                                found_temperature = True
                                break

                if not found_temperature:
                    raise ValueError(
                        f"The temperature was not detected in the GROMACS output file, {file}"
                    )

            if temperatures[0] != temperatures[-1]:
                raise ValueError("The temperatures at the endstates don't match!")

            return Relative._analyse_internal(
                files, temperatures, lambdas, "GROMACS", estimator, **kwargs
            )

        # For the older gromacs versions and native use the gmx bar analysis.
        elif method == "NATIVE":
            if _gmx_exe is None:
                raise _MissingSoftwareError(
                    "Cannot use native gmx bar analysis as GROMACS is not installed!"
                )

            _warnings.warn(
                "using 'native' for GROMACS does not return an overlap/dHdl."
            )
            _warnings.warn("using 'native' for GROMACS uses BAR.")

            # Create the command.
            glob_path = _pathlib.Path(work_dir)
            xvg_files = sorted(glob_path.glob("**/[!bar]*.xvg"))
            xvg_files = [str(file.absolute()) for file in xvg_files]
            command = "%s bar -f %s -o %s/bar.xvg" % (
                _gmx_exe,
                " ".join(xvg_files),
                work_dir,
            )

            # Run the command.
            proc = _subprocess.run(
                _Utils.command_split(command),
                shell=False,
                stdout=_subprocess.PIPE,
                stderr=_subprocess.PIPE,
            )
            if proc.returncode != 0:
                raise _AnalysisError("native GROMACS free-energy analysis failed!")

            # Initialise list to hold the data.
            data = []

            # Extract the data from the output files.
            with open("%s/bar.xvg" % work_dir) as file:
                # Read all of the lines into a list.
                lines = []
                for line in file:
                    # Ignore comments and xmgrace directives.
                    if line[0] != "#" and line[0] != "@":
                        lines.append(line.rstrip())

                # Store the initial free energy reading.
                data.append(
                    (
                        0.0,
                        0.0 * _Units.Energy.kcal_per_mol,
                        0.0 * _Units.Energy.kcal_per_mol,
                    )
                )

                # Zero the accumulated error.
                total_error = 0

                # Zero the accumulated free energy difference.
                total_freenrg = 0

                # Process the BAR data.
                for x, line in enumerate(lines):
                    # Extract the data from the line.
                    records = line.split()

                    # Update the total free energy difference.
                    total_freenrg += float(records[1])

                    # Extract the error.
                    error = float(records[2])

                    # Update the accumulated error.
                    total_error = _math.sqrt(total_error * total_error + error * error)

                    # Append the data.
                    data.append(
                        (
                            (x + 1) / (len(lines)),
                            (total_freenrg * _Units.Energy.kt).kcal_per_mol(),
                            (total_error * _Units.Energy.kt).kcal_per_mol(),
                        )
                    )

            return (data, None)

    @staticmethod
    def _analyse_somd(work_dir=None, estimator="MBAR", method="alchemlyb", **kwargs):
        """
        Analyse SOMD free energy data.

        Parameters
        ----------

        work_dir : str
            The path to the working directory.

        estimator : str
            The estimator ('MBAR' or 'TI') used.

        method : str
            The method used to analyse the data. Options are "alchemlyb" or "native".

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap or dHdl : numpy.matrix or alchemlyb.estimators.ti_.TI
            For MBAR, this returns the overlap matrix for the overlap between
            each lambda window. For TI, this returns None.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        estimator = estimator.replace(" ", "").upper()
        if estimator not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        if not isinstance(method, str):
            raise TypeError("'method' must be of type 'str'.")
        method = method.replace(" ", "").upper()
        if not method in ["ALCHEMLYB", "NATIVE"]:
            raise ValueError("'method' must be either 'alchemlyb' or 'native'.")

        if method == "ALCHEMLYB":
            # Glob the data files and work out the lambda values.
            glob_path = _pathlib.Path(work_dir)
            files = sorted(glob_path.glob("**/simfile.dat"))
            lambdas = []
            for file in files:
                for part in file.parts:
                    if "lambda" in part:
                        lambdas.append(float(part.split("_")[-1]))

            temperatures = []
            for file in files:
                found_temperature = False
                with open(file, "r") as f:
                    for line in f.readlines():
                        temp = None
                        start = "#Generating temperature is"
                        if start in line:
                            split_line = (line.split(start)[1]).strip().split(" ")
                            temp = split_line[0]
                            unit = split_line[-1]
                            if unit.upper() == "C":
                                temp = float(temp) + 273.15
                            else:
                                temp = float(temp)
                            temperatures.append(temp)
                            if temp is not None:
                                found_temperature = True
                                break

                if not found_temperature:
                    raise ValueError(
                        f"The temperature was not detected in the SOMD output file: {file}"
                    )

            if temperatures[0] != temperatures[-1]:
                raise ValueError("The temperatures at the endstates don't match!")

            return Relative._analyse_internal(
                files, temperatures, lambdas, "SOMD", estimator, **kwargs
            )

        elif method == "NATIVE":
            # Create the command.
            command = (
                "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar.txt --overlap --percent 5"
                % (_analyse_freenrg, work_dir, work_dir)
            )

            # Run the first command.
            proc = _subprocess.run(
                _Utils.command_split(command),
                shell=False,
                stdout=_subprocess.PIPE,
                stderr=_subprocess.PIPE,
            )
            if proc.returncode != 0:
                raise _AnalysisError("Native SOMD free-energy analysis failed!")

            # Re-run without subsampling if the subsampling has resulted in less than 50 samples.
            with open("%s/mbar.txt" % work_dir) as file:
                for line in file:
                    if (
                        "#WARNING SUBSAMPLING ENERGIES RESULTED IN LESS THAN 50 SAMPLES"
                        in line
                    ):
                        _warnings.warn(
                            "Subsampling resulted in less than 50 samples, "
                            f"re-running without subsampling for '{work_dir}'"
                        )
                        command = (
                            "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar.txt --overlap"
                            % (_analyse_freenrg, work_dir, work_dir)
                        )
                        proc = _subprocess.run(
                            _Util.command_split(command),
                            shell=False,
                            stdout=_subprocess.PIPE,
                            stderr=_subprocess.PIPE,
                        )
                        if proc.returncode != 0:
                            raise _AnalysisError("SOMD free-energy analysis failed!")
                        break

            # Initialise list to hold the data.
            data = []

            # Initialise list to hold the overlap matrix.
            overlap = []

            # Extract the data from the output files.
            with open("%s/mbar.txt" % work_dir) as file:
                # Process the MBAR data.
                for line in file:
                    # Process the overlap matrix.
                    if "#Overlap matrix" in line:
                        # Get the next row.
                        row = next(file)

                        # Loop until we hit the next section.
                        while not row.startswith("#DG"):
                            # Extract the data for this row.
                            records = [float(x) for x in row.split()]

                            # Append to the overlap matrix.
                            overlap.append(records)

                            # Get the next line.
                            row = next(file)

                    # Process the PMF.
                    elif "PMF from MBAR" in line:
                        # Get the next row.
                        row = next(file)

                        # Loop until we hit the next section.
                        while not row.startswith("#TI"):
                            # Split the line.
                            records = row.split()

                            # Append the data.
                            data.append(
                                (
                                    float(records[0]),
                                    float(records[1]) * _Units.Energy.kcal_per_mol,
                                    float(records[2]) * _Units.Energy.kcal_per_mol,
                                )
                            )

                            # Get the next line.
                            row = next(file)

        return (data, _np.matrix(overlap))

    @staticmethod
    def _analyse_somd2(work_dir=None, estimator="MBAR", method="alchemlyb", **kwargs):
        """
        Analyse SOMD2 free energy data.

        Parameters
        ----------

        work_dir : str
            The path to the working directory.

        estimator : str
            The estimator ('MBAR' or 'TI') used.

        method : str
            The method used to analyse the data. Options are "alchemlyb" or "native".

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap or dHdl : numpy.matrix or alchemlyb.estimators.ti_.TI
            For MBAR, this returns the overlap matrix for the overlap between
            each lambda window. For TI, this returns None.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(estimator, str):
            raise TypeError("'estimator' must be of type 'str'.")
        estimator = estimator.replace(" ", "").upper()
        if estimator not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        if not isinstance(method, str):
            raise TypeError("'method' must be of type 'str'.")
        method = method.replace(" ", "").upper()
        if method != "ALCHEMLYB":
            raise ValueError(
                "Only 'alchemlyb' analysis 'method' is supported for SOMD2."
            )

        # Glob the data files.
        glob_path = _pathlib.Path(work_dir)
        files = sorted(glob_path.glob("**/*.parquet"))

        # Loop over each file and try to extract the metadata to work out
        # the lambda value and temperature for each window.

        lambdas = []
        temperatures = []

        for file in files:
            try:
                metadata = _json.loads(
                    _pq.read_metadata(file).metadata["somd2".encode()]
                )
                lambdas.append(float(metadata["lambda"]))
                temperatures.append(float(metadata["temperature"]))
            except:
                raise IOError(f"Unable to parse metadata from SOMD2 file: {file}")

        # Sort the lists based on the lambda values.
        temperatures = [x for _, x in sorted(zip(lambdas, temperatures))]
        lambdas = sorted(lambdas)

        # Check that the temperatures at the end states match.
        if temperatures[0] != temperatures[-1]:
            raise ValueError("The temperatures at the endstates don't match!")

        return Relative._analyse_internal(
            files, temperatures, lambdas, "SOMD2", estimator, **kwargs
        )

    def _checkOverlap(self, estimator="MBAR", threshold=0.03):
        """
        Check the overlap of an FEP simulation leg.

        Parameters
        ----------

        estimator : str
            The estimator used for the free-energy analysis. ("MBAR" or "TI")

        threshold : float
            The threshold value used to check the off-diagonals.

        Returns
        -------

        is_okay : boolean
             True if the overlap is okay, False if any off-diagonals are less
             than the threshold value.

        num_low : int
            The number of off-diagonals that are less than the threshold value.
        """

        # Calculate the overlap for this object.
        _, overlap = self.analyse(estimator=estimator)

        if overlap:
            return Relative.checkOverlap(overlap, threshold=threshold)
        else:
            raise ValueError("Overlap matrix isn't supported for this estimator.")

    def _difference(self, pmf_ref):
        """
        Compute the relative free-energy difference between two perturbation
        legs.

        Parameters
        ----------

        pmf_ref : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The reference potential of mean force (PMF). The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The relative free-energy difference and its associated error.
        """

        # Calculate the PMF for this object.
        pmf, _ = self.analyse()

        # Now call the staticmethod passing in both PMFs.
        return Relative.difference(pmf, pmf_ref)

    @staticmethod
    def checkOverlap(overlap, estimator="MBAR", threshold=0.03):
        """Check the overlap of an FEP leg. 

           Parameters
           ----------

           overlap : [ [ float, float, ... ] ], numpy.ndarray
               The overlap matrix. This gives the overlap between lambda windows.

           estimator : str
               Must be "MBAR" for checking the overlap matrix.

           threshold : float
               The threshold value used to check the off-diagonals. Default is 0.03 .

           Returns
           -------

           overlap_okay : boolean
                True if the overlap is okay, False if any off-diagonals are less than the threshold value.

        """
        if not isinstance(overlap, _np.ndarray):
            raise TypeError("'overlap' must be of type 'numpy.matrix'.")

        # estimator must be MBAR for overlap matrix or TI for dhdl plot.
        if estimator not in ['MBAR']:
            raise ValueError("'estimator' must be 'MBAR'.")
        
        if not isinstance(threshold, float):
            raise TypeError("'threshold' must be of type 'float'.")

        if estimator == "MBAR":
            # check the overlap
            # get all off diagonals
            off_diagonal = (_np.diagonal(overlap, 1)).tolist()
            for a in (_np.diagonal(overlap, -1)).tolist():
                off_diagonal.append(a)

            # check if the off diagonals are less than the threshold value or larger.
            too_small = 0
            overlap_okay = False
            for o in off_diagonal:
                if o < threshold:
                    too_small += 1
            if too_small > 0:
                _warnings.warn(f"Overlap matrix is bad - {too_small} off-diagonals are less than {threshold}.")
            else:
                overlap_okay = True

        return overlap_okay, too_small

    def _check_overlap(self):
        """Check the overlap of an FEP leg. 

           Parameters
           ----------

           Returns
           -------

           overlap_okay : boolean
                True if the overlap is okay, False if any off-diagonals are less than 0.03.

        """

        # Calculate the overlap for this object.
        _, overlap = self.analyse()

        # Now call the staticmethod passing in the overlap and the work_dir of the run.
        return Relative.checkOverlap(overlap, estimator=self._estimator)

    @staticmethod
    def plot(overlap_dhdl, estimator=None, work_dir=None, file_name=None):
        """Plot either the overlap or the dhdl of the transformation. 

           Parameters
           ----------

           overlap_dhdl : numpy.ndarray or alchemlyb.estimators.ti_.TI
               For MBAR, this is the overlap matrix for the overlap between each lambda window.
               For TI, this the dHdl gradients from the alchemlyb analysis.

           estimator : str
               The estimator ('MBAR' or 'TI') used.

           work_dir : str
               The working directory for the free-energy perturbation simulation.
               If this is specified, the plot will be saved there.

           file_name : str
               The name for the saved file.
               If this is None, the file name will be either "overlap_MBAR.png" or "dHdl_TI.png".

           Returns
           -------
           the plot : matplotlib.axes._subplots.AxesSubplot

        """

        if work_dir:
            if not isinstance(work_dir, str):
                raise TypeError("'work_dir' must be of type 'str'.")
            if not _os.path.isdir(work_dir):
                raise ValueError("'work_dir' doesn't exist!")
        else:
            pass

        # estimator must be MBAR for overlap matrix or TI for dhdl plot.
        if estimator not in ['MBAR', 'TI', None]:
            raise ValueError("'estimator' must be 'MBAR' or 'TI'. If 'None, data type will be inferred.")

        if estimator is None:
            if isinstance(overlap_dhdl, _np.ndarray):
                estimator = "MBAR"
            elif isinstance(overlap_dhdl, _alchemlyb.estimators.ti_.TI):
                estimator = "TI"
            else:
                raise TypeError("Data type for estimator = 'None' could not be inferred / does not match allowed data types.")

        if file_name:
            if not isinstance(file_name, str):
                raise TypeError("'file_name' must be of type 'str'.")
        else:
            if estimator == "MBAR":
                file_name = "overlap_MBAR"
            elif estimator == "TI":
                file_name = "dHdl_TI"

        if estimator == "MBAR":
            if not isinstance(overlap_dhdl, _np.ndarray):
                        raise TypeError("'overlap' must be of type 'numpy.ndarray' for 'MBAR'.\
                            This is obtained from running analysis using estimator='MBAR'.")

            # use the alchemlyb functionality to plot the overlap matrix
            ax = _plot_mbar_overlap_matrix(overlap_dhdl)
            ax.set_title(f"overlap matrix")

            if work_dir is not None:
                ax.figure.savefig(
                    f"{work_dir}/{file_name}.png", bbox_inches='tight', pad_inches=0.0)      

        elif estimator == 'TI':
            if not isinstance(overlap_dhdl, _alchemlyb.estimators.ti_.TI):
                raise TypeError("'overlap' must be of type 'alchemlyb.estimators.ti_.TI' for 'TI'.\
                                    This is obtained from running analysis using estimator='TI'.")

            # use the alchemlyb functionality to plot the dhdl
            ax = _plot_ti_dhdl(overlap_dhdl)
            ax.set_title(f"dhdl plot")
            
            if work_dir is not None:
                ax.figure.savefig(
                    f"{work_dir}/{file_name}.png")
            
        return(ax)

    def _plot(self):
        """Plot either the overlap or the dhdl of the transformation.
           Saves the plot in the working directory of the run.

           Parameters
           ----------

           Returns
           -------
           the plot : matplotlib.axes._subplots.AxesSubplot

        """

        # Calculate the overlap for this object.
        _, overlap = self.analyse()

        # Now call the staticmethod passing in the overlap, estimator and the work_dir of the run.
        return Relative.plot(overlap, estimator=self._estimator, work_dir=self._work_dir)

    def _initialise_runner(self, system):
        """
        Internal helper function to initialise the process runner.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.
        """

        # Initialise list to store the processe
        processes = []

        # Convert to an appropriate water topology.
        if self._engine == "SOMD" or self._engine == "AMBER":
            system._set_water_topology("AMBER", property_map=self._property_map)
        elif self._engine == "GROMACS":
            system._set_water_topology("GROMACS", property_map=self._property_map)

        # Setup all of the simulation processes for each leg.

        # Get the lambda values from the protocol for the first leg.
        lam_vals = self._protocol.getLambdaValues()

        # Create a process for the first lambda value.
        lam = lam_vals[0]

        # Update the protocol lambda values.
        self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

        # Create and append the required processes for each leg.
        # Nest the working directories inside self._work_dir.

        # Name the first directory.
        first_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)

        # SOMD.
        if self._engine == "SOMD":
            # Check for GPU support.
            if "CUDA_VISIBLE_DEVICES" in _os.environ:
                platform = "CUDA"
            else:
                platform = "CPU"

            first_process = _Process.Somd(system,
                                          self._protocol,
                                          platform=platform,
                                          work_dir=first_dir,
                                          extra_options=self._extra_options,
                                          extra_lines=self._extra_lines,
                                          property_map=self._property_map,
                                          )

        # GROMACS.
        elif self._engine == "GROMACS":
            first_process = _Process.Gromacs(system,
                                          self._protocol,
                                          work_dir=first_dir,
                                          ignore_warnings=self._ignore_warnings,
                                          show_errors=self._show_errors,
                                          extra_options=self._extra_options,
                                          extra_lines=self._extra_lines,
                                          property_map=self._property_map,
                                          )

        # AMBER.
        elif self._engine == "AMBER":
            first_process = _Process.Amber(system,
                                           self._protocol,
                                           exe=self._exe,
                                           work_dir=first_dir,
                                           extra_options=self._extra_options,
                                           extra_lines=self._extra_lines,
                                           property_map=self._property_map,
                                           )

        if self._setup_only:
            del(first_process)
        else:
            processes.append(first_process)

        # Loop over the rest of the lambda values.
        for x, lam in enumerate(lam_vals[1:]):
            # Name the directory.
            new_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)

            # Use absolute path.
            if not _os.path.isabs(new_dir):
                new_dir = _os.path.abspath(new_dir)

            # Delete any existing directories.
            if _os.path.isdir(new_dir):
                _shutil.rmtree(new_dir, ignore_errors=True)

            # Copy the first directory to that of the current lambda value.
            _shutil.copytree(first_dir, new_dir)

            # Update the protocol lambda values.
            self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

            # Now update the lambda values in the config files.

            # SOMD.
            if self._engine == "SOMD":
                new_config = []
                with open(new_dir + "/somd.cfg", "r") as f:
                    for line in f:
                        if "lambda_val" in line:
                            new_config.append("lambda_val = %s\n" % lam)
                        else:
                            new_config.append(line)
                with open(new_dir + "/somd.cfg", "w") as f:
                    for line in new_config:
                        f.write(line)

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/somd.out"
                    process._std_err_file = new_dir + "/somd.err"
                    process._rst_file = new_dir + "/somd.rst7"
                    process._top_file = new_dir + "/somd.prm7"
                    process._traj_file = new_dir + "/traj000000001.dcd"
                    process._restart_file = new_dir + "/latest.rst"
                    process._config_file = new_dir + "/somd.cfg"
                    process._pert_file = new_dir + "/somd.pert"
                    process._gradients_file = new_dir + "/gradients.dat"
                    process._input_files = [
                        process._config_file,
                        process._rst_file,
                        process._top_file,
                        process._pert_file,
                    ]
                    processes.append(process)

            # GROMACS.
            elif self._engine == "GROMACS":
                # Delete the existing binary run file.
                _os.remove(new_dir + "/gromacs.tpr")
                new_config = []
                with open(new_dir + "/gromacs.mdp", "r") as f:
                    for line in f:
                        if "init-lambda-state" in line:
                            new_config.append("init-lambda-state = %d\n" % (x + 1))
                        else:
                            new_config.append(line)
                with open(new_dir + "/gromacs.mdp", "w") as f:
                    for line in new_config:
                        f.write(line)

                mdp = new_dir + "/gromacs.mdp"
                mdp_out = new_dir + "/gromacs.out.mdp"
                gro = new_dir + "/gromacs.gro"
                top = new_dir + "/gromacs.top"
                tpr = new_dir + "/gromacs.tpr"

                # Use grompp to generate the portable binary run input file.
                command = "%s grompp -f %s -po %s -c %s -p %s -r %s -o %s" % (
                    _gmx_exe,
                    mdp,
                    mdp_out,
                    gro,
                    top,
                    gro,
                    tpr,
                )

                # Run the command. If this worked for the first lambda value,
                # then it should work for all others.
                proc = _subprocess.run(
                    _Utils.command_split(command),
                    shell=False,
                    text=True,
                    stdout=_subprocess.PIPE,
                    stderr=_subprocess.PIPE,
                )

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/gromacs.out"
                    process._std_err_file = new_dir + "/gromacs.err"
                    process._gro_file = new_dir + "/gromacs.gro"
                    process._top_file = new_dir + "/gromacs.top"
                    process._traj_file = new_dir + "/gromacs.trr"
                    process._config_file = new_dir + "/gromacs.mdp"
                    process._tpr_file = new_dir + "/gromacs.tpr"
                    process._input_files = [
                        process._config_file,
                        process._gro_file,
                        process._top_file,
                        process._tpr_file,
                    ]
                    processes.append(process)

            # AMBER.
            elif self._engine == "AMBER":
                new_config = []
                with open(new_dir + "/amber.cfg", "r") as f:
                    for line in f:
                        if "clambda" in line:
                            new_config.append("   clambda=%s,\n" % lam)
                        else:
                            new_config.append(line)
                with open(new_dir + "/amber.cfg", "w") as f:
                    for line in new_config:
                        f.write(line)

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/amber.out"
                    process._std_err_file = new_dir + "/amber.err"
                    process._rst_file = new_dir + "/amber.rst7"
                    process._top_file = new_dir + "/amber.prm7"
                    process._traj_file = new_dir + "/amber.nc"
                    process._config_file = new_dir + "/amber.cfg"
                    process._nrg_file = new_dir + "/amber.nrg"
                    process._input_files = [process._config_file,
                                            process._rst_file,
                                            process._top_file]
                    processes.append(process)

        if not self._setup_only:
            # Initialise the process runner. All processes have already been nested
            # inside the working directory so no need to re-nest.
            self._runner = _Process.ProcessRunner(processes)

    def _update_run_args(self, args):
        """
        Internal function to update run arguments for all subprocesses.

        Parameters
        ----------

        args : dict, collections.OrderedDict
            A dictionary which contains the new command-line arguments
            for the process executable.
        """

        if not isinstance(args, dict):
            raise TypeError("'args' must be of type 'dict'")

        for process in self._runner.processes():
            process.setArgs(args)


def getData(name="data", file_link=False, work_dir=None):
    """
    Return a link to a zip file containing the data files required for
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
    return Relative.getData(name=name, file_link=file_link, work_dir=work_dir)

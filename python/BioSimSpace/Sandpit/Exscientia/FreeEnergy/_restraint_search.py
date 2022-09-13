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

__author__ = "Finlay Clark"
__email__ = "finlay.clark@ed.ac.uk"

__all__ = ["RestraintSearch"]

import sys as _sys
import os as _os
import tempfile as _tempfile
import warnings as _warnings

# Import tqdm for progress bars for BSS restraints derivation. 
# The import needs to be different depending on
# whether we are running in a notebook or not

try:
    shell = get_ipython().__class__.__name__
    if shell == 'ZMQInteractiveShell': # This is a notebook
        from tqdm.notebook import tqdm 
    else:
        from tqdm import tqdm
except NameError: # Import progress bars to work in terminal
    from tqdm import tqdm

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

try:
    import MDRestraintsGenerator
    from MDRestraintsGenerator import search as _search
    from MDRestraintsGenerator.restraints import FindBoreschRestraint as _FindBoreschRestraint
    is_MDRestraintsGenerator = True
except:
    print('Please install MDRestraintsGenerator for analysis using it.')
    is_MDRestraintsGenerator = False

import numpy as _np
import MDAnalysis as _mda 
from MDAnalysis.analysis.distances import dist as _dist
from MDAnalysis.lib.distances import calc_angles as _calc_angles
from MDAnalysis.lib.distances import calc_dihedrals as _calc_dihedrals
from scipy.stats import circmean as _circmean
import matplotlib.pyplot as _plt

from Sire.Base import getBinDir as _getBinDir
from Sire.Base import getShareDir as _getShareDir
from Sire.Units import k_boltz # kcal / (mol K)

from .. import _gmx_exe
from .. import _is_notebook
from .._Exceptions import AnalysisError as _AnalysisError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from ..Trajectory._trajectory import Trajectory as _Trajectory
from .. import Process as _Process
from .. import Protocol as _Protocol
from ..Types import Temperature as _Temperature
from .. import Units as _Units
from ..Units.Length import angstrom as _angstrom
from ..Units.Angle import radian as _radian
from ..Units.Angle import degree as _degree
from ..Units.Energy import kcal_per_mol as _kcal_per_mol

from ..MD._md import _find_md_engines

if _is_notebook:
    from IPython.display import FileLink as _FileLink

from ._restraint import Restraint

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
    _engines = ["GROMACS", "SOMD"] 

    def __init__(self, system, protocol=None, work_dir=None, engine=None,
            gpu_support=False, ignore_warnings=False,
            show_errors=True, extra_options=None, extra_lines=None,
            property_map={}, **kwargs):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the ABFE simulation. This must contain
               a single decoupled molecule and is assumed to have already
               been equilibrated.

           protocol : :class:`Protocol.Production <BioSimSpace.Protocol.Production>`, \
               The simulation production protocol.

           work_dir : str
               The working directory for the ABFE restraint generation
                simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "AMBER", "GROMACS", or "SOMD". If this argument is omitted
               then BioSimSpace will choose an appropriate engine for you.

           gpu_support : bool
               Whether the engine must have GPU support.

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

           kwargs :
               Keyword arguments to be passed to the BSS.Process.
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
            self._protocol = _Protocol.Production()

        self._extra_options = extra_options if extra_options is not None else {}
        self._extra_lines = extra_lines if extra_lines is not None else []

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
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

        # Initialise the process.
        self._initialise_process(self._system, gpu_support, **kwargs)

    def start(self):
        """Start the simulation."""
        self._process.start()

    def wait(self):
        """Wait for the simulation to finish."""
        self._process.wait()

    def kill(self):
        """Kill the process."""
        self._process.kill()

    def isRunning(self):
        """Check if the process is running."""
        return self._process.isRunning()

    def workDir(self):
        """Return the working directory.

           Returns
           -------

           work_dir : str
               The path of the working directory.
        """
        return self._work_dir

    def _analyse(self, rest_type='Boresch',
                method='MDRestraintsGenerator',
                append_to_lig_selection="",
                recept_selection_str='protein and name CA C N',
                cutoff=5, # In Angstrom
                restraint_idx=0,
                block='AUTO'):
        """Analyse trajectory and select restraints which best mimic strongest
           receptor-ligand interactions.

            Parameters
            ----------

            rest_type: str
                The type of restraints to select (currently only Boresch is available).
                Default is 'Boresch'.

            method: str
                The method to use to derive the restraints. Currently only 'MDRestraintsGenerator'
                is supported.

            append_to_lig_selection: str
                Appends the supplied string to the default atom selection which chooses
                the atoms in the ligand to consider as potential anchor points. The default
                atom selection is f'resname {ligand_resname} and not name H*'. Uses the
                mdanalysis atom selection language. For example, 'not name O*' will result
                in an atom selection of f'resname {ligand_resname} and not name H* and not
                name O*'.

            recept_selection_str: str
                The selection string for the atoms in the receptor to consider
                as potential anchor points. The default atom selection is
                'protein and name CA C N'. Uses the mdanalysis atom selection
                language.

            cutoff: float
                The greatest distance between ligand and receptor anchor atoms, in
                Angstrom. Receptor anchors further than cutoff Angstroms from the closest
                ligand anchors will not be included in the search for potential anchor points.

            restraint_idx: int
                The index of the restraint from a list of candidate restraints ordered by
                suitability. restraint_idx != 0 is only valid if method == 'BSS'.

            block : bool
                Whether to block until the process has finished running.

            Returns
            -------

            restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
                The restraints of `rest_type` which best mimic the strongest receptor-ligand
                interactions.
            """
        # Wait for the process to finish.
        if block is True or block == "AUTO":
            self.wait()

        # Return the result of calling the staticmethod, passing in the working
        # directory of this object.
        return RestraintSearch.analyse(self._work_dir, self._system,
                self._process.getTrajectory(),
                self._protocol.getTemperature(),
                rest_type=rest_type,
                method=method,
                append_to_lig_selection=append_to_lig_selection,
                recept_selection_str=recept_selection_str,
                cutoff=cutoff,
                restraint_idx=restraint_idx)

    def _initialise_process(self, system, gpu_support, **kwargs):
        """Internal helper function to initialise the process.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           gpu_support : bool
               Whether the engine must have GPU support.

           kwargs :
               Keyword arguments to be passed to the BSS.Process.

        """

        # Convert to an appropriate AMBER topology. (Required by SOMD for its
        # FEP setup.)
        if self._engine == "SOMD":
            system._set_water_topology("AMBER", property_map=self._property_map)

        # Setup all of the simulation processes.

        # SOMD.
        if self._engine == "SOMD":
            # Check for GPU support.
            if "CUDA_VISIBLE_DEVICES" in _os.environ:
                platform = "CUDA"
            else:
                if gpu_support:
                    raise ValueError("gpu_support cannot be True if CUDA_VISIBLE_DEVICES is not set.")
                platform = "CPU"

            self._process = _Process.Somd(system, self._protocol,
                platform=platform, work_dir=self._work_dir,
                property_map=self._property_map, extra_options=self._extra_options,
                extra_lines=self._extra_lines, **kwargs)

        # GROMACS.
        elif self._engine == "GROMACS":
            self._process = _Process.Gromacs(system, self._protocol,
                work_dir=self._work_dir, ignore_warnings=self._ignore_warnings,
                show_errors=self._show_errors, extra_options=self._extra_options,
                extra_lines=self._extra_lines, **kwargs)
            if gpu_support:
                self._process.setArg("-update", "gpu")

        # AMBER.
        elif self._engine == "AMBER":
            self._process = _Process.Amber(system, self._protocol, exe=self._exe,
                work_dir=self._work_dir, extra_options=self._extra_options,
                extra_lines=self._extra_lines, **kwargs)


    @staticmethod
    def analyse(work_dir, system, traj, temperature, rest_type='Boresch',
                method='MDRestraintsGenerator',
                append_to_lig_selection="",
                recept_selection_str='protein and name CA C N',
                cutoff=5,
                restraint_idx=0): # In Angstrom
        """Analyse existing trajectory from a simulation working directory and
        select restraints which best mimic the strongest receptor-ligand
        interactions.

           Parameters
           ----------

           work_dir : str
               The working directory for the simulation.

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the ABFE restraint calculation. This
               must contain a single decoupled molecule and is assumed to have
               already been equilibrated.

           traj : :class:`System <BioSimSpace.Trajectory._trajectory.Trajectory>`
               The trajectory for analysis.

           temperature : :class:`System <BioSimSpace.Types.Temperature>`
               The temperature of the system

           rest_type: str
               The type of restraints to select (currently only Boresch is available).
               Default is ``Boresch``.
           
           method: str
                The method to use to derive the restraints. 'BSS' or 'MDRestraintsGenerator'. 
                BSS uses the native BioSimSpace derivation.

           append_to_lig_selection: str
               Appends the supplied string to the default atom selection which chooses
               the atoms in the ligand to consider as potential anchor points. The default
               atom selection is f'resname {ligand_resname} and not name H*'. Uses the
               mdanalysis atom selection language. For example, 'not name O*' will result
               in an atom selection of f'resname {ligand_resname} and not name H* and not 
               name O*'.

           recept_selection_str: str
               The selection string for the atoms in the receptor to consider
               as potential anchor points. The default atom selection is
               'protein and name CA C N'. Uses the mdanalysis atom selection
               language.

           cutoff: float
               The greatest distance between ligand and receptor anchor atoms, in
               Angstrom. Receptor anchors further than cutoff Angstroms from the closest
               ligand anchors will not be included in the search for potential anchor points.

            restraint_idx: int
                The index of the restraint from a list of candidate restraints ordered by
                suitability. restraint_idx != 0 is only valid if method == 'BSS'.

           Returns
           -------

           restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
               The restraints of `rest_type` which best mimic the strongest receptor-ligand
               interactions.

        """
        # Check all inputs

        if not isinstance(work_dir, str):
            raise TypeError(f"work_dir: {work_dir} must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError(f"work_dir: {work_dir} doesn't exist!")

        if not isinstance(system, _System):
            raise TypeError(f"system {type(system)} must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of solvated system.
            _system = system.copy()

        if not isinstance(traj, _Trajectory):
            raise TypeError(f"traj {type(traj)} must be of type 'BioSimSpace.Trajectory._trajectory.Trajectory'")

        if not isinstance(temperature, _Temperature):
            raise ValueError(
                f"temperature {type(temperature)} must be of type 'BioSimSpace.Types.Temperature'")

        if not isinstance(rest_type, str):
            raise TypeError(f"rest_type {type(rest_type)} must be of type 'str'.")
        if not rest_type.lower() == 'boresch':
            raise NotImplementedError("Only Boresch restraints are currently implemented")
        
        if not isinstance(method, str):
            raise TypeError(f"method {type(method)} must be of type 'str'.")
        if not method.lower() in ['mdrestraintsgenerator', 'bss']:
            raise NotImplementedError("Deriving restraints using 'MDRestraintsGenerator'"
                                      "or 'BSS' are the only options implemented.")
                            
        if not isinstance(append_to_lig_selection, str):
            raise TypeError(f"append_to_lig_selection {type(append_to_lig_selection)} must be of type 'str'.")

        if not isinstance(recept_selection_str, str):
            raise TypeError(f"append_to_recept_selection {type(recept_selection_str)} must be of type 'str'.")
        
        if not isinstance(cutoff, (int, float)):
            raise TypeError(f"cutoff {type(cutoff)} must be of type 'int' or 'float'.")
        
        # There must be a single molecule to be decoupled (or annihilated).
        if system.nDecoupledMolecules() != 1:
            raise ValueError("The system must contain a single molecule to be decoupled! "
                                "Use the 'BioSimSpace.Align.Decouple' function to mark a molecule"
                                " to be decoupled.")

        # Extract restraints from simulation
        # Get mdanalysis universe object
        u = traj.getTrajectory(format='mdanalysis')

        # Find decoupled molecule and use it to create ligand selection
        decoupled_mol = _system.getDecoupledMolecules()[0]
        decoupled_resname = decoupled_mol.getResidues()[0].name()

        lig_selection_str = f'((resname {decoupled_resname}) and (not name H*))'
        if append_to_lig_selection:
            lig_selection_str += ' and '
            lig_selection_str += append_to_lig_selection

        if rest_type.lower() == 'boresch':
            return RestraintSearch._boresch_restraint(
                u, system, temperature, lig_selection_str,
                recept_selection_str, method, work_dir, cutoff,
                restraint_idx=restraint_idx)

    @staticmethod
    def _boresch_restraint(u, system, temperature, lig_selection_str,
                           recept_selection_str, method, work_dir, cutoff,
                           restraint_idx=0):
        """Generate the Boresch Restraint.

           Parameters
           ----------

           u : MDAnalysis.Universe
               The trajectory for the ABFE restraint calculation as a
               MDAnalysis.Universe object.

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the ABFE restraint calculation. This
               must contain a single decoupled molecule and is assumed to have
               already been equilibrated.

           temperature : :class:`System <BioSimSpace.Types.Temperature>`
               The temperature of the system

           lig_selection_str: str
               The selection string for the atoms in the ligand to consider
               as potential anchor points.

           recept_selection_str: str
               The selection string for the atoms in the receptor to consider
               as potential anchor points. Uses the mdanalysis atom selection
               language.

           method: str
               The method to use to derive the restraints. 'BSS' or 'MDRestraintsGenerator'.
               BSS uses the native BioSimSpace derivation.

           work_dir : str
               The working directory for the simulation.

           cutoff: float
               The greatest distance between ligand and receptor anchor atoms, in
               Angstrom. Receptor anchors further than cutoff Angstroms from the closest
               ligand anchors will not be included in the search for potential anchor points.

            restraint_idx: int
                The index of the restraint from a list of candidate restraints ordered by
                suitability. restraint_idx != 0 is only valid if method == 'BSS'.

           Returns
           -------

           restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
               The restraints of `rest_type` which best mimic the strongest receptor-ligand
               interactions.

        """
        if method == "MDRestraintsGenerator":
            if restraint_idx != 0:
                raise ValueError("restraint_idx must be 0 for MDRestraintsGenerator.")
            if is_MDRestraintsGenerator:
                return RestraintSearch._boresch_restraint_MDRestraintsGenerator(
                    u, system, temperature, lig_selection_str,
                    recept_selection_str, work_dir)
            else:
                raise ImportError('MDRestraintsGenerator not available.')

        elif method == "BSS":
            return RestraintSearch._boresch_restraint_BSS(
                u, system, temperature, lig_selection_str,
                recept_selection_str, work_dir, cutoff,
                restraint_idx=restraint_idx)

    @staticmethod
    def _boresch_restraint_MDRestraintsGenerator(u, system, temperature, lig_selection_str,
                           recept_selection_str, work_dir):
        """Generate the Boresch Restraint using MDRestraintsGenerator.

       Parameters
       ----------

       u : MDAnalysis.Universe
           The trajectory for the ABFE restraint calculation as a
           MDAnalysis.Universe object.

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The molecular system for the ABFE restraint calculation. This
           must contain a single decoupled molecule and is assumed to have
           already been equilibrated.

       temperature : :class:`System <BioSimSpace.Types.Temperature>`
           The temperature of the system

       lig_selection_str: str
           The selection string for the atoms in the ligand to consider
           as potential anchor points.

       recept_selection_str: str
           The selection string for the protein in the ligand to consider
           as potential anchor points.

       work_dir : str
           The working directory for the simulation.

       Returns
       -------

       restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
           The restraints of `rest_type` which best mimic the strongest receptor-ligand
           interactions.
        """
        ligand_atoms = _search.find_ligand_atoms(
            u,
            l_selection=lig_selection_str,
            p_align=recept_selection_str)
        # find protein atoms
        atom_set = []
        for l_atoms in ligand_atoms:
            psearch = _search.FindHostAtoms(
                u,
                l_atoms[0],
                p_selection=recept_selection_str)
            psearch.run()
            atom_set.extend(
                [(l_atoms, p) for p in psearch.host_atoms])
        # Create the boresch finder analysis object
        boresch = _FindBoreschRestraint(u, atom_set)
        # Run the restraint analysis
        boresch.run()

        # Save the analysis results
        boresch.restraint.plot(path=work_dir)
        # Write out the intermolecular section to a topology
        boresch.restraint.write(path=work_dir)
        dG_off = boresch.restraint.standard_state()
        with open(f'{work_dir}/dG_off.dat', 'w') as writer:
            writer.write(f'{dG_off}')

        # This just shows how is the data being recorded, so the
        # index gets reassigned multiple times.
        # https://github.com/IAlibay/MDRestraintsGenerator/blob/master/MDRestraintsGenerator/datatypes.py#L810-L825
        l1_idx, r1_idx = boresch.restraint.bond.atomgroup.atoms.ix
        l2_idx, l1_idx, r1_idx = boresch.restraint.angles[0].atomgroup.atoms.ix
        l1_idx, r1_idx, r2_idx = boresch.restraint.angles[1].atomgroup.atoms.ix
        l3_idx, l2_idx, l1_idx, r1_idx = boresch.restraint.dihedrals[
            0].atomgroup.atoms.ix
        l2_idx, l1_idx, r1_idx, r2_idx = boresch.restraint.dihedrals[
            1].atomgroup.atoms.ix
        l1_idx, r1_idx, r2_idx, r3_idx = boresch.restraint.dihedrals[
            2].atomgroup.atoms.ix

        # The index of the best frame
        index = boresch.restraint.min_frame
        # r1-l1 (r0, kr)
        r0 = boresch.restraint.bond.values[index] * _angstrom
        kr = 10 * _kcal_per_mol / (_angstrom ** 2)
        # r2-r1-l1 (thetaA0, kthetaA)
        thetaA0 = boresch.restraint.angles[1].values[index] * _degree
        kthetaA = 10 * _kcal_per_mol / (_radian ** 2)
        # r1-l1-l2 (thetaB0, kthetaB)
        thetaB0 = boresch.restraint.angles[0].values[index] * _degree
        kthetaB = 10 * _kcal_per_mol / (_radian ** 2)
        # r3-r2-r1-l1 (phiA0, kphiA)
        phiA0 = boresch.restraint.dihedrals[2].values[index] * _degree
        kphiA = 10 * _kcal_per_mol / (_radian ** 2)
        # r2-r1-l1-l2 (phiB0, kphiB)
        phiB0 = boresch.restraint.dihedrals[1].values[index] * _degree
        kphiB = 10 * _kcal_per_mol / (_radian ** 2)
        # r1-l1-l2-l3 (phiC0, kphiC)
        phiC0 = boresch.restraint.dihedrals[0].values[index] * _degree
        kphiC = 10 * _kcal_per_mol / (_radian ** 2)

        restraint_dict = {
            # The default index is in the format of numpy.int64
            # So need to convert to int
            "anchor_points": {"r1": system.getAtom(int(r1_idx)),
                              "r2": system.getAtom(int(r2_idx)),
                              "r3": system.getAtom(int(r3_idx)),
                              "l1": system.getAtom(int(l1_idx)),
                              "l2": system.getAtom(int(l2_idx)),
                              "l3": system.getAtom(int(l3_idx))},
            "equilibrium_values": {"r0": r0,
                                   "thetaA0": thetaA0,
                                   "thetaB0": thetaB0,
                                   "phiA0": phiA0,
                                   "phiB0": phiB0,
                                   "phiC0": phiC0},
            "force_constants": {"kr": kr,
                                "kthetaA": kthetaA,
                                "kthetaB": kthetaB,
                                "kphiA": kphiA,
                                "kphiB": kphiB,
                                "kphiC": kphiC
                                }}
        # TODO: extract the best frame and feed it into Restraint
        # Waiting for the BSS to fix the getFrames
        # best_frame = traj.getFrames(index)
        best_frame = system
        restraint = Restraint(best_frame, restraint_dict,
                              temperature,
                              rest_type='Boresch')
        return restraint

    @staticmethod
    def _boresch_restraint_BSS(u, system, temperature, lig_selection_str,
                           recept_selection_str, work_dir, cutoff, restraint_idx=0):
        """Generate the Boresch Restraint.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE restraint calculation. This
            must contain a single decoupled molecule and is assumed to have
            already been equilibrated.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        lig_selection_str: str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        recept_selection_str: str
            The selection string for the protein in the ligand to consider
            as potential anchor points.

        work_dir : str
            The working directory for the simulation.

        cutoff: float
            The greatest distance between ligand and receptor anchor atoms, in
            Angstrom. Receptor anchors further than cutoff Angstroms from the closest
            ligand anchors will not be included in the search for potential anchor points.

        restraint_idx: int
            The index of the restraint from a list of candidate restraints ordered by
            increasing configurational volume - a lower index gives a stronger restraint.

        Returns
        -------

        restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
            The restraints of `rest_type` which best mimic the strongest receptor-ligand
            interactions.
        """
        

        def findOrderedPairs(u, lig_selection_str, recept_selection_str, cutoff):
            """Return a list of receptor-ligand anchor atoms pairs in the form
            (lig atom index, receptor atom index), where the pairs are ordered
            from low to high variance of distance over the trajectory.

            Parameters
            ----------

            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.

            lig_selection_str: str
                The selection string for the atoms in the ligand to consider
                as potential anchor points.

            recept_selection_str: str
                The selection string for the protein in the ligand to consider
                as potential anchor points.

            cutoff: float
                The greatest distance between ligand and receptor anchor atoms, in
                Angstrom. Receptor anchors further than cutoff Angstroms from the closest
                ligand anchors will not be included in the search for potential anchor points.

            Returns
            -------

            pairs_ordered_sd : list of tuples
                List of receptor-ligand atom pairs ordered by increasing variance of distance over
                the trajectory.
            """

            lig_selection = u.select_atoms(lig_selection_str)
            pair_variance_dict = {}

            # Get all receptor atoms within specified distance of cutoff
            for lig_atom in lig_selection:
                for prot_atom in u.select_atoms(
                        f"{recept_selection_str} and (around {cutoff} index {lig_atom.index})"):
                    pair_variance_dict[(lig_atom.index, prot_atom.index)] = {}
                    pair_variance_dict[(lig_atom.index, prot_atom.index)]["dists"] = []

            # Compute Average Distance and SD
            for frame in tqdm(u.trajectory, desc="Searching for low variance pairs. Frame no: "):
                for lig_atom_index, prot_atom_index in pair_variance_dict.keys():
                    distance = _dist(_mda.AtomGroup([u.atoms[lig_atom_index]]),
                                    _mda.AtomGroup([u.atoms[prot_atom_index]]),
                                    box=frame.dimensions)[2][0]
                    pair_variance_dict[(lig_atom_index, prot_atom_index)][
                        "dists"].append(distance)

            # change lists to numpy arrays
            for pair in pair_variance_dict.keys():
                pair_variance_dict[pair]["dists"] = _np.array(
                    pair_variance_dict[pair]["dists"])

            # calculate SD
            for pair in pair_variance_dict.keys():
                pair_variance_dict[pair]["sd"] = pair_variance_dict[pair]["dists"].std()

            # get n pairs with lowest SD
            pairs_ordered_sd = []
            for item in sorted(pair_variance_dict.items(),
                            key=lambda item: item[1]["sd"]):
                pairs_ordered_sd.append(item[0])

            return pairs_ordered_sd


        def getAnchorAts(a1_idx, u):
            """Takes in index of anchor atom 1 (in either the receptor or ligand)
            and universe and returns list of all three anchor atoms, which are chosen
            to be contiguous and not H.

            Parameters
            ----------

            a1_idx : int
                Index of the first anchor atom

            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.

            Returns
            -------

            inds : List of ints
                The indices of all three selected anchor atoms to be used in 
                Boresch restraints.
            """
            a1_at = u.atoms[a1_idx]
            bonded_heavy_at = a1_at.bonded_atoms.select_atoms("not name H*")
            a2_idx = bonded_heavy_at[0].index

            if len(bonded_heavy_at) > 1:
                # not at end of chain
                a3_idx = bonded_heavy_at[1].index
                # Might be better to return all possible combinations
            else:
                # at end of chain, get next heavy atom along
                a3_idx = \
                bonded_heavy_at[0].bonded_atoms.select_atoms("not name H*")[
                    0].index

            return a1_idx, a2_idx, a3_idx


        def getDistance(idx1, idx2, u):
            """ Distance in Angstrom"""
            distance = _dist(_mda.AtomGroup([u.atoms[idx1]]),
                             _mda.AtomGroup([u.atoms[idx2]]),
                             box=u.dimensions)[2][0]
            return distance


        def getAngle(idx1, idx2, idx3, u):
            """Angle in rad"""
            C = u.atoms[idx1].position
            B = u.atoms[idx2].position
            A = u.atoms[idx3].position
            angle = _calc_angles(C, B, A, box=u.dimensions)
            return angle


        def getDihedral(idx1, idx2, idx3, idx4, u):
            """Dihedral in rad"""
            positions = [u.atoms[idx].position for idx in
                         [idx1, idx2, idx3, idx4]]
            dihedral = _calc_dihedrals(positions[0], positions[1],
                                       positions[2], positions[3],
                                       box=u.dimensions)
            return dihedral


        def getBoreschDof(l1, l2, l3, r1, r2, r3, u):
            """Calculate Boresch degrees of freedom from indices of anchor atoms"""
            # Ordering of connection of anchors is r3,r2,r1,l1,l2,l3
            r = getDistance(r1, l1, u)
            thetaA = getAngle(r2, r1, l1, u)
            thetaB = getAngle(r1, l1, l2, u)
            phiA = getDihedral(r3, r2, r1, l1, u)
            phiB = getDihedral(r2, r1, l1, l2, u)
            phiC = getDihedral(r1, l1, l2, l3, u)
            # Not restrained but distance from collinearity must be checked
            thetaR = getAngle(r3, r2, r1, u)  # Receptor internal angle
            thetaL = getAngle(l1, l2, l3, u)  # Ligand internal angle
            return r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL


        def getConfigVol(equil_vals, force_consts, temp):
            """Find the configurational volume accessible to the 
            decoupled restrained ligand based on the Boresch restraint.
            Based on part of Eqn 32 of J. Phys. Chem. B 2003, 107, 35, 9535–9551.
            
            Parameters
            ----------

            equil_vals : dict
                Dictionary of equilibrium values for the Boresch restraint. Must have units
                of Angstrom and radians. Of the form {"r": r0, "thetaA": thetaA0, 
                "thetaB": thetaB0, "phiA": phiA0, "phiB": phiB0, "phiC": phiC0}.

            force_consts : dict
                Dictionary of force constants for the Boresch restraint. Force constants
                must have units of kcal /mol. Of the form {"r": kr, "thetaA": kthetaA,
                "thetaB": kthetaB, "phiA": kphiA0, "phiB": kphiB, "phiC": kphiC}.

            temp : float
                The temperature, in K.

            Returns
            -------

            config_vol : float
                The configurational volume accessible to the restrained decoupled ligand,
                in Angstrom^3.
            """
            RT = k_boltz * temp # in kcal / mol
            numerator1 = (equil_vals["r"] ** 2) * _np.sin(equil_vals["thetaA"]) * \
                            _np.sin(equil_vals["thetaB"]) # Units: A**2
            numerator2 = (2 * _np.pi * RT) **3 # Units: (kcal / mol )**3
            denominator = _np.sqrt(_np.array([val for val in force_consts.values()]).prod()) 
            # Units: (kcal / mol)**3 A**-1
            config_vol = numerator1 * numerator2 / denominator # Units: A**3
            return config_vol


        def findOrderedBoresch(u, pair_list, temp, no_pairs=50):
            """Calculate a list of Boresch restraints and associated 
            statistics over the trajectory.

            Parameters
            ----------

            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.

            pair_list : List of tuples
                List of receptor-ligand atom pairs to be used as the r1 and l1
                anchor points in candidate Boresch restraints.

            temp : float
                The temperature, in K.

            no_pairs : int
                Number of pairs to be used in the calculation. Pairs in pair_list
                after the index supplied will be ignored.

            Returns
            -------

            pairs_ordered_boresch : List of tuples
                List of receptor-ligand atom pairs from pair_list ordered by priority
                the Boresch restraints they generate. pairs_ordered_boresch[0] labels
                the optimal restraint according to the algorithm implemented.

            boresch_dof_data: dict
                Dictionary of statistics for the Boresch restraints obtained over the
                trajectory. Keys are the pair tuples supplied in pair_list.
            """
            boresch_dof_list = ["r", "thetaA", "thetaB", "phiA", "phiB",
                                "phiC", "thetaR", "thetaL"] #thetaR and thetaL are the internal
                                                            # angles of the receptor and ligand

            # get values of degrees of freedom for lowest SD pairs across whole trajectory

            boresch_dof_data = {}
            for pair in tqdm(pair_list[:no_pairs], desc="Scoring candidate Boresch anchor points. Anchor set no: "):
                boresch_dof_data[pair] = {}
                l1_idx, r1_idx = pair
                _, l2_idx, l3_idx = getAnchorAts(l1_idx, u)
                _, r2_idx, r3_idx = getAnchorAts(r1_idx, u)
                boresch_dof_data[pair]["anchor_ats"] = [l1_idx, l2_idx, l3_idx,
                                                        r1_idx, r2_idx, r3_idx]


                # Add sub dictionaries for each Boresch degree of freedom
                for dof in boresch_dof_list:
                    boresch_dof_data[pair][dof] = {}
                    boresch_dof_data[pair][dof]["values"] = []

                for i, _ in enumerate(
                        u.trajectory):  # TODO: Use MDA.analysis.base instead?
                    r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL = getBoreschDof(
                        l1_idx, l2_idx, l3_idx, r1_idx, r2_idx, r3_idx, u)
                    boresch_dof_data[pair]["r"]["values"].append(r)
                    boresch_dof_data[pair]["thetaA"]["values"].append(thetaA)
                    boresch_dof_data[pair]["thetaB"]["values"].append(thetaB)
                    boresch_dof_data[pair]["phiA"]["values"].append(phiA)
                    boresch_dof_data[pair]["phiB"]["values"].append(phiB)
                    boresch_dof_data[pair]["phiC"]["values"].append(phiC)
                    boresch_dof_data[pair]["thetaR"]["values"].append(thetaR)
                    boresch_dof_data[pair]["thetaL"]["values"].append(thetaL)

                # Calculate statistics for each Boresch degree of freedom
                for dof in boresch_dof_list:
                    boresch_dof_data[pair][dof]["values"] = _np.array(
                        boresch_dof_data[pair][dof]["values"])
                    # Check not dihedral
                    if not dof[:3] == "phi":
                        boresch_dof_data[pair][dof]["avg"] = \
                        boresch_dof_data[pair][dof]["values"].mean()
                        boresch_dof_data[pair][dof]["var"] = \
                        boresch_dof_data[pair][dof]["values"].var()
                    # If dihedral, have to calculate circular stats
                    else:
                        circmean = _circmean(boresch_dof_data[pair][dof]["values"], 
                                    high=_np.pi, low=-_np.pi)
                        boresch_dof_data[pair][dof]["avg"] = circmean

                        # Cannot use scipy's circvar as later than v 1.8
                        # as this is calculated in the range 0 - 1
                        corrected_values = []
                        for val in boresch_dof_data[pair][dof]["values"]:
                            dtheta = abs(val - circmean)
                            corrected_values.append(
                                min(dtheta, 2 * _np.pi - dtheta))
                        corrected_values = _np.array(
                            corrected_values)
                        boresch_dof_data[pair][dof][
                            "var"] = corrected_values.var()

                    # Assume Gaussian distributions and calculate force constants for harmonic potentials
                    # so as to reproduce these distributions at 298 K
                    boresch_dof_data[pair][dof]["k"] = k_boltz * temp / (
                                boresch_dof_data[pair][dof][
                                    "var"])  # Force constants in kcal mol-1 A-2 [rad-2]

                # Calculate the configurational volume accessible based on each restraint
                equil_vals = {dof:boresch_dof_data[pair][dof]["avg"] for dof in boresch_dof_list}
                force_consts = {dof:boresch_dof_data[pair][dof]["k"] for dof in boresch_dof_list}
                boresch_dof_data[pair]["config_vol"] = getConfigVol(equil_vals, force_consts, temp)

            # Order pairs according to configurational volume - smaller volume indicates stronger
            # restraints mimicking stronger native interactions
            pairs_ordered_boresch_var = []
            for item in sorted(boresch_dof_data.items(),
                            key=lambda item: item[1]["config_vol"]):
                pairs_ordered_boresch_var.append(item[0])

            # Filter out r < 1, theta >150 or < 30
            pairs_ordered_boresch = []
            for pair in pairs_ordered_boresch_var:
                # Might be an improvement to filter by kT to collinearity, rather than distance,
                # although this also could be problematic when the distrubutions are far
                # from Gaussian
                cond_dist = boresch_dof_data[pair]["r"]["avg"] > 1
                avg_angles = []
                angles = ["thetaA",
                        "thetaB"]  # May also be good to check internal angles, although will be much stiffer
                for angle in angles:
                    avg_angles.append(boresch_dof_data[pair][angle]["avg"])
                cond_angles = list(
                    map(lambda x: (x < 2.62 and x > 0.52), avg_angles)) # 150 and 30 degrees
                if cond_dist and all(cond_angles):
                    pairs_ordered_boresch.append(pair)

            if len(pairs_ordered_boresch) == 0:
                raise _AnalysisError(
                    "No candidate sets of Boresch restraints are suitable. Please expand "
                    "search criteria.")
            
            return pairs_ordered_boresch, boresch_dof_data


        def plotDOF(ordered_restraint_labels, dof_data, restraint_idx=0,
                    dof_to_plot=["r", "thetaA", "thetaB", "phiA", "phiB","phiC"]):
            """Plot historgrams and variation with time of DOF over a trajectory.

            Parameters
            ----------

            ordered_restraint_labels : list
                List of labels for candidate restraints, which will be used to 
                look up the restraints in dof_data. Presumed to be ordered by
                some measure of desirability for use.

            dof_data : dict
                Dictionary of data (obtained from the trajectory) for each restraint
                with a label given in ordered_restraint_labels.

            restraint_idx : int
                Index of the restraint in ordered_restraint_labels to use.

            dof_to_plot : list
                List of DOF to plot.
            """
            # The labels for each DOF to use on the plots
            dof_labels = {"r": r"$r$ / $\mathrm{\AA}$", "thetaA": r"$\theta_A$ / rad", "thetaB": r"$\theta_B$ / rad",
                          "phiA": r"$\phi_A$ / rad", "phiB": r"$\phi_B$ / rad", "phiC": r"$\phi_C$ / rad"}

            n_dof = len(dof_to_plot)
            label = ordered_restraint_labels[restraint_idx]

            # Plot histograms
            fig, axs = _plt.subplots(1, n_dof, figsize=(16, 4), dpi=500)
            for i, dof in enumerate(dof_to_plot):
                axs[i].hist(dof_data[label][dof]["values"], bins=10)
                axs[i].axvline(x=dof_data[label][dof]["avg"], color='r',
                            linestyle='dashed', linewidth=2, label="mean")
                axs[i].set_xlabel(dof_labels[dof])
                axs[i].set_ylabel("Num Vals")
                axs[i].legend()
            fig.tight_layout()
            fig.savefig(f'{work_dir}/restraint_idx{restraint_idx}_dof_hist.png', facecolor="white")

            # Plot variation with time to see if there are slow DOF
            fig, axs = _plt.subplots(1, n_dof, figsize=(16, 4), dpi=500)
            for i, dof in enumerate(dof_to_plot):
                axs[i].plot([x for x in range(len(dof_data[label][dof]["values"]))],
                            dof_data[label][dof]["values"])
                axs[i].set_ylabel(dof_labels[dof])
                axs[i].set_xlabel("Frame No")
            fig.tight_layout()
            fig.savefig(f'{work_dir}/restraint_idx{restraint_idx}_dof_time.png', facecolor="white")


        def getBoreschRestraint(pair, boresch_dof_data):
            """Get the Boresch restraints for a specified pair in a form compatible
            with BSS.

            Parameters
            ----------

            pair : tuple
                The (receptor_idx, ligand_idx) anchor atom pair labelling
                the restraint in boresch_dof_data

            restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
                The restraint defined by boresch_dof_data and labelled by pair.

            """
            anchor_idxs = {'l1': boresch_dof_data[pair]["anchor_ats"][0],
                           'l2': boresch_dof_data[pair]["anchor_ats"][1],
                           'l3': boresch_dof_data[pair]["anchor_ats"][2],
                           'r1': boresch_dof_data[pair]["anchor_ats"][3],
                           'r2': boresch_dof_data[pair]["anchor_ats"][4],
                           'r3': boresch_dof_data[pair]["anchor_ats"][5]}

            anchor_ats = {k: system.getAtom(int(v)) for k, v in anchor_idxs.items()}

            # Check we have found all anchors
            if not len(anchor_ats) == 6:
                raise _AnalysisError(
                    "Could not find all anchor atoms in system")

            # Get remaining parameters
            r0 = boresch_dof_data[pair]["r"]["avg"]
            thetaA0 = boresch_dof_data[pair]["thetaA"]["avg"]
            thetaB0 = boresch_dof_data[pair]["thetaB"]["avg"]
            phiA0 = boresch_dof_data[pair]["phiA"]["avg"]
            phiB0 = boresch_dof_data[pair]["phiB"]["avg"]
            phiC0 = boresch_dof_data[pair]["phiC"]["avg"]
            kr = boresch_dof_data[pair]["r"]["k"]
            kthetaA = boresch_dof_data[pair]["thetaA"]["k"]
            kthetaB = boresch_dof_data[pair]["thetaB"]["k"]
            kphiA = boresch_dof_data[pair]["phiA"]["k"]
            kphiB = boresch_dof_data[pair]["phiB"]["k"]
            kphiC = boresch_dof_data[pair]["phiC"]["k"]

            restraint_dict = {
                "anchor_points": {"r1": anchor_ats["r1"],
                                  "r2": anchor_ats["r2"],
                                  "r3": anchor_ats["r3"],
                                  "l1": anchor_ats["l1"],
                                  "l2": anchor_ats["l2"],
                                  "l3": anchor_ats["l3"]},
                "equilibrium_values": {"r0": r0 * _Units.Length.angstrom,
                                       "thetaA0": thetaA0 * _radian,
                                       "thetaB0": thetaB0 * _radian,
                                       "phiA0": phiA0 * _radian,
                                       "phiB0": phiB0 * _radian,
                                       "phiC0": phiC0 * _radian},
                "force_constants": {"kr": kr * _kcal_per_mol / _angstrom ** 2,
                                    "kthetaA": kthetaA * _kcal_per_mol / (
                                                _radian * _radian),
                                    "kthetaB": kthetaB * _kcal_per_mol / (
                                                _radian * _radian),
                                    "kphiA": kphiA * _kcal_per_mol / (
                                                _radian * _radian),
                                    "kphiB": kphiB * _kcal_per_mol / (
                                                _radian * _radian),
                                    "kphiC": kphiC * _kcal_per_mol / (
                                                _radian * _radian)}}

            restraint =  Restraint(system, restraint_dict, temperature=temperature, rest_type='Boresch')
            return restraint


        # Find pairs with lowest SD
        pairs_ordered_sd = findOrderedPairs(u, lig_selection_str, recept_selection_str, cutoff)

        # Convert to Boresch anchors, order by correction, and filter
        pairs_ordered_boresch, boresch_dof_data = findOrderedBoresch(u, pairs_ordered_sd, temperature.value())

        # Plot
        plotDOF(pairs_ordered_boresch, boresch_dof_data, restraint_idx = restraint_idx)

        # Convert to BSS compatible dictionary
        restraint = getBoreschRestraint(pairs_ordered_boresch[restraint_idx], boresch_dof_data)

        return restraint  # TODO: implement normal frame - waiting for traj.getFrames() to be fixed.

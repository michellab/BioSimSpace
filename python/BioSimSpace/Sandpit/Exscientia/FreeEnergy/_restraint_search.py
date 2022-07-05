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
from MDAnalysis.lib.distances import calc_dihedrals as _calc_dihedrals
from numpy.linalg import norm as _norm
import matplotlib.pyplot as _plt

from Sire.Base import getBinDir as _getBinDir
from Sire.Base import getShareDir as _getShareDir

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
    _engines = ["GROMACS", ] # TODO: "AMBER", "SOMD"

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
        self._process.isRunning()

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
                cutoff=10, # In Angstrom
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
                cutoff=cutoff)

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
                cutoff=10): # In Angstrom
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
                recept_selection_str, method, work_dir, cutoff)

    @staticmethod
    def _boresch_restraint(u, system, temperature, lig_selection_str,
                           recept_selection_str, method, work_dir, cutoff):
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

           Returns
           -------

           restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
               The restraints of `rest_type` which best mimic the strongest receptor-ligand
               interactions.

        """
        if method == "MDRestraintsGenerator":
            if is_MDRestraintsGenerator:
                return RestraintSearch._boresch_restraint_MDRestraintsGenerator(
                    u, system, temperature, lig_selection_str,
                    recept_selection_str, work_dir)
            else:
                raise ImportError('MDRestraintsGenerator not available.')
        # No need to review this part as Finlay is still working on it.
        elif method == "BSS":
            return RestraintSearch._boresch_restraint_BSS(
                u, system, temperature, lig_selection_str,
                recept_selection_str, work_dir, cutoff)

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
                           recept_selection_str, work_dir, cutoff):
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

       Returns
       -------

       restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
           The restraints of `rest_type` which best mimic the strongest receptor-ligand
           interactions.
        """
        # Please don't review this part yet as Finlay is working on it.
        # TODO Tidy this up and improve. Change algorithm to fit to match that discussed.
        lig_selection = u.select_atoms(lig_selection_str)

        # anchors dict of dict. For each ligand heavy atom there is a dictionary of protein heavy atoms,
        # for each of which there is a dictionary of average distances and standard deviation

        anchors_dict = {}
        for lig_atom in lig_selection:
            for prot_atom in u.select_atoms(
                    f"{recept_selection_str} and (around {cutoff} index {lig_atom.index})"):
                anchors_dict[(lig_atom.index, prot_atom.index)] = {}
                anchors_dict[(lig_atom.index, prot_atom.index)]["dists"] = []

        ### Compute Average Distance and SD

        for frame in u.trajectory:
            for lig_atom_index, prot_atom_index in anchors_dict.keys():
                distance = _dist(_mda.AtomGroup([u.atoms[lig_atom_index]]),
                                 _mda.AtomGroup([u.atoms[prot_atom_index]]),
                                 box=frame.dimensions)[2][0]
                anchors_dict[(lig_atom_index, prot_atom_index)][
                    "dists"].append(distance)

        # change lists to numpy arrays
        for pair in anchors_dict.keys():
            anchors_dict[pair]["dists"] = _np.array(
                anchors_dict[pair]["dists"])

        # calculate average and SD
        for pair in anchors_dict.keys():
            anchors_dict[pair]["avg_dist"] = anchors_dict[pair]["dists"].mean()
            anchors_dict[pair]["sd_dist"] = anchors_dict[pair]["dists"].std()

        # get n pairs with lowest SD
        pairs_ordered_sd = []
        for item in sorted(anchors_dict.items(),
                           key=lambda item: item[1]["sd_dist"]):
            pairs_ordered_sd.append(item[0])
            # print(f'Pair: {item[0]}, av distance: {item[1]["avg_dist"]:.2f}, SD: {item[1]["sd_dist"]:.2f}')

        # Print out pairs with lowest SD
        # print("The ligand-protein atom pairs with the lowest SD in distance are:")
        # for i in range(5):
        #    print(f"{u.atoms[pairs_ordered_sd[i][0]]} and {u.atoms[pairs_ordered_sd[i][1]]}")

        ### For Pairs with Lowest Pairwise RMSDs, find Adjacent Heavy Atoms

        def getAnchorAts(a1_idx, u):
            """Takes in index of anchor atom 1 and universe and returns
            list of all three anchor atoms, which are chosen to be bonded
            and not H"

            Args:
                a1_idx (int): Index of the first anchor atom
                u (mda universe): The mda universe

            Returns:
                ints: The indices of all three anchor points
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

        ### Use These As Anchors and Plot Variance of Associated Degrees of Freedom

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
            BA = A - B
            BC = C - B
            angle = _np.arccos(_np.dot(BA, BC) / (_norm(BA) * _norm(BC)))
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
            # Not restrained but distance from coolinearity must be checked
            thetaR = getAngle(r3, r2, r1, u)  # Receptor internal angle
            thetaL = getAngle(l1, l2, l3, u)  # Ligand internal angle
            return r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL

        lig_anchors = getAnchorAts(pairs_ordered_sd[0][0], u)
        prot_anchors = getAnchorAts(pairs_ordered_sd[0][1], u)
        getBoreschDof(lig_anchors[0], lig_anchors[1], lig_anchors[2],
                      prot_anchors[0], prot_anchors[1], prot_anchors[2], u)

        # get values of degrees of freedom for lowest SD pairs across whole trajectory

        boresch_dof_dict = {}
        for pair in pairs_ordered_sd[:200]:  # Check top 200 pairs
            boresch_dof_dict[pair] = {}
            l1_idx, r1_idx = pair
            _, l2_idx, l3_idx = getAnchorAts(l1_idx, u)
            _, r2_idx, r3_idx = getAnchorAts(r1_idx, u)
            boresch_dof_dict[pair]["anchor_ats"] = [l1_idx, l2_idx, l3_idx,
                                                    r1_idx, r2_idx, r3_idx]

            boresch_dof_list = ["r", "thetaA", "thetaB", "phiA", "phiB",
                                "phiC", "thetaR", "thetaL"]

            # Add sub dictionaries for each Boresch degree of freedom
            for dof in boresch_dof_list:
                boresch_dof_dict[pair][dof] = {}
                boresch_dof_dict[pair][dof]["values"] = []

            # Populate these dictionaries with values from trajectory
            n_frames = len(u.trajectory)

            for i, frame in enumerate(
                    u.trajectory):  # TODO: Use MDA.analysis.base instead?
                r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL = getBoreschDof(
                    l1_idx, l2_idx, l3_idx, r1_idx, r2_idx, r3_idx, u)
                boresch_dof_dict[pair]["r"]["values"].append(r)
                boresch_dof_dict[pair]["thetaA"]["values"].append(thetaA)
                boresch_dof_dict[pair]["thetaB"]["values"].append(thetaB)
                boresch_dof_dict[pair]["phiA"]["values"].append(phiA)
                boresch_dof_dict[pair]["phiB"]["values"].append(phiB)
                boresch_dof_dict[pair]["phiC"]["values"].append(phiC)
                boresch_dof_dict[pair]["thetaR"]["values"].append(thetaR)
                boresch_dof_dict[pair]["thetaL"]["values"].append(thetaL)

                if i == n_frames - 1:
                    boresch_dof_dict[pair]["tot_var"] = 0
                    for dof in boresch_dof_list:
                        boresch_dof_dict[pair][dof]["values"] = _np.array(
                            boresch_dof_dict[pair][dof]["values"])
                        boresch_dof_dict[pair][dof]["avg"] = \
                        boresch_dof_dict[pair][dof]["values"].mean()
                        # For dihedrals, compute variance and mean based on list of values corrected for periodic boundary at
                        # pi radians, because there is no problem with dihedrals in this region
                        if dof[:3] == "phi":
                            avg = boresch_dof_dict[pair][dof]["avg"]

                            # correct variance - fully rigorous
                            corrected_values_sd = []
                            for val in boresch_dof_dict[pair][dof]["values"]:
                                dtheta = abs(val - avg)
                                corrected_values_sd.append(
                                    min(dtheta, 2 * _np.pi - dtheta))
                            corrected_values_sd = _np.array(
                                corrected_values_sd)
                            boresch_dof_dict[pair][dof][
                                "sd"] = corrected_values_sd.std()

                            # Correct mean (will fail if very well split above and below 2pi)
                            # get middle of interval based on current mean
                            # TODO: Should use circular mean instead, see scipy
                            corrected_values_avg = []
                            periodic_bound = avg - _np.pi
                            if periodic_bound < -_np.pi:
                                periodic_bound += 2 * _np.pi
                            # shift vals from below periodic bound to above
                            for val in boresch_dof_dict[pair][dof]["values"]:
                                if val < periodic_bound:
                                    corrected_values_avg.append(
                                        val + 2 * _np.pi)
                                else:
                                    corrected_values_avg.append(val)
                            corrected_values_avg = _np.array(
                                corrected_values_avg)
                            mean_corrected = corrected_values_avg.mean()
                            # shift mean back to normal range
                            if mean_corrected > _np.pi:
                                boresch_dof_dict[pair][dof][
                                    "avg"] = mean_corrected - 2 * _np.pi
                            else:
                                boresch_dof_dict[pair][dof][
                                    "avg"] = mean_corrected

                        else:
                            boresch_dof_dict[pair][dof]["sd"] = \
                            boresch_dof_dict[pair][dof]["values"].std()
                        # Exclude variance of internal angles as these are not restrained
                        if (dof != "thetaR" and dof != "thetaL"):
                            boresch_dof_dict[pair]["tot_var"] += \
                            boresch_dof_dict[pair][dof]["sd"] ** 2
                        # Assume Gaussian distributions and calculate force constants for harmonic potentials
                        # so as to reproduce these distributions
                        boresch_dof_dict[pair][dof]["k"] = 0.593 / (
                                    boresch_dof_dict[pair][dof][
                                        "sd"] ** 2)  # RT at 289 K is 0.593 kcal mol-1

        ### Filter, Pick Optimum Degrees of Freedom, and Select Force Constants Based on Variance, Select Equilibrium Values

        # Order pairs according to variance
        pairs_ordered_boresch_var = []
        for item in sorted(boresch_dof_dict.items(),
                           key=lambda item: item[1]["tot_var"]):
            pairs_ordered_boresch_var.append(item[0])
            # print(f'Pair: {item[0]}, total variance: {boresch_dof_dict[item[0]]["tot_var"]}')

        # Filter out r < 1, theta >150 or < 30
        selected_pairs_boresch = []
        for pair in pairs_ordered_boresch_var:
            cond_dist = boresch_dof_dict[pair]["r"]["avg"] > 1
            avg_angles = []
            # angles = ["thetaA", "thetaB", "thetaR","thetaL"] # also check internal angles
            angles = ["thetaA",
                      "thetaB"]  # May also be good to check internal angles, although will be much stiffer
            for angle in angles:
                avg_angles.append(boresch_dof_dict[pair][angle]["avg"])
            cond_angles = list(
                map(lambda x: (x < 2.62 and x > 0.52), avg_angles))
            if cond_dist and all(cond_angles):
                selected_pairs_boresch.append(pair)

        # Plotting
        dof_to_plot = ["r", "thetaA", "thetaB", "phiA", "phiB",
                       "phiC"]  # Can also plot thetaR and thetaL
        n_dof = len(dof_to_plot)
        pair = selected_pairs_boresch[0]

        # Plot histograms
        fig, axs = _plt.subplots(1, n_dof, figsize=(16, 4), dpi=500)
        for i, dof in enumerate(dof_to_plot):
            axs[i].hist(boresch_dof_dict[pair][dof]["values"], bins=10)
            axs[i].axvline(x=boresch_dof_dict[pair][dof]["avg"], color='r',
                           linestyle='dashed', linewidth=2, label="mean")
            if dof == "r":
                axs[i].set_xlabel("r ($\AA$)")
            else:
                axs[i].set_xlabel(f"{dof} (rad)")
            axs[i].set_ylabel("Num Vals")
            axs[i].legend()
        fig.tight_layout()
        fig.savefig(f'{work_dir}/boresch_dof_hist.png', facecolor="white")

        # Plot variation with time to see if there are slow DOF
        fig, axs = _plt.subplots(1, n_dof, figsize=(16, 4), dpi=500)
        for i, dof in enumerate(dof_to_plot):
            axs[i].plot([x for x in range(300)],
                        boresch_dof_dict[pair][dof]["values"])
            if dof == "r":
                axs[i].set_ylabel("r ($\AA$)")
            else:
                axs[i].set_ylabel(f"{dof} (rad)")
            axs[i].set_xlabel("Frame No")
        fig.tight_layout()
        fig.savefig(f'{work_dir}/boresch_dof_time.png', facecolor="white")

        # Print out Boresch parameters
        def getBoreschRestraint(pair):
            anchor_idxs = {'l1': boresch_dof_dict[pair]["anchor_ats"][0],
                           'l2': boresch_dof_dict[pair]["anchor_ats"][1],
                           'l3': boresch_dof_dict[pair]["anchor_ats"][2],
                           'r1': boresch_dof_dict[pair]["anchor_ats"][3],
                           'r2': boresch_dof_dict[pair]["anchor_ats"][4],
                           'r3': boresch_dof_dict[pair]["anchor_ats"][5]}

            anchor_ats = {}

            # Get atoms from system, excluding water and ions
            for mol in system.search(
                    "not (resname WAT) or (resname CL) or (resname NA)"):  # TODO: Expand this
                for anchor in list(
                        anchor_idxs.keys()):  # Change to list to avoid error due to changing dict size
                    search = mol.search(f'atomidx {anchor_idxs[anchor]}')
                    # see if we have found the anchor
                    if len(search) == 1:
                        anchor_ats[anchor] = search[0]
                        del anchor_idxs[anchor]

            # Check we have found all anchors
            if not len(anchor_ats) == 6:
                raise _AnalysisError(
                    "Could not find all anchor atoms in system")

            # Get remaining parameters
            r0 = boresch_dof_dict[pair]["r"]["avg"]
            thetaA0 = boresch_dof_dict[pair]["thetaA"]["avg"]
            thetaB0 = boresch_dof_dict[pair]["thetaB"]["avg"]
            phiA0 = boresch_dof_dict[pair]["phiA"]["avg"]
            phiB0 = boresch_dof_dict[pair]["phiB"]["avg"]
            phiC0 = boresch_dof_dict[pair]["phiC"]["avg"]
            kr = boresch_dof_dict[pair]["r"]["k"]
            kthetaA = boresch_dof_dict[pair]["thetaA"]["k"]
            kthetaB = boresch_dof_dict[pair]["thetaB"]["k"]
            kphiA = boresch_dof_dict[pair]["phiA"]["k"]
            kphiB = boresch_dof_dict[pair]["phiB"]["k"]
            kphiC = boresch_dof_dict[pair]["phiC"]["k"]

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

            return Restraint(system, restraint_dict, type='Boresch')

        restraint = getBoreschRestraint(selected_pairs_boresch[0])
        normal_frame = None

        return restraint  # TODO: implement normal frame

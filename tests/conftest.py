from pathlib import Path

collect_ignore_glob = ["*/out_test*.py"]

import os

import BioSimSpace as BSS

from BioSimSpace._Utils import _try_import, _have_imported

# Store the tutorial URL.
url = BSS.tutorialUrl()

# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None

# Make sure AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

# Make sure NAMD is installed.
try:
    from sire.legacy.Base import findExe

    findExe("namd2")
    has_namd = True
except:
    has_namd = False

# Check whether AMBER parameterisation executables are installed.
has_tleap = BSS.Parameters._Protocol._amber._tleap_exe is not None
has_antechamber = BSS.Parameters._Protocol._amber._antechamber_exe is not None

# Check if openff is installed.
_openff = _try_import("openff")
has_openff = _have_imported(_openff)

# Check for MDAnalysis.
mda = _try_import("MDAnalysis")
has_mdanalysis = _have_imported(mda)

# Check for MDTraj.
mdtraj = _try_import("mdtraj")
has_mdtraj = _have_imported(mdtraj)


# Check for MDRestraintsGenerator.
MDRestraintsGenerator = _try_import(
    "MDRestraintsGenerator",
    install_command="pip install MDRestraintsGenerator",
)

has_mdrestraints_generator = _have_imported(MDRestraintsGenerator)

root_fp = Path(__file__).parent.resolve()

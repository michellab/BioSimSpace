import BioSimSpace as BSS
import os
import shutil
from toPertFile import _toPertFile

if os.path.exists("runs/"):
	shutil.rmtree("runs/")


ethane = BSS.IO.readMolecules(["ethane.prm7", "ethane.rst7"])[0]
methanol = BSS.IO.readMolecules(["methanol.prm7", "methanol.rst7"])[0]

mapping = BSS.Align.matchAtoms(ethane, methanol)

# Align ethane to methanol based on the mapping.
ethane = BSS.Align.rmsdAlign(ethane, methanol, mapping)

print("Merging..")
# Merge the ethane and methanol based on the mapping.
merged = BSS.Align.merge(ethane, methanol, mapping)

print("Writing pert file..")
# write pert file:
for pert_type in ["discharge_soft",
				"vanish_soft",
				#"change_bonds",
				"change_hard",
				"grow_dummies",
				"charge_dummies"]:
	print(pert_type)
	_toPertFile(merged, pert_type=pert_type)

# print("Solvating..")
# solvated = BSS.Solvent.tip3p(molecule=merged, box=3*[3*BSS.Units.Length.nanometer])

# print("Setting protocol..")
# protocol = BSS.Protocol.FreeEnergy(
# 			timestep=2*BSS.Units.Time.femtosecond, 
# 			runtime=4*BSS.Units.Time.nanosecond, 
# 			num_lam=5
# 			)
# print("Initiating environment..")
# freenrg = BSS.FreeEnergy.Solvation(solvated, protocol, work_dir="runs/proto_1", multistep=True)



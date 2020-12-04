import BioSimSpace as BSS
import Sire
import os
import shutil
from toPertFile import _toPertFile

# if os.path.exists("runs/"):
# 	shutil.rmtree("runs/")


ethane = BSS.IO.readMolecules(["ethane.prm7", "ethane.rst7"])[0]
methanol = BSS.IO.readMolecules(["methanol.prm7", "methanol.rst7"])[0]

mapping = BSS.Align.matchAtoms(ethane, methanol)

# Align ethane to methanol based on the mapping.
ethane = BSS.Align.rmsdAlign(ethane, methanol, mapping)

print("Merging..")
# Merge the ethane and methanol based on the mapping.
merged = BSS.Align.merge(ethane, methanol, mapping)

print("Writing pert files..")
# write pert file:
for pert_type in [
				"standard",
				"discharge_soft",
				"vanish_soft",
				"change_bonds",
				"change_hard",
				"grow_soft",
				"charge_soft"
				]:

	pert_mol = _toPertFile(merged, pert_type=pert_type)

	# Create a composite system
	system1 = BSS._SireWrappers.System(BSS._SireWrappers.Molecule(pert_mol))

	print("Solvating", pert_type+"..")
	solvated = BSS.Solvent.tip3p(molecule=system1, box=3*[3*BSS.Units.Length.nanometer])
	
	BSS.IO.saveMolecules("SOMD_INPUTS/"+pert_type, solvated, "PRM7")
	BSS.IO.saveMolecules("SOMD_INPUTS/"+pert_type, solvated, "RST7")


	# if pert_type == "standard":
	# 	num_lam = 18
	# else:
	# 	num_lam = 3

	# print("Setting protocol..")
	# protocol = BSS.Protocol.FreeEnergy(
	# 			timestep=2*BSS.Units.Time.femtosecond, 
	# 			runtime=4*BSS.Units.Time.nanosecond, 
	# 			num_lam=num_lam
	# 			)
	# print("Initiating environment..")
	# freenrg = BSS.FreeEnergy.Solvation(solvated, protocol, work_dir="SOMD/ethane_methanol/"+pert_type)


# ethane = BSS.IO.readMolecules(["ethane.prm7", "ethane.rst7"])[0]
# methanol = BSS.IO.readMolecules(["methanol.prm7", "methanol.rst7"])[0]

# mapping = BSS.Align.matchAtoms(methanol, ethane)

# # Align ethane to methanol based on the mapping.
# methanol = BSS.Align.rmsdAlign(methanol, ethane, mapping)

# print("Merging..")
# # Merge the ethane and methanol based on the mapping.
# merged = BSS.Align.merge(methanol, ethane, mapping)

# print("Writing pert file..")
# # write pert file:
# for pert_type in ["discharge_soft",
# 				"vanish_soft",
# 				"change_bonds",
# 				"change_hard",
# 				"grow_soft",
# 				"charge_soft"]:
# 	print(pert_type)
# 	_toPertFile(merged, pert_type=pert_type)

# print("Solvating..")
# solvated = BSS.Solvent.tip3p(molecule=merged, box=3*[3*BSS.Units.Length.nanometer])

# print("Setting protocol..")
# protocol = BSS.Protocol.FreeEnergy(
# 			timestep=2*BSS.Units.Time.femtosecond, 
# 			runtime=1*BSS.Units.Time.nanosecond, 
# 			num_lam=3
# 			)
# print("Initiating environment..")
# freenrg = BSS.FreeEnergy.Solvation(solvated, protocol, work_dir="SOMD/mol-eth_onestep_3w")
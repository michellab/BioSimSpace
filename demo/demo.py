
import BioSimSpace as BSS

print("Loading molecules...")
complex = BSS.readMolecules( ["proteinbox.crd", "proteinbox.top"] )

print("Extracting the old ligand...")
old_ligand = complex[ BSS.MolWithResName("ZAN") ]
complex.remove(old_ligand.number())

print("Reading in the new ligand...")
new_ligand = BSS.readMolecule( ["ligand.crd", "ligand.top"] )

print("Finding the maximum common substructure between the ligands...")
print("...and using this to align the ligands")
new_ligand = new_ligand.move().align(old_ligand,BSS.MCSMatcher()).commit()

print("Adding the new ligand to the complex...")
complex.add(new_ligand)

print("Pretend that we are doing MD now to equilibrate...")

print("Writing the complex to a set of output files, using the ")
print("same fileformat as the original files.")
filenames = BSS.saveMolecules( "smallest", complex, complex.fileFormat() )

print("Names of saved files are %s" % filenames)

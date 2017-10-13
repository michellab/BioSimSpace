
import BioSimSpace as BSS

print("\nLoading molecules...")
complex = BSS.readMolecules( ["proteinbox.crd", "proteinbox.top"] )

print("\nExtracting the old ligand...")
old_ligand = complex.take( BSS.MolWithResName("ZAN") )

print("\nReading in the new ligand...")
new_ligand = BSS.readMolecule( ["ligand.crd", "ligand.top"] )

print("\nFinding the maximum common substructure between the ligands...")
print("...and using this to align the ligands")
new_ligand = new_ligand.move().align(old_ligand,BSS.MCSMatcher()).commit()

print("\nAdding the new ligand to the complex...")
complex.add(new_ligand)

print("\nPretend that we are doing MD now to equilibrate...")

print("\nWriting the complex to a set of output files, using the ")
print("same fileformat as the original files.")
filenames = BSS.saveMolecules( "smallest", complex, complex.fileFormat() )

print("\nNames of saved files are %s" % filenames)

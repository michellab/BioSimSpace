# BioSimSpace.Parameters

This package provides functionality for parameterising molecular systems.
Parameterisation is achieved by wrapping the `tLEaP` package from
[AmberTools](http://ambermd.org/AmberTools.php) and `pbd2gmx` from
[GROMACS](http://www.gromacs.org).

As an example:

```python
import BioSimSpace as BSS

# Load a molecular system from file.
system = BSS.IO.readMolecules("molecule.pdb")

# Extract the first molecule from the system.
molecule = system.getMolecules()[0]

# Get a list of the available force fields.
BSS.Parameters.forceFields()

# Parameterise using the Amber03 force field.
ff03_mol = BSS.Parameters.ff03(molecule).getMolecule()

# Parameterise using the Generalized Amber Force Field.
gaff_mol = BSS.Parameters.gaff(molecule).getMolecule()
```

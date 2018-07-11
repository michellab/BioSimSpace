# BioSimSpace.Solvent

This package provides functionality for solvating molecular systems. This is
achieved by wrapping the `gmx solvate` program from [GROMACS](http://www.gromacs.org).

As an example:

```python
import BioSimSpace as BSS

# Load a molecular system from file.
system = BSS.IO.readMolecules("molecule.pdb")

# Extract the first molecule from the system.
molecule = system.getMolecules()[0]

# Solvate the molecule, centred in a 5 nanometer box.
solvated_system = BSS.Solvent.tip3p(molecule, box=3 * [5 * BSS.Units.Length.nanometer])

# Solvate the molecule, centred in a 5 nanometer box with a 1 nanometer
# shell of solvent around the molecule.
solvated_system = BSS.Solvent.tip3p(molecule, box=3 * [5 * BSS.Units.Length.nanometer],
                                              shell=BSS.Units.Length.nanometer)


# Solvate the molecule, centred in a 5 nanometer box with an ionic strength
# of 0.1 mol / litre and neutralised.
solvated_system = BSS.Solvent.tip3p(molecule, box=3 * [5 * BSS.Units.Length.nanometer],
                                              ion_conc=0.1, is_neutral=True)
```

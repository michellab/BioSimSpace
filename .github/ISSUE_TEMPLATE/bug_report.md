---
name: Bug report
about: Create a report to help us improve
title: ''
labels: 'bug'
assignees: ''
---

# Prerequisites

If this is your first time posting an issue please take time to read the
information in this page.

Please answer the following questions for yourself before submitting an issue.

**YOU MAY DELETE THE PREREQUISITES SECTION.**

- [ ] I am running the latest version (devel branch, or dev Conda package).
- [ ] I checked the documentation and found no answer.
- [ ] I checked to make sure that this issue has not already been filed.

# Expected Behavior

Please describe the behavior you are expecting. (Where relevant, use
screenshots and upload files to help illustrate this.)

# Current Behavior

What is the current behavior?

# Failure Information

Please help provide information about the failure. (Where relevant, use screenshots and upload files to help illustrate this.)

## Steps to Reproduce

Please provide detailed steps for reproducing the issue. Make use of markdown code blocks, e.g., if BioSimSpace failed to parameterise a molecule:

```python
import BioSimSpace as BSS

# Load the system and extract the only molecule.
system = BSS.IO.readMolecules(["molecule.prm7", "molecule.rst7"])
molecule = system[0]

# Create a background process to parameterise the molecule.
process = BSS.Parameters.gaff(molecule)

# Get the parameterised molecule from the process.
molecule = process.getMolecule()
```

## Failure logs

BioSimSpace provides built in functionality for extracting log files and error messages from background processes. Please examine these for pertinent information.

For the example above, we could retrieve all of the intermediate log files from the parameterisation process as follows:

```python
import BioSimSpace as BSS

# Load the system and extract the only molecule.
system = BSS.IO.readMolecules(["molecule.prm7", "molecule.rst7"])
molecule = system[0]

# Create a background process to parameterise the molecule.
process = BSS.Parameters.gaff(molecule)

# Get the parameterised molecule from the process.
molecule = process.getMolecule()

# Get the intermediate files from the background process.
process.getOutput(filename="param_eror.zip")
```

Alternatively, the process could be run in a user directory. (By default, background processes are run in a temporary directory which is destroyed when the Python kernel exits.)

```python
import BioSimSpace as BSS

# Load the system and extract the only molecule.
system = BSS.IO.readMolecules(["molecule.prm7", "molecule.rst7"])
molecule = system[0]

# Create a background process to parameterise the molecule.
# Run the process in a directory called "my_directory".
process = BSS.Parameters.gaff(molecule, work_dir="my_directory")

# Get the parameterised molecule from the process.
molecule = process.getMolecule()
```

## Context

Please provide any relevant information about your setup. This is important in case the issue is not reproducible except for under certain conditions.

* Operating System: Linux, macOS.
* Installation method: conda, binary, source. (If installed from source, please
also give the [Sire](https://github.com/michellab/Sire) version number.)
* If you're running on a cluster: What is the scheduler? What type of nodes are you using, e.g. number of CPUs, GPUS, amount of RAM?
* What molecular system are you simulating?. If possible, please attach topology and coordinate files to help us debug. If you can't, please give basic information, such as the system size, which will help us estimate the expected memory footprint.

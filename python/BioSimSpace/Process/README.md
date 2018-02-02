# BioSimSpace.Process

This sub-package provides functionality for running different molecular
simulation processes.

The `process.py` module defines a set of common helper functions and a
base class, `Process`, that defines common properties and methods for
all processes. Individual modules, such as `namd.py`, defined functionality
for running process with a particular software package. At present we
provide support for [NAMD](http://www.ks.uiuc.edu/Research/namd) and
[AMBER](http://ambermd.org).

## Object instantiation

All process classes must take at least two arguments to their constructor:

* `system`: A Sire molecular system.

* `protocol`: A [`BioSimSpace.Protocol`](../Protocol) object defining the
protocol for the simulation process, e.g. an equilibration protocol. A
user who wishes to run a custom process can simply replace the protocol
argument with the path to an appropriate configuration file.

For example, to initialise an object to run a default minimistion protocol
using NAMD:

```python
import BioSimSpace as BSS

# Create a list of the input files.
files = ["alanin.psf", "alanin.pdb", "alanin.params"]

# Create a molecular system.
system = BSS.readMolecules(files)

# Create a default minimisation protocol.
protocol = BSS.Protocol.Minimisation()

# Initialise the NAMD process.
proc = BSS.Process.Namd(system, protocol)
```

To use a custom protocol, the constructor could be called as follows:

```python
# Initialise the NAMD process using a custom protocol.
proc = BSS.Process.Namd(system, protocol="config.namd")
```

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
using AMBER:

```python
import BioSimSpace as BSS

# Create a molecular system.
system = BSS.readMolecules(["ala.crd", "ala.top"])

# Create a default minimisation protocol.
protocol = BSS.Protocol.Minimisation()

# Initialise the AMBER process.
process = BSS.Process.Amber(system, protocol)
```

To use a custom protocol, the constructor could be called as follows:

```python
# Initialise the AMBER process using a custom protocol.
process = BSS.Process.Amber(system, protocol="config.amber")
```

By default, each process is run in a temporary workspace. To specify
the working directory the user can pass an appropriate keyword argument:

```python
# Initialise the AMBER process using a custom working directory.
process = BSS.Process.Amber(system, work_dir="/my/custom/path")
```

The directory will be created if it doesn't already exist (assuming write
priviledes on the path).

Once initialised, the process object will have set up all of the appropriate
input and configuration files needed to run the desired simulation protocol
using AMBER.

To get a list of the auto-generated input files:

```python
# Return a list of names for the input files.
files = process.inputFiles()
```

To get a list of the configuration file options:

```python
# Return the contents of the configuration file as a list of strings.
config = process.getConfig()
```

BioSimSpace uses a set of well chosen configuration parameters as defaults for
each simulation protocol. However, we provide lots of flexibility for overriding
these defaults. For example:

```python
# Add a single additional configuration string.
param = "some-parameter = some-value"
process.addToConfig(param)

# Add a list of additional configuration parameters.
params = ["some-parameter = some-value", "some-other-parameter = "some-other-value"]
process.addToConfig(params)

# Add some additional parameters from a file.
process.addToConfig("params.txt")

# Overwrite the entire configuration using a new set of parameters.
process.setConfig(params)         # Using a list of parameter strings.
process.setConfig("params.txt")   # Using a parameter file.

# Write the current configuration parameters to a file.
process.writeConfig("params.txt")
```

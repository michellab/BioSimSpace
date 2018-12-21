# BioSimSpace

[![Build Status](https://travis-ci.org/michellab/BioSimSpace.svg?branch=testing)](https://travis-ci.org/michellab/BioSimSpace)

Code and resources for the [EPSRC](https://epsrc.ukri.org)
[BioSimSpace](https://biosimspace.org) project.

BioSimSpace is an interoperable framework for biomolecular simulation. With it you
can:

- Write robust and portable biomolecular workflow components that work on
different hardware, with different software packages, and that can be run
in different ways, e.g. command-line, [Jupyter](http://jupyter.org).
- Interact with running molecular simulation processes in real-time within
a Jupyter notebook.

## Installation

The [`install.sh`](install.sh) script will automatically install and configure
a clean BioSimSpace application on Linux or macOS. Following a successful
installation, the following aliases will be set:

* `bss_python` - The BioSimSpace python interpreter.
* `bss_ipython` - The BioSimSpace interactive python shell.
* `bss_jupyter` - The BioSimSpace Jupyter executable.
* `bss_test` - Run the BioSimSpace unit tests.

(You'll need to run `source ~/.biosimspacerc` or restart your shell for these
to take effect.)

Note that this defaults to installing the latest stable version of BioSimSpace
from the [`master`](https://github.com/michellab/BioSimSpace/tree/master)
branch. To install the development version, run:

```bash
./install.sh --branch=devel
```

Since the Sire binary is quite large (~700 MB), it is only re-downloaded
whenever a new version is released. To force a clean install with a fresh
download, run:

```bash
./install.sh --clean
```

BioSimSpace runs many subprocess in temporary working directories. To control
the directory location, set the `TMPDIR` environment variable.

Full installation details are given below for those who wish to install
BioSimSpace manually, or already have built their own Sire application.

---

BioSimSpace is built on top of the [Sire](https://siremol.org) molecular
simulations framework. At present you will need to need to install BioSimSpace
into an existing [Sire](https://siremol.org/pages/download.html) application:

```bash
cd python
python setup.py install
```

(You can safely ignore warnings regarding `mmtf-python`, `mkl-random`, and `mkl-fft`.)

Note that this assumes that Sire's `bin` directory, typically `$HOME/sire.app/bin`,
has been added to your path. Alternatively, run:

```bash
$HOME/sire.app/bin/python setup.py install
```
At present we recommend using the `devel` branch of [Sire](https://github.com/michellab/Sire)
to ensure that you have access to the latest features and bug-fixes.

If you experience problems with [`matplotlib`](https://matplotlib.org) when
importing BioSimSpace on macOS, e.g.

```python
RuntimeError**: Python is not installed as a framework.
```

simply add the following to `~/.matplotlib/matplotlibrc`

```bash
backend: TkAgg
```

Note that plotting functionality will be disabled if you are using BioSimSpace
within IPython on a remote server without X forwarding.

Several additional packages are required for full access to all of BioSimSpace's
functionality. Please download and install these packages according to their
recommended installation instructions.

* [Amber / AmberTools](http://ambermd.org): _Dynamics / Parameterisation / Solvation_
* [Gromacs](http://www.gromacs.org/): _Dynamics / Parameterisation / Solvation_
* [Somd](https://siremol.org/tutorials/somd/): _Dynamics_ (bundled)
* [Namd](http://www.ks.uiuc.edu/Research/namd/): _Dynamics_ (optional)

(Installation and use of these packages is subject to the terms of their
respective licenses.)

For Somd support you'll need to export the following (Bash style) environment
variable:

```bash
# Use the appropriate path for your Sire installation.
OPENMM_PLUGIN_DIR=$HOME/sire.app/lib/plugins/
```

## Documentation

Each package has its own README page:

- [BioSimSpace.Align](python/BioSimSpace/Align)
- [BioSimSpace.FreeEnergy](python/BioSimSpace/FreeEnergy)
- [BioSimSpace.Gateway](python/BioSimSpace/Gateway)
- [BioSimSpace.IO](python/BioSimSpace/IO)
- [BioSimSpace.MD](python/BioSimSpace/MD)
- [BioSimSpace.Notebook](python/BioSimSpace/Notebook)
- [BioSimSpace.Parameters](python/BioSimSpace/Parameters)
- [BioSimSpace.Process](python/BioSimSpace/Process)
- [BioSimSpace.Protocol](python/BioSimSpace/Protocol)
- [BioSimSpace.Solvent](python/BioSimSpace/Solvent)
- [BioSimSpace.Trajectory](python/BioSimSpace/Trajectory)
- [BioSimSpace.Types](python/BioSimSpace/Types)
- [BioSimSpace.Units](python/BioSimSpace/Units)

Full API documentation and examples can be found at [biosimspace.org](https://biosimspace.org).

## Dependencies

BioSimSpace makes use of several external python packages. These should be
automatically installed into your Sire app by the [setup.py](python/setup.py)
script.

## Demos

A collection of example scripts and notebooks are included in the [demo](demo)
directory.

## Developing

BioSimSpace is written in Python but uses an object-oriented style with C++
naming conventions. More details can be found [here](python).

Please create a feature branch for development work. Branches related to a
particular feature should be prefixed with `feature`, i.e. for the
BioSimSpace.Process module:

```bash
feature-process
```

When you are happy with your feature, merge into the
[`testing`](https://github.com/michellab/BioSimSpace/tree/testing) branch to
check that the continuous integration tests pass. When ready, create a pull
request so that the feature can be merged into devel. Once the merge is
successful, please delete the redundant feature branch.

## Tests

We currently don't run continuous integration (CI) against all BioSimSpace branches.
This is because BioSimSpace is built on top of [Sire](https://siremol.org) meaning
that a binary version of Sire must be downloaded and unpacked on the CI server, into
which BioSimSpace is installed.

Since the Sire binary is large, executing a fresh download for every commit is
currently prohibitively expensive. Instead, we use a
[`testing`](https://github.com/michellab/BioSimSpace/tree/testing) branch for CI.
Feature branches are periodically merged into `testing` to validate the build
process and check that all unit tests pass.

For those who wish to run things locally, a test suite can be found
[here](https://github.com/michellab/BioSimSpaceUnitTests).

## Issues

Please report bugs and other issues using the GitHub [issue tracker](https://github.com/michellab/BioSimSpace/issues).
When reporting issues please try to include a minimal code snippet that reproduces
the problem. Additional files can be also be uploaded as an archive, e.g. a zip
file. Please also report the branch on which you are experiencing the issue,
along with the BioSimSpace version number. This can be found by running:

```python
import BioSimSpace as BSS
print(BSS.__version__)
```

## Versioning

BioSimSpace uses the following versioning convention for releases, `YYYY.MAJOR.MINOR`,
e.g. `2018.1.0`. For a summary of the version history, please refer to the [CHANGELOG](CHANGELOG.md).

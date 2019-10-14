---
title: 'BioSimSpace: An interoperable Python framework for biomolecular simulation.'
tags:
  - biomolecular-simulation
  - computational-chemistry
  - computational-physics
  - computational-biology
  - interoperability
  - molecular dynamics
  - reproducibility
authors:
  - name: Lester Hedges
    orcid: 0000-0002-5624-0500
    affiliation: 1
  - name: Antonia Mey
    orcid: 0000-0001-7512-5252
    affiliation: 2
  - name: Christopher Woods
    orcid: 0000-0001-6563-9903
    affiliation: 1
  - name: Julien Michel
    orcid: 0000-0003-0360-1760
    affiliation: 2
affiliations:
 - name: Advanced Computing Research Centre, University of Bristol
   index: 1
 - name: Department of Chemistry, University of Edinburgh
   index: 2
date: \today
bibliography: paper.bib
---

# Summary

In most research communities there is not a single, unified software-framework. Instead, researchers are presented with a collection of competing packages from which they pick and choose the functionality that they desire. Interoperability between packages is often poor, with incompatibilities between file formats, hardware, etc. This leads to brittle workflows, poor reproducibility, and lock in to specific software. For the biomolecular simulation community, our solution has been the introduction of an interoperable framework that collects together the core functionality of many packages and exposes it through a simple Python API. By not choosing to reinvent the wheel, we can take advantage of all the fantastic software that has already been written, and can easily plug into new software packages as they appear. Our software can convert between many common molecular file formats and automatically find packages available within the environment on which it is run. This allows allows the user to write portable workflow components that can be run with different input, on different environments, and in completely different ways, e.g. from the command-line, or within a [Jupyter](https://jupyter.org) notebook running on a cloud server.

# Molecular dynamics

One of the core features of BioSimSpace is the ability to set up and run molecular dynamics (MD) simulations. There are a large number of packages that can run molecular dynamics for biomolecules and BioSimSpace supports several of these: [AMBER](http://ambermd.org), [GROMACS](http://www.gromacs.org), and [NAMD](https://www.ks.uiuc.edu/Research/namd). BioSimSpace also comes with a bundled GPU accelerated molecular dynamics engine, SOMD, so there is always a fall back in case no other packages are installed.

While, broadly speaking, the different molecular dynamics engines offer a similar range of features, their interfaces are quite different. This makes it hard to take expertise in one package and immediately apply it to another. At the heart of this problem is the incompatibility between the molecular file formats used by the different packages. While they all contain the same information, i.e. how atoms are laid out in space and how they interact with each other, the structure of the files is very different. In order to provide interoperability betwen packages we will need to be able to read and write many different file formats, and be able to interconvert between them too.

# Features

## Parsers

At its core, BioSimSpace is built around a powerful set of file parsers which allow reading and writing of a wide range of molecular file formats. File input/output is provided via the `BioSimSpace.IO` package using parsers from the [Sire](https://siremol.org) molecular simulation framework, on top of which BioSimSpace is built. Unlike many other programs, we take the approach that it is the _contents_ of the file that defines it format, not the _extension_. As such, we attempt to parse a file with all of our parsers in parallel. Any parser for which the contents of the file is incompatible will be rejected early, with the eventual format of the file determined by the parser that completed without error.

Typically, the information needed to construct a molecular system is split across multiple files, e.g. a _coordinate_ file containing the atomic coordinates, and a _topology_ file that describes how the atoms within each molecule are bonded together, along with parameters for the potential of the molecular model. To handle this, each of our parsers are assigned as being able to _lead_, or _follow_, or both. Lead parsers are able to initialise a molecular system (typically by constructing the topology), whereas those that follow can add additional information to an existing molecular system. Lead parsers may also be able to follow, such that when multiple lead parsers are associated with a set of files then the one that ultimately leads will be determined by which lead parser is unable to follow. This approach allows us to easily parse molecular information from multiple files, even if those formats aren't typically associated with each other. As long as the molecular topology corresponding the information in the files is consistent, then they can be read.

As files are parsed, records in those files are assigned to a set of _properties_ that are associated with molecules in the system, e.g. `charge`, `coordinates`, `element`, etc. While some of these properties are unique to particular parsers, others are shared across formats and are converted to a consistent set of internal units on read. Those properties which represent mathematical expressions are stored using Sire's built in computer algebra system. On write, each parser expects molecules in the system to contain a specific set of properties, which are then extracted and converted in order to generate the appropriate records for the format in question. In this way, a bond record from an [AMBER](http://ambermd.org) format file can be read into an internal bond expression, which could then be converted to the appropriate [GROMACS](http://www.gromacs.org) bond record on write.

Another feature of our parsers is guaranteed read/write self-consistency. Any file that can be read can also be written, and vice-versa. In addition, when expected molecular information is missing from a file, unless obvious, we don't attempt to guess what it may have been. In this sense our parsers don't attempt to be _too_ clever, which can lead to unexpected behaviour, particulars when information is modified or supplemented behind the user's back.

The code below shows how to load a set of AMBER format files from a directory:

```python
import BioSimSpace as BSS

system = BSS.IO.readMolecules(BSS.IO.glob("amber/ala/*"))
```

## Protocols and Processes

BioSimSpace simplifies the set-up and running of molecular simulations through an abstraction of simulation _protocols_ and _processes_. Protocols define what a user wants to do to a molecular system, e.g. performing a _minimisation_ to find the local potential energy minimum. Processes are used to apply a protocol to a system using a specific molecular simulation engine.

The `BioSimSpace.Protocol` package provides a set of high-level objects for several common molecular simulation protocols. Each protocol offers as set of configurable options that are handled by all of the molecular simulation engines that we support. `BioSimSpace.Process` provides objects for configuring, running, and interacting with simulation processes for each of the supported engines. When a process is created by passing in a system and protocol, BioSimSpace automatically writes all of the input files required by the specific simulation engine and configures any command-line options required by its executable. (Expert users of a particular engine are free to fully override any of the configuration options if desired.)

The example below shows how to configure and run a default energy minimisation protocol for the molecular system that was loaded earlier. Here we use AMBER as the molecular dynamics engine:

```python
protocol = BSS.Protocol.Minimisation()
process = BSS.Process.Amber(system, protocol)
process.start()
```

## Interoperability

While it is useful to be able to configure and run simulation processes using specific engines, any script written in this way would not be portable since we can't guarantee what software will be available on a different computer. To this end, the `BioSimSpace.MD` package provides functionality for automatically configuring a simulation process using _any_ suitable molecular dynamics engine that is installed on the host system. For example, the AMBER specific example in the previous section can be translated to an interoperable alternative as follows:

```python
protocol = BSS.Protocol.Minimisation()
process = BSS.MD.run(system, protocol)
# By default MD.run starts the process automatically.
```

The `BSS.MD.run` function searches the system for suitable packages that support the chosen protocol, then chooses the most appropriate one to run the simulation. For example, if AMBER was installed then the process returned by `BSS.MD.run` would be of type `BSS.Process.Amber`, if not then the input files could be converted to a different format allowing the use of a different process such as `BSS.Process.Gromacs`.

## Robust and flexible workflow components

The building blocks described above can be used to write interoperable workflow components, or _nodes_. Typically, a node will perform a single, well-defined, unit of work with clear inputs and outputs. The `BioSimSpace.Gatway` package acts as a bridge between BioSimSpace and the outside world, allowing a user to construct a node and define the input and output requirements, along with restrictions on their types and values. As an example, the following code snippet shows how the minimisation example described above can be translated into a node.

```python
import BioSimSpace as BSS

# Initialise the Node object.
node = BSS.Gateway.Node("Minimise a molecular system and save to file.")

# Set the node author and license.
node.addAuthor(name="Lester Hedges",
               email="lester.hedges@bristol.ac.uk",
               affiliation="University of Bristol")
node.setLicense("GPLv3")

# Set the node inputs.
node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
node.addInput("steps", BSS.Gateway.Integer(help="The number of minimisation steps.",
                                           minimum=0, maximum=1000000, default=10000))

# Set the node outputs.
node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system."))

# Show the graphical user interface (GUI) to allow the user to set the inputs.
# This will only happen if running interactively, i.e. in a Jupyter notebook.
node.showControls()

# Load the molecular system using the user defined input "files".
system = BSS.IO.readMolecules(node.getInput("files"))

# Define the minimisation protocol using the user defined number of "steps".
protocol = BSS.Protocol.Minimisation(steps=node.getInput("steps"))

# Execute the process using any available molecular dynamics engine.
process = BSS.MD.run(system, protocol)

# Set the node output to the final configuration of the minimisation process.
# Note that the pass block=True to the getSystem call to ensure that the
# process finished before getting the final configuration. (It is possible
# to query the running process in real time when running interactively.)
# Note that the original file format of the system is preserved on write.
node.setOutput("minimised", BSS.IO.saveMolecules("minimised",
    process.getSystem(block=True), system.fileFormat()))

# Finally, validate the node to make sure that outputs are set correctly
# and no errors have been raised. If running interactively, this will
# generate a download link to a zip file containing the node outputs.
node.validate()
```

BioSimSpace nodes are flexible in the way in which they can be used, with the same script working seamlessly from within a Jupyter notebook or on the command-line. Typically, a user would a write a node as a fully documented, interactive Jupyter notebook, then save it as a regular Python script to run from the command-line. (For inclusion here we simply include the Python script representation of the node, which could be re-converted to a notebook using, e.g., [p2j](https://pypi.org/project/p2j).) Any purely interactive elements included in the node, e.g.  visualisations and plots, are simply ignored when the script is run in a non-interactive mode. To facilitate this dual-use the `node.addInput` method generates a custom [ipywidgets](https://ipywidgets.readthedocs.io/en/latest) based graphical user interface for interative use in Jupyter, or a custom [argparse](https://docs.python.org/3/library/argparse.html) parser for handling command-line arguments. Figure 1 shows the example node above running within a Jupyter notebook (top) and from the command-line (bottom).

![BioSimSpace nodes can be run within a Jupyter notebook (top) or
from the command-line (bottom)](figures/fig1.png)

When working interactively, BioSimSpace also provides functionality for interacting with processes while they are running. This allows the user to monitor the progress of a simulation and generate near real-time plots and visualisations.

## Forwards compatibility

To ensure that BioSimSpace nodes are forwards compatible as new features are added all sub packages can query their own functionality and present this to the user. For example, calling `BioSimSpace.IO.fileFormats()` returns a list of the currently supported molecular file formats, `BioSimSpace.Solvent.waterModels` returns a list of the supported water models, etc. These values can be passed as the `allowed` keyword argument when setting an input requirement of a node, ensuring that the node supports the latest functionality of the package version that is installed. The following code snippet shows a node that can be used to convert to any supported molecular file format, which will continue to work as additional formats are added.

```python
import BioSimSpace as BSS

# Initialise the Node object.
node = BSS.Gateway.Node("Convert between molecular file formats.")

# Set the node author and license.
node.addAuthor(name="Lester Hedges",
               email="lester.hedges@bristol.ac.uk",
               affiliation="University of Bristol")
node.setLicense("GPLv3")

# Set the node inputs.
node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
node.addInput("file_format", BSS.Gateway.String(help="The format to convert to.",
                                                allowed=BSS.IO.fileFormats()))

# Set the node outputs.
node.addOutput("converted", BSS.Gateway.File(help="The converted file."))

# Show the graphical user interface to allow the user to set the inputs.
# This will only happen if running interactively, i.e. in a Jupyter notebook.
node.showControls()

# Load the molecular system using the user defined input "files".
system = BSS.IO.readMolecules(node.getInput("files"))

# Convert the system to the chosen format and set the output.
node.setOutput("converted",
    BSS.IO.saveMolecules("converted", system, node.getInput("file_format")))

# Validate the node.
node.validate()
```

Figure 2 shows how the `allowed=BSS.IO.fileFormats()` argument is translated into a dropdown menu for the Jupyter GUI (top), or using the _choices_ option of argparse on the command-line (bottom). This means that the script is adaptive to the support of additional file parsers in future without need for modification.

![BioSimSpace nodes are adaptive to new functionality without modification.](figures/fig2.png)

# Acknowledgments

This work was funded through an EPSRC Flagship Software grant: EP/P022138/1.

# References

(Probably need to cite AMBER/GROMACS. Not sure whether to use web sites, or papers.)

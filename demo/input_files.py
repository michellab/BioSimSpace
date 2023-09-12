#!/usr/bin/env python
# coding: utf-8

# Author: Lester Hedges<br>
# Email:&nbsp;&nbsp; lester@openbiosim.org
#
# # Input files
#
# In this notebook you'll learn how to use BioSimSpace to generate input files for different molecular dynamics engines.
#
# First we'll need to import BioSimSpace:

# In[ ]:


import BioSimSpace as BSS


# We begin by creating a `Node` object. This is the core of our molecular workflow component. It defines what it does, what input is needed, and the output that is produced.

# In[ ]:


node = BSS.Gateway.Node(
    "A node to create input files for molecular dynamics simulation."
)


# We'll now set the author and license:

# In[ ]:


node.addAuthor(
    name="Lester Hedges",
    email="lester@openbiosim.org",
    affiliation="OpenBioSim",
)
node.setLicense("GPLv3")


# Nodes require inputs. To specify inputs we use the `Gateway` module, which is used as a bridge between BioSimSpace and the outside world. This will allow us to document the inputs, define their type, and specify any constraints on their allowed values. Here we will need a set of files that define the molecular system, a string that will specify the type of simulation protocol that we'd like to perform, and another string to specify the molecular dynamics package. (Note that the strings take specific values, which can be automatically determined by querying BioSimSpace.)

# In[ ]:


node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
node.addInput(
    "protocol",
    BSS.Gateway.String(
        help="The molecular simulation protocol.",
        allowed=BSS.Protocol.protocols(),
        default="Minimisation",
    ),
)
node.addInput(
    "package",
    BSS.Gateway.String(
        help="The molecular dynamics package.",
        allowed=BSS.Process.engines(),
        default="Amber",
    ),
)


# Note that we don't allow the user to specify any further details of the protocol. In this case we use a _best practice_ protocol.
#
# We now need to define the output of the node. In this case we will return the set of input files for the chosen molecular dynamics package.

# In[ ]:


node.addOutput(
    "input_files",
    BSS.Gateway.FileSet(
        help="A zip file containing the molecular dynamics input files."
    ),
)


# When working interactively within a Jupyter notebook we need a way to allow users to set the input requirements. The `node.showControls` method will display a graphical user interface (GUI), from which inputs can be set.
#
# Note that the GUI requires active user input. All input requirements that don't have a default value _must_ be set before the node can proceed. If you try to query the node for one of the user values then an error will be raised. Use the dropdown button to choose options for the protocol and molecular dynamics package.
#
# When working interactively you will typically be running on a remote server where you won't have access to the local filesystem. In this case you'll need to upload files for any of the `File` or `FileSet` input requirements. The GUI below will provide buttons that allow you to browse your own filesystem and select files. Since Jupyter has a limit of 5MB for file transfers, we provide support for compressed formats, such as `.zip` or `.tar.gz`. (A single archive can contain a set of files, allowing you to set a single value for a `FileSet` requirement.) We've provided some example input files that can be used in the training notebooks, which are available to download from the links below. These can then be re-uploaded using the GUI.
#
# AMBER: [ala.crd](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/amber/ala/ala.crd), [ala.top](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/amber/ala/ala.top)
#
# GROMACS: [kigaki.gro](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/gromacs/kigaki/kigaki.gro), [kigaki.top](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/gromacs/kigaki/kigaki.top)
#
# NAMD: [alanin.pdb](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/namd/alanin/alanin.pdb), [alanin.psf](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/namd/alanin/alanin.psf), [alanin.params](https://raw.githubusercontent.com/openbiosim/biosimspace/devel/demo/namd/alanin/alanin.params)
#
# When uploading files the name of the current file(s) will replace the `Upload` button. If you need to change the file, simply click on the button again and choose a new file.

# In[ ]:


node.showControls()


# Once all requirements are set then we can access the values using the `node.getInput` method. The first time this is called the `node` will automatically validate all of the input and report the user if any errors were found.
#
# We'll now create a molecular system using the input files uploaded by the user. Note that we don't specify the format of the files, since this is automatically determined by BioSimSpace. (BioSimSpace has support for a wide range of formats and can convert between certain formats too.)

# In[ ]:


system = BSS.IO.readMolecules(node.getInput("files"))


# Next we setup a molecular dynamics process using the molecular system and the chosen simulation protocol. Note that an error will be thrown if the protocol isn't currently supported by the chosen molecular dynamics package.

# In[ ]:


process = BSS.Process.createProcess(
    system,
    BSS.Protocol.createProtocol(node.getInput("protocol")),
    node.getInput("package"),
)


# Now we extract a zip file containing the autogenerated input files from the process object and assign this to the output variable of the node.

# In[ ]:


node.setOutput("input_files", process.getInput())


# Finally, we validate that the node completed successfully. This will check that all output requirements are satisfied and that no errors were raised by the user. Any file outputs will be available for the user to download as a compressed archive.
#
# Note that the validation will fail until the cell above finishes running.

# In[ ]:


node.validate()


# Once we are satisfied with our node we can choose to download it as a regular Python script that can be run from the command-line.
#
# In JupyterHub, click on: `File/Download As/Python`\
# In JupyterLab, click on: `File/Export Notebook As/Export Notebook to Executable Script`

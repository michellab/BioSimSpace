#!/usr/bin/env python
# coding: utf-8

# Author: Lester Hedges<br>
# Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk
# 
# # Minimisation
# 
# In this notebook you'll learn how to use BioSimSpace to write a robust and interoperable workflow node to perform energy minimisation on a molecular system.
# 
# First we'll need to import BioSimSpace:

# In[ ]:


import BioSimSpace as BSS


# We begin by creating a `Node` object. This is the core of our molecular workflow component. It defines what it does, what input is needed, and the output that is produced.

# In[ ]:


node = BSS.Gateway.Node("A node to perform energy minimisation and save the minimised molecular configuration to file.")


# We'll now set the author and license:

# In[ ]:


node.addAuthor(name="Lester Hedges", email="lester.hedges@bristol.ac.uk", affiliation="University of Bristol")
node.setLicense("GPLv3")


# Nodes require inputs. To specify inputs we use the `Gateway` module, which is used as a bridge between BioSimSpace and the outside world. This will allow us to document the inputs, define their type, and specify any constraints on their allowed values. Here we will need a set of files that define the molecular system, and an integer that indicates the number of minimisation steps to perform.

# In[ ]:


node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
node.addInput("steps", BSS.Gateway.Integer(help="The number of minimisation steps.", minimum=0, maximum=1000000, default=10000))


# Note that the input requirement `steps` has a default value, so it is optional.
# 
# We now need to define the output of the node. In this case we will return a set of files representing the minimised molecular system.

# In[ ]:


node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system."))


# When working interactively within a Jupyter notebook we need a way to allow users to set the input requirements. The `node.showControls` method will display a graphical user interface (GUI), from which inputs can be set.
# 
# Note that the GUI requires active user input. All input requirements that don't have a default value _must_ be set before the node can proceed. If you try to query the node for one of the user values then an error will be raised. For bounded integer inputs you can use a slider to set the value, or type in the input box and press enter.
# 
# When working interactively you will typically be running on a remote server where you won't have access to the local filesystem. In this case you'll need to upload files for any of the `File` or `FileSet` input requirements. The GUI below will provide buttons that allow you to browse your own filesystem and select files. Since Jupyter has a limit of 5MB for file transfers, we provide support for compressed formats, such as `.zip` or `.tar.gz`. (A single archive can contain a set of files, allowing you to set a single value for a `FileSet` requirement.) We've provided some example input files that can be used in the training notebooks, which are available to download from the links below. These can then be re-uploaded using the GUI.
# 
# AMBER: [ala.crd](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/amber/ala/ala.crd), [ala.top](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/amber/ala/ala.top)
# 
# GROMACS: [kigaki.gro](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/gromacs/kigaki/kigaki.gro), [kigaki.top](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/gromacs/kigaki/kigaki.top)
# 
# NAMD: [alanin.pdb](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/namd/alanin/alanin.pdb), [alanin.psf](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/namd/alanin/alanin.psf), [alanin.params](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/namd/alanin/alanin.params)
# 
# 
# When uploading files the name of the current file will replace the `Browse` button. If you need to change the file, simply click on the button again and choose a new file. For `FileSet` requirements, a new `Browse` button will appear whenever an additional file is uploaded.

# In[ ]:


node.showControls()


# Once all requirements are set then we can acces the values using the `node.getInput` method. The first time this is called the `node` will automatically validate all of the input and report the user if any errors were found.
# 
# We'll now create a molecular system using the input files uploaded by the user. Note that we don't specify the format of the files, since this is automatically determined by BioSimSpace. (BioSimSpace has support for a wide range of formats and can convert between certain formats too.)

# In[ ]:


system = BSS.IO.readMolecules(node.getInput("files"))


# In order to run a minimisation we need to define a protocol. This can be done using the `BioSimSpace.Protocol` module. Here we will create a "best practice" minimisation protocol, overriding the number of steps with the input from the user.

# In[ ]:


protocol = BSS.Protocol.Minimisation(steps=node.getInput("steps"))


# We now have everything that is required to run a minimisation. To do so, we use `BioSimSpace.MD` to find an appropriate molecular dynamics package on our current environment. What package is found will depend upon both the system and protocol, as well as the hardware that is available to the user.

# In[ ]:


process = BSS.MD.run(system, protocol)


# When the process has finished running we can get the minimised molecular configuration. We will save this to file using the same format as the original system, and set the output requirement to the list of file names that were written.
# 
# Note that we pass `block=True` to `getSystem`. Since we're working interactively the process is running in the background and `getSystem` will grab and return the _latest_ molecular configuration. By passing `block=True` we wait for the process to finish before getting the _final_ configuration and writing it to file. (You will see that an asterisk appeak next to the cell until the process completes.)

# In[ ]:


node.setOutput("minimised", BSS.IO.saveMolecules("minimised", process.getSystem(block=True), system.fileFormat()))


# Finally, we validate that the node completed succesfully. This will check that all output requirements are satisfied and that no errors were raised by the user. Any file outputs will be available for the user to download as a compressed archive.
# 
# Note that the validation will fail until the cell above finishes running.

# In[ ]:


node.validate()


# Once we are satisfied with our node we can choosed to download it as a regular Python script that can be run from the command-line.
# 
# Click on: `File/Download As/Python`

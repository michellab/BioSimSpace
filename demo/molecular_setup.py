#!/usr/bin/env python
# coding: utf-8

# Author: Lester Hedges<br>
# Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk
# 
# # Molecular setup
# 
# In this notebook you'll learn how to use BioSimSpace to write a robust and interoperable workflow node to generate a molecular system ready for simulation with AMBER. To do so, we'll read in a molecule from a file, e.g. a PDB, parameterise the molecule with a user specified force field, then solvate the molecule with a chosen water model.
# 
# First we'll need to import BioSimSpace:

# In[ ]:


import BioSimSpace as BSS


# We begin by creating a `Node` object. This is the core of our molecular workflow component. It defines what it does, what input is needed, and the output that is produced.

# In[ ]:


node = BSS.Gateway.Node("A node to parameterise and solvate a molecule ready for molecular simulation with AMBER.")


# We'll now set the author and license:

# In[ ]:


node.addAuthor(name="Lester Hedges", email="lester.hedges@bristol.ac.uk", affiliation="University of Bristol")
node.setLicense("GPLv3")


# Nodes require inputs. To specify inputs we use the `Gateway` module, which is used as a bridge between BioSimSpace and the outside world. This will allow us to document the inputs, define their type, and specify any constraints on their allowed values.

# In[ ]:


node.addInput("file", BSS.Gateway.File(help="A molecular input file, e.g. a PDB file."))

node.addInput("forcefield", BSS.Gateway.String(help="The name of the force field to use for parameterisation.",
                                               allowed=BSS.Parameters.forceFields(), default="ff14SB"))

node.addInput("water", BSS.Gateway.String(help="The name of the water model to use for solvation.",
                                          allowed=BSS.Solvent.waterModels(), default="tip3p"))

node.addInput("box_size", BSS.Gateway.Length(help="The base length of the cubic simulation box.", unit="nanometer"))

node.addInput("ion_conc", BSS.Gateway.Float(help="The ionic concentration in mol/litre.",
                                            minimum=0, maximum=1, default=0))


# We also need to define the output of the node. In this case we will return a set of AMBER format files for the parameterised and solvated system.

# In[ ]:


node.addOutput("system", BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER format."))


# When working interactively within a Jupyter notebook we need a way to allow users to set the input requirements. The `node.showControls` method will display a graphical user interface (GUI), from which inputs can be set.
# 
# Note that the GUI requires active user input. All input requirements that don't have a default value _must_ be set before the node can proceed. If you try to query the node for one of the user values then an error will be raised. For bounded  inputs you can use a slider to set the value, or type in the input box and press enter.
# 
# When working interactively you will typically be running on a remote server where you won't have access to the local filesystem. In this case you'll need to upload files for the _file_ input requirement. The GUI below will provide buttons that allow you to browse your own filesystem and select files. Since Jupyter has a limit of 5MB for file transfers, we provide support for compressed formats, such as `.zip` or `.tar.gz`. For simplicity, some example molecule files are available to download from the links below. These can then be re-uploaded using the GUI.
# 
# Example files: [2JJC.pdb](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/molecules/2JJC.pdb), [benzene.pdb](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/molecules/benzene.pdb)
# 
# When uploading files the name of the current file will replace the `Browse` button. If you need to change the file, simply click on the button again and choose a new file.

# In[ ]:


node.showControls()


# Once all requirements are set then we can acces the values using the `node.getInput` method. The first time this is called the `node` will automatically validate all of the input and report the user if any errors were found.
# 
# We'll now create a molecular system using the input file uploaded by the user. Note that we don't specify the format of the file, since this is automatically determined by BioSimSpace. (BioSimSpace has support for a wide range of formats and can convert between certain formats too.)

# In[ ]:


system = BSS.IO.readMolecules(node.getInput("file"))


# We now need to grab the molecule from the system so that we can parameterise it. The system object can be indexed by molecule. Here we grab the first (and only) molecule from the system.

# In[ ]:


molecule = system[0]


# Now let's paramterise the molecule. We do so by calling the `parameterise` function from the `BSS.Parameters` package, passing the molecule and force field name as arguments. Since parameterisation can be slow, the function returns a handle to a process that runs the parameterisation in the background. To get the parameterised molecule from the process we need to call the `getMolecule` method. This is a blocking operation which waits for the process to finish before grabbing the molecule and returning it. 

# In[ ]:


molecule = BSS.Parameters.parameterise(molecule, node.getInput("forcefield")).getMolecule()


# Next we will solvate our molecule in a box of water using the `solvate` function from the `BSS.Solvent` package. This will centre the molecule in a cubic box of a specified size and surround it by water molecules. Here we allow the user to specify the size of the box and the ionic strength (by default, the resulting system will also be neutralised too).
# 
# Note that the molecule is an optional `keyword` argument to the solvate function. This is because its also possible to create a pure water box, i.e. without any molecules in it.

# In[ ]:


system = BSS.Solvent.solvate(node.getInput("water"), molecule=molecule,
                                                     box=3 * [node.getInput("box_size")],
                                                     ion_conc=node.getInput("ion_conc"))


# We now need to define the output of the node. In this case we will return a set of AMBER format files representing the parameterised and solvated molecular system.

# In[ ]:


node.setOutput("system", BSS.IO.saveMolecules("system", system, ["prm7", "rst7"]))


# Finally, we validate that the node completed succesfully. This will check that all output requirements are satisfied and that no errors were raised by the user. Any file outputs will be available for the user to download as a compressed archive.
# 
# Note that the validation will fail until the cell above finishes running.

# In[ ]:


node.validate()


# Once we are satisfied with our node we can choosed to download it as a regular Python script that can be run from the command-line.
# 
# Click on: `File/Download As/Python`

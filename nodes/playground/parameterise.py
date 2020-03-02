#!/usr/bin/env python
# coding: utf-8

# Author: Julien Michel<br>
# Email:&nbsp;&nbsp; julien.michel@ed.ac.uk
# 
# # Parameterise
# 
# This notebook parameterises a molecule and writes a topology usable for a molecular simulation. Based on the Molecular Seup node written by Lester Hedges

# In[ ]:


import BioSimSpace as BSS


# In[ ]:


node = BSS.Gateway.Node("A node to parameterise a molecule ready for molecular simulation.")


# In[ ]:


node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edimnburgh")
node.setLicense("GPLv3")


# In[ ]:


node.addInput("input", BSS.Gateway.File(help="A molecular input file, e.g. a PDB file."))

node.addInput("forcefield", BSS.Gateway.String(help="The name of the force field to use for parameterisation.",
                                               allowed=BSS.Parameters.forceFields(), default="ff14SB"))
node.addInput("output", BSS.Gateway.String(help="The root name of the output files."))


# In[ ]:


node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER format."))


# In[ ]:


node.showControls()


# In[ ]:


system = BSS.IO.readMolecules(node.getInput("input"))


# In[ ]:


molecule = system[0]


# In[ ]:


print (molecule)


# In[ ]:


molecule = BSS.Parameters.parameterise(molecule, node.getInput("forcefield")).getMolecule()


# In[ ]:


BSS.IO.saveMolecules(node.getInput("output"), molecule, ["prm7", "rst7"])


# In[ ]:


node.setOutput("nodeoutput", BSS.IO.saveMolecules(node.getInput("output"), molecule, ["prm7", "rst7"]))


# In[ ]:


node.validate()


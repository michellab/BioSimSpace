
# coding: utf-8

# Author: Julien Michel<br>
# Email:&nbsp;&nbsp; julien.michel@ed.ac.uk
# 
# # Parameterise
# 
# This notebook parameterises a molecule and writes a topology usable for a molecular simulation. Based on the Molecular Seup node written by Lester Hedges

# In[1]:


import BioSimSpace as BSS


# In[2]:


node = BSS.Gateway.Node("A node to parameterise a molecule ready for molecular simulation.")


# In[3]:


node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edimnburgh")
node.setLicense("GPLv3")


# In[4]:


node.addInput("input", BSS.Gateway.File(help="A molecular input file, e.g. a PDB file."))

node.addInput("forcefield", BSS.Gateway.String(help="The name of the force field to use for parameterisation.",
                                               allowed=BSS.Parameters.forceFields(), default="ff14SB"))
node.addInput("output", BSS.Gateway.String(help="The root name of the output files."))


# In[5]:


node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER format."))


# In[6]:


node.showControls()


# In[7]:


system = BSS.IO.readMolecules(node.getInput("input"))


# In[8]:


molecule = system.getMolecules()[0]


# In[9]:


molecule = BSS.Parameters.parameterise(molecule, node.getInput("forcefield")).getMolecule()


# In[ ]:


BSS.IO.saveMolecules(node.getInput("output"), molecule, ["prm7", "rst7"])


# In[10]:


node.setOutput("nodeoutput", BSS.IO.saveMolecules(node.getInput("output"), molecule, ["prm7", "rst7"]))


# In[11]:


node.validate()



# coding: utf-8

# Author: Julien Michel
# email: julien.michel@ed.ac.uk
# 
# # Combine
# 
# Merge two input molecules into a single system

# In[1]:


import BioSimSpace as BSS


# In[2]:


node = BSS.Gateway.Node("A node to solvate a molecule ready for molecular simulation with AMBER.")


# In[3]:


node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edinburgh")
node.setLicense("GPLv3")


# In[4]:


node.addInput("system1", BSS.Gateway.FileSet(help="A topology and coordinates file describing system1"))

node.addInput("system2", BSS.Gateway.FileSet(help="A topology and coordinates file describing system2"))

node.addInput("output", BSS.Gateway.String(help="The root name of the combined system."))


# In[5]:


node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="The combined system"))


# In[6]:


node.showControls()


# In[7]:


system1 = BSS.IO.readMolecules(node.getInput("system1"))


# In[8]:


system2 = BSS.IO.readMolecules(node.getInput("system2"))


# In[9]:


combined = system1 + system2


# In[ ]:


#print (combined)


# In[ ]:


#BSS.IO.saveMolecules(node.getInput("output"), combined, ["prm7", "rst7"])


# In[10]:


node.setOutput("nodeoutput", BSS.IO.saveMolecules(node.getInput("output"), combined, ["prm7", "rst7"]))


# In[11]:


node.validate()


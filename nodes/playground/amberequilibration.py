
# coding: utf-8

# Author: Julien Michel<br>
# Email:&nbsp;&nbsp; julien.michel@ed.ac.uk
# 
# # equilibration
# 
# This notebook runs a 3 step equilibration using AMBER.
# 
# TODO) Make this MD package agnostic

# In[1]:


import BioSimSpace as BSS


# In[2]:


node = BSS.Gateway.Node("A node to equilibrate a solvated molecule.")


# In[3]:


node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edimnburgh")
node.setLicense("GPLv3")


# In[4]:


node.addInput("input", BSS.Gateway.FileSet(help="Topology and coordinate files describing the solvated system"))

node.addInput("output", BSS.Gateway.String(help="The root name of the output files."))


# In[5]:


node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER format."))


# In[6]:


node.showControls()


# # Load input

# In[7]:


system = BSS.IO.readMolecules(node.getInput("input"))


# # Minimisation

# In[8]:


# Initialise a short minimisation protocol.
protocol = BSS.Protocol.Minimisation(steps=1000)


# In[9]:


process = BSS.Process.Amber(system, protocol, name="minimise")


# In[10]:


process.getConfig()


# In[11]:


process.start()


# In[12]:


process.wait()


# In[13]:


minimised = process.getSystem()


# # NVT equilibration (heat up to 300 K with positional restraints on non water molecules)

# In[14]:


# TODO Extend API so can restrain everything but water
# Use loop to run until target T achieved
protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.01, "nanosecond"),
                                      temperature_start=BSS.Types.Temperature(0, "kelvin"), 
                                      temperature_end=BSS.Types.Temperature(300, "kelvin"), 
                                      restrain_backbone=True)


# In[15]:


process = BSS.Process.Amber(minimised, protocol, name="equilibrate-nvt")


# In[16]:


process.getConfig()


# In[17]:


process.start()


# In[18]:


process.wait()


# In[19]:


temps = process.getTemperature(time_series=True)


# In[20]:


# Generate a plot of volume vs temperature.
plotT = BSS.Notebook.plot(process.getTime(time_series=True),
    temps, xlabel="Time (ns)", ylabel="Temperature (K)")


# In[21]:


nvt_equilibrated = process.getSystem()


# # NPT equilibration (300 K, no positional restraints)

# In[48]:


# TODO 
# Use loop to run until target density achieved
# Allow setting of target pressure from constructor
protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.01, "nanosecond"),
                                      temperature=BSS.Types.Temperature(300, "kelvin"), 
                                      ensemble="NPT")


# In[49]:


process = BSS.Process.Amber(nvt_equilibrated, protocol, name="equilibrate-npt")


# ## Override default protocol to wrap particle coordinates in minimum image

# In[50]:


config = process.getConfig()
print (config)


# In[51]:


config[-1] = "iwrap=1,"
config.append(' /')


# In[52]:


print(config)


# In[53]:


process.setConfig(config)


# In[54]:


process.start()


# In[55]:


process.wait()


# In[ ]:


#process.getDensity()


# In[56]:


vols = process.getVolume(time_series=True)


# In[57]:


# Generate a plot of volume vs temperature.
plotV = BSS.Notebook.plot(process.getTime(time_series=True),
    vols, xlabel="Time (ns)", ylabel="Volume (A3^3)")


# In[58]:


npt_equilibrated = process.getSystem()


# # Terminate node

# In[59]:


node.setOutput("nodeoutput", BSS.IO.saveMolecules(node.getInput("output"), npt_equilibrated, ["rst7"]))


# In[60]:


node.validate()


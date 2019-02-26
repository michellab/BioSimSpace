.. _ref-MD:

BioSimSpace.MD
==============

The *MD* package provides functionality for automatically configuring and
running molecular dynamics simulations.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecular system from file.
   system = BSS.IO.readMolecules(BSS.IO.glob("amber/ala/*")

   # Create a default minimisation protocol.
   protocol = BSS.Protocol.Minimisation()

   # Find a molecular dynamics package on the host system that supports the
   # system and protocol defined above. If a package exists, BioSimSpace
   # will auto-generate all of the required input files and return a handle
   # to a process that can run the simulation.
   process = BSS.MD.run(system, protocol)

   # Now start the process in the background.
   process.start()

.. automodule:: BioSimSpace.MD

.. toctree::
   :maxdepth: 1

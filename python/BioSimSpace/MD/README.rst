
BioSimSpace.MD
==============

This package provides functionality for configuring and running various
types of molecular dynamics simulation processes.

At present, we provide support for finding molecular dynamics packages and
automatically configuring and starting processes for the user.
``BioSimSpace.MD`` acts as a wrapper around `BioSimSpace.Process <../Process>`_\ ,
meaning that the user doesn't need to specify the molecular dynamics
package that they wish to use. This will automatically be determined from
the file format of the molecular configuration, the type of simulation
protocol that was chosen, and the hardware resources that are available.

As an example:

.. code-block:: python

   import BioSimSpace as BSS

   # Create a molecular system.
   system = BSS.IO.readMolecules(["ala.crd", "ala.top"])

   # Create a default minimisation protocol.
   protocol = BSS.Protocol.Minimisation()

   # Automatically find an appropriate molecular dynamics package on the
   # host environment, set up the simulation and return a running process.
   # Since the system was generated from a CRD and TOP file, BSS.MD will
   # try to find an AMBER package.
   process = BSS.MD.run(system, protocol)


BioSimSpace.Notebook
====================

This package provides functionality for viewing molecules and creating graphs
when using BioSimSpace interactively within a `Jupyter <http://jupyter.org>`_
notebook.

BioSimSpace.Notebook.View
-------------------------

The ``View`` class provides a wrapper around `NGLView <https://github.com/arose/nglview>`_
to enable visualisation of molecular systems. ``View`` can directly display a ``Sire.System``\ ,
or can attach to a running `Process <../Process>`_ to allow real-time visualisation
of molecular configurations while a process is running.

.. code-block:: python

   import BioSimSpace as BSS

   # Create a molecular system.
   system = BSS.IO.readMolecules(["ala.crd", "ala.top"])

   # Create a view object for the system.
   view = BSS.NoteBook.View(system)

   # View all of the molecules in the system.
   view.molecules()

   # View the first molecule in the system.
   view.molecule(0)

   # View a list of molecules.
   view.molecule([0, 5, 10])

   # Create a default equilibration protocol.
   protocol = BSS.Protocol.Equilibration()

   # Start a simulation process.
   process = BSS.MD.run(system, protocol)

   # Attach a view to the running process.
   view = BSS.Notebook.View(process)

   # View the latest molecular system.
   view.molecules()

BioSimSpace.Notebook.plot
-------------------------

This is a simple wrapper around `Matplotlib <https://matplotlib.org>`_ for
plotting time-series data.

For example:

.. code-block:: python

   # Generate a plot of time vs total energy.
   BSS.Notebook.plot(process.getTime(time_series=True), process.getTotalEnergy(time_series=True),
       xlabel="Time (ns)", ylabel="Energy (kcal/mol)")


Demo
====

This directory contains example scripts that illustrate various ways of using
BioSimSpace. If you have already installed BioSimSpace into a Sire package then
you can run any of the scripts directly, e.g.:

.. code-block:: bash

   # Run an example AMBER workflow.
   $HOME/sire.app/bin/python amber.py

   # Run an example NAMD workflow.
   # You'll need to have installed NAMD on your computer.
   $HOME/sire.app/bin/python namd.py

   # Run a generic minimisation node using AMBER input files.
   $HOME/sire.app/bin/python minimisation.py --steps=1000 --files amber/ala/*

   # Run a generic minimisation node using NAMD input files.
   # You'll need to have installed NAMD on your computer.
   $HOME/sire.app/bin/python minimisation.py --steps=1000 --files namd/ala*/*

Alternatively, e.g. when developing, set the ``PYTHONPATH`` environment so that
the python interpreter can locate BioSimSpace:

.. code-block:: bash

   PYTHONPATH=../python $HOME/sire.app/bin/python amber.py

Notebooks
---------

We also provide several `Jupyter <http://jupyter.org>`_ notebooks that show how
to work with BioSimSpace in an interactive environment. Current notebooks are:


* `interactive_md.ipynb <interactive_md.ipynb>`_ - A short notebook showing how
  to interact with a running molecular simulation process in real-time.
* `minimisation.ipynb <minimisation.ipynb>`_ - An example BioSimSpace node
  showing how to write a fully documented, validated, and portable workflow
  component.
* `amber.ipynb <amber.ipynb>`_ - A detailed notebook showing how to use the
  low-level functionality of BioSimSpace to set up and configure
  `AMBER <http://ambermd.org>`_ simulations.
* `conversion.ipynb <conversion.ipynb>`_ - A simple example node to convert
  between common molecular file formats.
* `input_files.ipynb <input_files.ipynb>`_ - A node to auto-generate input
  files for different molecular simulation protocols using different molecular
  dynamics packages.
* `molecular_setup.ipynb <molecular_setup.ipynb>`_ - A node to parameterise
  and solvate a molecular system ready for simulation.

To run a notebook, simply launch the Jupyter notebook app:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter notebook

This will open a notebook dashboard in your browser, where you can click on
the notebook of your choice. You can also launch a specific notebook directly,
e.g.:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter notebook interactive_md.ipynb


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

To run a notebook, simply launch the Jupyter notebook app:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter notebook

This will open a notebook dashboard in your browser, where you can click on
the notebook of your choice. You can also launch a specific notebook directly,
e.g.:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter notebook interactive_md.ipynb

If you used the regular `setup.py <../python/setup.py>`_ installation, then you
should be good to go. If not, the following describes full installation
instructions for the required dependencies.

To enable plotting functionality in the example notebooks you will also need
to install `Matplotlib <https://matplotlib.org>`_\ :

.. code-block:: bash

   $HOME/sire.app/bin/pip install matplotlib

Molecular visualisation requires `NGLView <https://github.com/arose/nglview>`_\ ,
which can be a bit of a pain to configure correctly.

First install Jupyter:

.. code-block:: bash

   $HOME/sire.app/bin/pip install jupyter

Now install NGLView, overwriting any dependencies that were already installed above:

.. code-block:: bash

   $HOME/sire.app/bin/pip install --no-cache-dir --ignore-installed nglview

We now need to enable the notebook extension:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter-nbextension enable nglview --py --sys-prefix

Versions of NGLView earlier than 1.1.2 have an issues that causes intermittent
problems when reading molecular structures from file. This can be fixed by
running the provided helper script:

.. code-block:: bash

   ./nglview_patch.sh $HOME/sire.app

To find out what version of NGLView you have installed, run:

.. code-block:: bash

   $HOME/sire.app/bin/pip show nglview

In a similar manner we need to install and activate the Jupyter
`fileupload <https://pypi.python.org/pypi/fileupload>`_ widget.

.. code-block:: bash

   $HOME/sire.app/bin/pip install fileupload
   $HOME/sire.app/bin/jupyter-nbextension install fileupload --py --sys-prefix
   $HOME/sire.app/bin/jupyter-nbextension enable fileupload --py --sys-prefix

Finally, you may need to clean any existing Jupyter configuration folders, which
will likely be present if you already have Jupyter installed on your machine, i.e.
installed outside of your ``sire.app``. For simplicity we'll delete these files,
although you may want to back them up:

.. code-block:: bash

   rm -r ~/.jupyter ~/.local/share/jupyter

(In future, we probably want to set different ``JUPYTER_CONFIG_DIR`` and
``JUPYTER_PATH`` environment variables for BioSimSpace, or provide a custom
``jupyter`` executable that sets the variables before launching Jupyter.)

Phew!

You should now be able to launch Jupyter:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter notebook

   or

   PYTHONPATH=../python $HOME/sire.app/bin/jupyter notebook

This will open a tab in your browser, where you can choose a specific demo
notebook, such as ``amber.ipynb``. If you find that you have issues visualising
systems with a large number of molecules, e.g. solvated systems, then you
might need to launch Jupyter as follows:

.. code-block:: bash

   $HOME/sire.app/bin/jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000

Remember to shutdown the python kernel before closing any notebook tabs so that
any BioSimSpace processes that are running can be safely killed and cleaned up.

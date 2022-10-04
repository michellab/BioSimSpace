.. _ref_install:

============
Installation
============

Conda install
=============

The easiest way to install BioSimSpace is using our `conda channel <https://anaconda.org/michellab/repo>`__.
BioSimSpace is built using dependencies from `conda-forge <https://conda-forge.org/>`__,
so please ensure that the channel takes strict priority. We recommend using
`Miniforge <https://github.com/conda-forge/miniforge>`__.

To create a new environment:

.. code-block:: bash

    conda create -n biosimspace -c conda-forge -c michellab biosimspace
    conda activate biosimspace

To install the latest development version you can use:

.. code-block:: bash

    conda create -n biosimspace-dev -c conda-forge -c michellab/label/dev biosimspace
    conda activate biosimspace-dev

When updating the development version it is generally advised to update `Sire <https://github.com/michellab/Sire>`_
at the same time:

.. code-block:: bash

    conda update -c conda-forge -c michellab/label/dev biosimspace sire

If you plan on using BioSimSpace interactively via Jupyter, then you might also
need to enable the required notebook extensions within your Conda environment:

.. code-block:: bash

    jupyter-nbextension enable nglview --py --sys-prefix

Unless you add the required channels to your Conda configuration, then you'll
need to add them when updating, e.g., for the development package:

.. code-block:: bash

    conda update -c conda-forge -c michellab/label/dev biosimspace

If you find that Conda is particularly slow to install or upgrade BioSimSpace,
then we advise using `mamba <https://github.com/TheSnakePit/mamba>`__:

.. code-block:: bash

    conda install -c conda-forge mamba

You can then replace all ``conda`` commands with ``mamba``, e.g.:

.. code-block:: bash

    mamba create -n biosimspace -c conda-forge -c michellab biosimspace

Binary install
==============

Prior to the 2022.1.0 release we provided a self-extracting binary install
of BioSimSpace. Please visit our :ref:`binaries <ref_binaries>` page if you
wish to download any of these.

For developers
==============

The following documents a full installation of BioSimSpace from source. Before
starting, you'll need a working `Git <https://git-scm.com>`__ installation.

BioSimSpace is built on top of the `Sire <https://github.com/michellab/Sire>`__
molecular simulation framework. To download and install Sire, follow the
instructions `here <https://github.com/michellab/Sire#installation>`__, making
sure that BioSimSpace's dependencies are installed into the Sire conda
environment at the point at which Sire is installed.

Next you will need to download BioSimSpace and install it into your Sire
Conda environment.

.. code-block:: bash

   git clone https://github.com/michellab/BioSimSpace
   cd BioSimSpace/python
   python setup.py install

Once finished, you can test the installation by running:

.. code-block:: bash

   python

Then try importing the BioSimSpace package:

.. code-block:: python

   import BioSimSpace as BSS

When developing you may not wish to continually re-install BioSimSpace and its
associated dependencies. To avoid this, you can either make use of ``PYTHONPATH``,
e.g.

.. code-block:: bash

   PYTHONPATH=$HOME/Code/BioSimSpace/python python script.py

or use the ``develop`` argument when running the ``setup.py`` script, i.e.

.. code-block:: bash

   python setup.py develop

You can also skip installation of external dependencies by setting the
environment variable ``BSS_SKIP_DEPENDENCIES``, e.g.

.. code-block:: bash

   BSS_SKIP_DEPENDENCIES=True python setup.py install

Common issues
=============

* If you experience problems with `Matplotlib <https://matplotlib.org>`__ when
  importing BioSimSpace on macOS, e.g.

.. code-block:: bash

   RuntimeError**: Python is not installed as a framework.

simply add the following to ``~/.matplotlib/matplotlibrc``

.. code-block:: bash

   backend: TkAgg

Note that plotting functionality will be disabled if you are using
BioSimSpace on a remote server without X forwarding.

* If you experience problems with `Jupyter <https://jupyter.org>`__ permissions,
  try removing ``$HOME/.jupyter`` or ``$HOME/.local/share/jupyter``

External dependencies
=====================

Several additional packages are required for full access to all of BioSimSpace's
functionality. Please download and install these packages according to their
recommended installation instructions.

* `Amber / AmberTools <http://ambermd.org>`__ -- *Dynamics / Parameterisation*
* `Gromacs <http://www.gromacs.org>`__ -- *Dynamics / Parameterisation / Solvation*
* `Namd <http://www.ks.uiuc.edu/Research/namd>`__ -- *Dynamics*

For `Amber / AmberTools <http://ambermd.org>`__, we also recommend adding
``${AMBERHOME}/bin`` to your ``PATH`` to ensure that its binaries are
visible to third-party libraries, such as
`openff-toolkit <https://github.com/openforcefield/openff-toolkit>`__.

Please visit our :ref:`compatibility <ref_compatibility>` page to see which
versions of the external dependencies BioSimSpace has currently been tested
against.

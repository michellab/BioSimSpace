.. _ref_install:

============
Installation
============

Binary install
==============

The latest self-extracting binary for the development version of BioSimSpace
can be downloaded from the following link:

* `biosimspace_devel_latest_linux.run <https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/ZH4wscDHe59T28yVJtrMH8uqifI_ih0NL5IyqxXQjSo/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_linux.run>`_

One downloaded, the binary can be unpacked as follows:

.. code-block:: bash

   chmod +x biosimspace_devel_latest_linux.run
   ./biosimspace_devel_latest_linux.run

For developers
==============

The following documents a full installation of BioSimSpace from source. Before
starting, you'll need working `Git <https://git-scm.com>`_ and `CMake <https://cmake.org>`_
installations.

1. BioSimSpace is built on top of the `Sire <https://github.com/michellab/Sire>`_
   molecular simulation framework. To download and install Sire:

.. code-block:: bash

   git clone https://github.com/michellab/Sire
   cd Sire
   ./compile_sire.sh

Assuming the default installation path, this will install Sire into ``$HOME/sire.app``.

(Note that the installation is slow and can take in excess of an hour.)

2. Next you will need to download BioSimSpace and install it into your Sire
   application. (The following assumes the default Sire installation path.)

.. code-block:: bash

   git clone https://github.com/michellab/BioSimSpace
   cd BioSimSpace/python
   $HOME/sire.app/bin/python setup,py install

Once finished, you can test the installation by running:

.. code-block:: bash

   $HOME/sire.app/bin/ipython

Then try importing the BioSimSpace package:

.. code-block:: python

   import BioSimSpace as BSS

Common issues
=============

* If you experience problems with `Matplotlib <https://matplotlib.org>`_ when
  importing BioSimSpace on macOS, e.g.

.. code-block:: bash

   RuntimeError**: Python is not installed as a framework.

simply add the following to ``~/.matplotlib/matplotlibrc``

.. code-block:: bash

   backend: TkAgg

Note that plotting functionality will be disabled if you are using
BioSimSpace on a remote server without X forwarding.

* If you experience problems with `Jupyter <https://jupyter.org>`_ permissions,
  try removing ``$HOME/.jupyter`` or ``$HOME/.local/share/jupyter``

External dependencies
=====================

Several additional packages are required for full access to all of BioSimSpace's
functionality. Please download and install these packages according to their
recommended installation instructions.

* `Amber / AmberTools <http://ambermd.org>`_ -- *Dynamics / Parameterisation*
* `Gromacs <http://www.gromacs.org>`_ -- *Dynamics / Parameterisation / Solvation*
* `Namd <http://www.ks.uiuc.edu/Research/namd>`_ -- *Dynamics*

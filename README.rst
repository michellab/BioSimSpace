`BioSimSpace <http://biosimspace.org>`__
========================================

.. image:: https://dev.azure.com/michellab/BioSimSpace/_apis/build/status/michellab.BioSimSpace?branchName=devel
   :target: https://dev.azure.com/michellab/BioSimSpace/_build
   :alt: Build Status

.. image:: https://anaconda.org/michellab/biosimspace/badges/downloads.svg
   :target: https://anaconda.org/michellab/biosimspace
   :alt: Conda Downloads

.. image:: https://img.shields.io/badge/License-GPL%20v2-blue.svg
   :target: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
   :alt: License

About
-----

`BioSimSpace <https://biosimspace.org>`__ is an interoperable Python framework
for biomolecular simulation. With it you can:

* Write robust and portable biomolecular workflow components that work on
  different hardware, with different software packages, and that can be
  run in different ways, e.g. command-line, `Jupyter <https://jupyter.org>`__.
* Interact with molecular-simulation processes in real time.

Documentation
-------------

Full documentation can be found `here <https://biosimspace.org>`__.

Installation
------------

Conda package
^^^^^^^^^^^^^

The easiest way to install BioSimSpace is using our `conda channel <https://anaconda.org/michellab/repo>`__:

.. code-block:: bash

    conda install -c conda-forge -c omnia -c michellab biosimspace

To install the latest development version you can use:

.. code-block:: bash

    conda install -c conda-forge -c omnia -c michellab/label/dev biosimspace

If you plan on using BioSimSpace interactively via Jupyter, then you might also
need to enable the required notebook extensions within your Conda environment:

.. code-block:: bash

    jupyter-nbextension enable fileupload --py --sys-prefix
    jupyter-nbextension enable nglview --py --sys-prefix

Unless you add the required channels to your Conda configuration, then you'll
need to add them when updating, e.g., for the development package:

.. code-block:: bash

    conda update -c conda-forge -c omnia -c michellab/label/dev biosimspace

Note that because of Conda's peculiar scoring metrics you might not end up with
the latest version of BioSimSpace when performing a fresh install or update.
(It tries to minimise various things, such as the number of dependencies
installed, which is difficult when your package depends on many other packages.)
To see what packages are available, run:

.. code-block:: bash

    conda search -c michellab/label/dev biosimspace

You can then install the latest version by explicitly stating the full package
name, e.g.:

.. code-block:: bash

    conda install -c conda-forge -c omnia -c michellab/label/dev biosimspace=2019.1.0=py37h14c3975_85

Using the prebuilt binaries
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The latest self-extracting binary for the development version of BioSimSpace
can be downloaded from one of the following links:

* Linux: `biosimspace_devel_latest_linux.run <https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/ZH4wscDHe59T28yVJtrMH8uqifI_ih0NL5IyqxXQjSo/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_linux.run>`__
* Mac OS X: `biosimspace_devel_latest_osx.run <https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/whcwfvWfndjA4RxupM-4gsVsjcdR0w5I9aP1RJKPruQ/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_osx.run>`__

Once downloaded, the binary can be unpacked as follows:

.. code-block:: bash

   chmod +x biosimspace_devel_latest_linux.run
   ./biosimspace_devel_latest_linux.run

Unless a different installation path was given, BioSimSpace can be found in:
``$HOME/biosimspace.app``. BioSimSpace comes with a bundled with a Python
interpreter, an interactive Python (IPython) shell, and `Jupyter <https://jupyter.org>`__.

For example, to run a BioSimSpace Python script, use:

.. code-block:: bash

   $HOME/biosimspace.app/bin/python script.py

To launch an interactive BioSimSpace session:

.. code-block:: bash

   $HOME/biosimspace.app/bin/ipython

To run a BioSimSpace notebook:

.. code-block:: bash

   $HOME/biosimspace.app/bin/jupyter notebook notebook.ipynb

Installing from source
^^^^^^^^^^^^^^^^^^^^^^

Alternatively, to install BioSimSpace from source:

(Before starting, you'll need a working `Git <https://git-scm.com>`__ installation.)

BioSimSpace is built on top of the `Sire <https://github.com/michellab/Sire>`__
molecular simulation framework. To download and install Sire:

.. code-block:: bash

   git clone https://github.com/michellab/Sire
   cd Sire
   ./compile_sire.sh

Assuming the default installation path, this will install Sire into ``$HOME/sire.app``.

(Note that the installation is slow and can take in excess of an hour.)

Next you will need to download BioSimSpace and install it into your Sire
application. (The following assumes the default Sire installation path.)

.. code-block:: bash

   git clone https://github.com/michellab/BioSimSpace
   cd BioSimSpace/python
   $HOME/sire.app/bin/python setup.py install

Once finished, you can test the installation by running:

.. code-block:: bash

   $HOME/sire.app/bin/ipython

Then try importing the BioSimSpace package:

.. code-block:: python

   import BioSimSpace as BSS

Docker images
-------------

If you don't want to build or install, you can also run BioSimSpace via one of
our docker images. The easy way to run the latest development image of
BioSimSpace is via:

.. code-block:: bash

   docker run -it biosimspace/biosimspace-devel:latest

This will download the latest BioSimSpace development container, and will run
it, giving you a bash prompt inside the container.

OpenMM compatibility
--------------------

Some BioSimSpace functionality requires `OpenMM <http://openmm.org>`__. Although
a bundled version is provided as part of the installation, this may not
be appropriate for your GPU drivers. To automatically detect and install
a suitable version of OpenMM, simply run the following command post-install::

    optimise_openmm

(Note that, depending on your installation method, ``optimise_openmm`` may
be located in ``$HOME/sire.app/bin``.)

Alternatively, to manually install a particular version of OpenMM you can
use a specific Conda label, e.g.::

    conda install -c omnia/label/cuda90 openmm

If you have compiled Sire against a custom OpenMM installation, then you'll
need to set the ``OPENMM_PLUGIN_DIR`` environment variable to point to the
correct plugin location. By default this variable is set to the plugin
directory of the bundled OpenMM package.

Developers
----------

Please follow the `developer's guide <https://biosimspace.org/development.html>`__.

Issues
------

Please report bugs and other issues using the GitHub `issue tracker <https://github.com/michellab/BioSimSpace/issues>`__.
When reporting issues please try to include a minimal code snippet that reproduces
the problem. Additional files can be also be uploaded as an archive, e.g. a zip
file. Please also report the branch on which you are experiencing the issue,
along with the BioSimSpace version number. This can be found by running:

.. code-block:: python

   import BioSimSpace as BSS
   print(BSS.__version__)

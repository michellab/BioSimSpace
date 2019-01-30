
BioSimSpace
===========

.. image:: https://dev.azure.com/michellab/BioSimSpace/_apis/build/status/michellab.BioSimSpace?branchName=devel
   :target: https://dev.azure.com/michellab/BioSimSpace/_build
   :alt: Build Status

Code and resources for the `EPSRC <https://epsrc.ukri.org>`_
`BioSimSpace <https://michellab.github.io/BioSimSpaceWebsite>`_ project.

What is it?
-----------

BioSimSpace is an interoperable Python framework for biomolecular simulation.
With it you can:

* Write robust and portable biomolecular workflow components that work on
  different hardware, with different software packages, and that can be
  run in different ways, e.g. command-line, `Jupyter <https://jupyter.org>`_.
* Interact with running molecular simulation processes in real-time.

Documentation
-------------

Full documentation can be found `here <https://michellab.github.io/BioSimSpaceWebsite>`_.

Installation
------------

The latest self-extracting binary for the development version of BioSimSpace
can be downloaded from the following link:

* `biosimspace_devel_latest_linux.run <https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/ZH4wscDHe59T28yVJtrMH8uqifI_ih0NL5IyqxXQjSo/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_linux.run>`_

One downloaded, the binary can be unpacked as follows:

.. code-block:: bash

   chmod +x biosimspace_devel_latest_linux.run
   ./biosimspace_devel_latest_linux.run

Alternatively, to install BioSimSpace from source:

(Before starting, you'll need working `Git <https://git-scm.com>`_ and
`CMake <https://cmake.org>`_ installations.)

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

Issues
------

Please report bugs and other issues using the GitHub `issue tracker <https://github.com/michellab/BioSimSpace/issues>`_.
When reporting issues please try to include a minimal code snippet that reproduces
the problem. Additional files can be also be uploaded as an archive, e.g. a zip
file. Please also report the branch on which you are experiencing the issue,
along with the BioSimSpace version number. This can be found by running:

.. code-block:: python

   import BioSimSpace as BSS
   print(BSS.__version__)

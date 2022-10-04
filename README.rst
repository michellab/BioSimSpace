`BioSimSpace <http://biosimspace.org>`__
========================================

.. image:: https://github.com/michellab/BioSimSpace/workflows/Build/badge.svg
   :target: https://github.com/michellab/BioSimSpace/actions?query=workflow%3ABuild)
   :alt: Build status

.. image:: https://anaconda.org/michellab/biosimspace/badges/downloads.svg
   :target: https://anaconda.org/michellab/biosimspace
   :alt: Conda Downloads

.. image:: https://img.shields.io/badge/License-GPL%20v2-blue.svg
   :target: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
   :alt: License

.. image:: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225/status.svg
   :target: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225
   :alt: Paper

About
-----

`BioSimSpace <https://biosimspace.org>`__ is an interoperable Python framework
for biomolecular simulation. With it you can:

* Write robust and portable biomolecular workflow components that work on
  different hardware, with different software packages, and that can be
  run in different ways, e.g. command-line, `Jupyter <https://jupyter.org>`__.
* Interact with molecular-simulation processes in real time.

Citation |DOI for Citing BioSimSpace|
=====================================

If you use BioSimSpace in any scientific software, please cite the following paper: ::

    @article{Hedges2019,
      doi = {10.21105/joss.01831},
      url = {https://doi.org/10.21105/joss.01831},
      year = {2019},
      publisher = {The Open Journal},
      volume = {4},
      number = {43},
      pages = {1831},
      author = {Lester Hedges and Antonia Mey and Charles Laughton and Francesco Gervasio and Adrian Mulholland and Christopher Woods and Julien Michel},
      title = {BioSimSpace: An interoperable Python framework for biomolecular simulation},
      journal = {Journal of Open Source Software}
    }

.. |DOI for Citing BioSimSpace| image:: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225/status.svg
   :target: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225

Documentation
-------------

Full documentation can be found `here <https://biosimspace.org>`__.

Installation
------------

Conda package
^^^^^^^^^^^^^

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

Installing from source
^^^^^^^^^^^^^^^^^^^^^^

Alternatively, to install BioSimSpace from source:

(Before starting, you'll need a working `Git <https://git-scm.com>`__ installation.)

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

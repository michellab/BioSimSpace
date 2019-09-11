.. _ref-Metadynamics:

BioSimSpace.Metadynamics
========================

The *Metadynamics* package contains tools automatically configure, run, and
analyse metadynamics simulations.

Metadynamics support requires `PLUMED <https://www.plumed.org>`__. We provide
support for version 2.5 and above. Once installed, make sure that the
``plumed`` binary is visible in your ``PATH`` so that it can be found by
BioSimSpace. You will also need to ensure that the PLUMED kernel is
visible, and likely need to update your ``LD_LIBRARY_PATH`` too, e.g:

.. code-block:: bash

    export PLUMED_KERNEL=/usr/local/lib/libplumedKernel.so
    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

To use PLUMED you will also have to patch an appropriate molecular dynamics
engine, following the instructions
`here <https://www.plumed.org/doc-v2.5/user-doc/html/_installation.html>`__.
At present we only support metadynamics simulations using
`GROMACS <http://www.gromacs.org>`__.

.. automodule:: BioSimSpace.Metadynamics

Packages
--------

.. toctree::
   :maxdepth: 1

   index_Metadynamics_CollectiveVariable

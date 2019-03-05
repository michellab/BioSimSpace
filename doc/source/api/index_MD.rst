.. _ref-MD:

BioSimSpace.MD
==============

The *MD* package provides functionality for automatically configuring and
running molecular dynamics simulations. The function will choose the most
appropriate molecular dynamics driver based on the software available
on the host, the available hardware, the molecular
:class:`System <BioSimSpace._SireWrappers.System>`, and
:class:`Protocol <BioSimSpace.Protocol>`, returning the user a handle to
a :class:`Process <BioSimSpace.Process>` object to run the simulation.

.. automodule:: BioSimSpace.MD

.. toctree::
   :maxdepth: 1

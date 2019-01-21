
BioSimSpace.Protocol
====================

This package provides functionality for defining protocols for molecular
simulation.

The ``_protocol.py`` module provides a ``Protocol`` base class, containing
properties common to all protocol types. Classes for specific protocols are
defined in their own module files, e.g. ``_minimisation.py``.

Protocols are used by modules within the ``BioSimSpace.Process`` package to
define what type of simulation to run. The information contained in the
protocol object is used to create the input and configuration files needed
to run a particular process.

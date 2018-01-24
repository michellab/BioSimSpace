# BioSimSpace.Protocol

This sub-package provides functionality for defining protocols for molecular
simulation.

The `ProtocolType` class provides a simple enumerated type to store the
currently supported protocols. At present we prodvide support for minimisation
(`ProtocolType.MINIMISATION`), equilibration (`ProtocolType.EQUILIBRATION`),
and production (`ProtocolType.PRODUCTION`) simulations.

The `protocol.py` module also provides a `Protocol` base class, containing
properties common to all protocol types. Classes for specific protocols are
defined in their own module files, e.g. `minimisation.py`.

Protocols are used by modules within the `BioSimSpace.Process` sub-package to
define what type of simulation to run. The information contained in the
protocol object is used to create the input and configuration files needed
to run a particular process.

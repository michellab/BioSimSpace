.. _ref_compatibility:

=============
Compatibility
=============

BioSimSpace has been tested against the versions of its
external dependencies listed below. Please let us know if
you encounter problems by raising an issue on our
`GitHub <https://github.com/openbiosim/biosimspace/issues>`__
page. (It may also work with more recent versions and we
will update these lists as they are validated.)

BioSimSpace is built to be compatible with the
`ambertools <https://anaconda.org/conda-forge/ambertools>`__
and `gromacs <https://anaconda.org/conda-forge/gromacs>`__
packages from conda-forge, but this aren't required as
a hard run-time dependency. This is because users will
likely choose to install external versions of the packages
that are optimised for their particular computing environment
and use case.


`AMBER <http://ambermd.org>`__
==============================

* AmberTools22
* AmberTools21
* AmberTools20
* AmberTools19
* AmberTools18

Note that AmberTools >= 20 is required for Open Force Field support.
We currently have support for force fields with CMAP terms, e.g. ff19SB.

`GROMACS <http://www.gromacs.org>`__
====================================

* 2021 series
* 2020 series
* 2019 series
* 2018 series

`NAMD <https://www.ks.uiuc.edu/Research/namd>`__
================================================

* 2.14
* 2.13
* 2.12

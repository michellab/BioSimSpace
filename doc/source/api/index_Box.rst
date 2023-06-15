.. _ref-Box:

BioSimSpace.Box
===============

The *Box* package contains tools for generating parameters for different
simulation boxes. This is particularly useful when needing box magnitudes and
angles for generating :ref:`ref-Solvent` boxes.

To support triclinic boxes that work across the range of molecular simulation
engines that support, we represent the triclinic space in reduced form, using
the approach documented in Appendix A of Chapter 3 from
`"Molecular dynamics of sense and sensibility in processing and analysis of data" <https://research.rug.nl/files/2839528/01_c1.pdf>`__
by Tsjerk A. Wassenaar. Due to the fixed-width format used to represent box
dimensions and angles in the various molecular input files, repeated reading
and writing can lead to oscillation of the box angles on reduction due to
rounding precision errors. To account for this, we add a small bias so that
we always round in a consistent direction.

.. automodule:: BioSimSpace.Box

.. toctree::
   :maxdepth: 1

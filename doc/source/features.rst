========
Features
========

:mod:`BioSimSpace` is an interoperable Python framework for biomolecular
simulation.

* :data:`Read and write <BioSimSpace.IO>` molecules from a number of molecular
  file formats. This includes AMBER, CHARMM, and GROMACS files.
* :data:`Parameterise <BioSimSpace.Parameters>` molecules using force fields from
  `AmberTools <https://ambermd.org/AmberTools.php>`__ and the
  `Open Force Field Initiative <https://openforcefield.org>`__.
* :data:`Solvate <BioSimSpace.Solvent>` molecules using a range of standard water
  models. We support both orthorhombic and triclinic simulation boxes.
* Configure :data:`protocols <BioSimSpace.Protocol>` for a range of biomolecular
  simulations and create :data:`processes <BioSimSpace.Process>` to run them using
  our supported molecular dynamics (MD) engines (
  `AMBER <https://ambermd.org>`__, `GROMACS <https://www.gromacs.org>`__,
  `OpenMM <https://openmm.org>`__, `NAMD <http://www.ks.uiuc.edu/Research/namd>`__, SOMD).
* Setup and run :data:`free energy perturbation <BioSimSpace.FreeEnergy>` simulations
  with SOMD and GROMACS. (See our :doc:`hydration free energy <tutorials/hydration_freenrg>`
  tutorial for a complete example.)
* Configure and run :data:`metadynamics <BioSimSpace.Metadynamics>` simulations
  for a range of supported collective variables. (See our :doc:`tutorial <tutorials/metadynamics>`
  for a complete example.)
* Write interoperable workflow components or :doc:`nodes <guides/nodes>` that
  work with any supported MD engine and can be run from the command-line,
  within a `Jupyter <https://jupyter.org>`__ notebook, or exported to a
  `Common Workflow Language <https://www.commonwl.org>`__ Tool wrapper.
* Interact with molecular simulation process within interactive Python
  environments. This includes: starting and stopping simulations, monitoring
  their progress, getting time series data of thermodynamic state information,
  plotting data, analysing trajectories, and visualising molecules.

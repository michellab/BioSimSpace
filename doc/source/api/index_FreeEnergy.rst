.. _ref-FreeEnergy:

BioSimSpace.FreeEnergy
======================

The *FreeEnergy* package contains tools to configure, run, and analyse
*relative* free energy simulations.

Free-energy perturbation simulations require a
:class:`System <BioSimSpace._SireWrappers.System>` containing a *merged*
molecule that can be *perturbed* between two molecular end states by use
of an *alchemical potential*. To create merged molecules, please use the
:ref:`ref-Align` package.

Relative free-energy calculations require the simulation of two perturbations,
typically referred to as *legs*. A potential of mean force (PMF) is computed
for each leg, which can then be used to computed the relative free-energy
difference. For generality and flexibility, BioSimSpace decouples the two legs,
allowing the use of difference molecular simulation engines,
:class:`protocols <BioSimSpace.Protocol.FreeEnergy>`, and for legs to be
re-used in different calculations.

Simulations are typically used to compute solvation (currently hydration only)
or binding free-energies. In the examples that follow, ``merged`` refers to a
perturbable molecule created by merging two ligands, ``ligA`` and ``ligB``,
``merged_sol`` refers to the same perturbable molecule in solvent, and
``complex_sol`` is a solvated protein-ligand complex containing the same
perturbable molecule. We assume that each molecule/system has been
appropriately minimised and equlibrated.

To setup, run, and analyse a binding free-energy calculation:

.. code-block:: python

   import BioSimSpace as BSS

   ...

   # Create two a protocol for the two legs of a binding free-energy simulation.
   # Use more lambda windows for the "bound" leg.
   protocol_bound = BSS.Protocol.FreeEnergy(num_lam=20)
   protocol_free  = BSS.Protocol.FreeEnergy(num_lam=12)

   # Setup the perturbations for each leg, using the SOMD engine. This will
   # create all of the input files and simulation processes that are required.
   fep_bound = BSS.FreeEnergy.Relative(
       complex_sol,
       protocol_bound,
       engine="somd",
       work_dir="ligA_ligB/bound"
   )
   fep_free  = BSS.FreeEnergy.Relative(
       merged_sol,
       protocol_bound,
       engine="somd",
       work_dir="ligA_ligB/free"
   )

   # Run all simulations for each leg. Note that the lambda windows are run
   # sequentially, so this is a sub-optimal way of executing the simulation
   # if you have access to HPC resources.

   # Bound leg.
   fep_bound.run()
   fep_bound.wait()

   # Free leg.
   fep_free.run()
   fep_free.wait()

   # Analyse the simulation data from each leg, returning the PMF and overlap
   # matrix.
   pmf_bound, overlap_bound = fep_bound.analyse()
   pmf_free,  overlap_free  = fep_free.analyse()

   # Compute the relative free-energy difference.
   free_nrg_binding = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)

Similarly, for a solvation free-energy calculation:

.. code-block:: python

   # Here we are assuming that we are using the same ligands, so will re-use
   # the free leg from the previous example.

   # Setup the perturbation for the vacuum leg using a default protocol.
   fep_vacuum = BSS.FreeEnergy.Relative(
       merged.toSystem(),
       engine="somd",
       work_dir="ligA_ligB/vacuum"
   )

   # Run the simulations for the perturbation.
   fep_vacuum.run()
   fep_vacuum.wait()

   # Analyse the simulation data.
   pmf_vacuum, overlap_vacuum = fep_vacuum.analyse()

   # Compute the relative free-energy difference.
   free_nrg_solvation = BSS.FreeEnergy.Relative.difference(pmf_free, pmf_vacuum)

Since it is usually preferable to run simulations intensive simulation such as
these on external HPC resources, the ``BioSimSpace.FreeEnergy`` package also
provides support for only creating the input files that are needed by passing
the ``setup_only=True`` argument. This saves the overhead of creating
:class:`Process <BioSimSpace.Process>` objects. The input files can then be
copied to a remote server, with the indivual simulations curated in a job
submission script. (We don't yet provide support for configuring and writing
submission scripts for you.)

To just setup the vacuum leg input files:

.. code-block:: python

   # Setup the input for the vacuum leg. No processes are created so the .run()
   # method won't do anything.
   fep_vacuum = BSS.FreeEnergy.Relative(
       merged.toSystem(),
       engine="somd",
       work_dir="ligA_ligB/vacuum",
       setup_only=True
   )


It is also possible to analyse existing simulation output directly by passing
the path to a working directory to :class:`FreeEnergy.Relative.analyse <BioSimSpace.FreeEnergy.Relative.analyse>`:

.. code-block:: python

   pmf_vacuum, overlap_vacuum = BSS.FreeEnergy.Relative.analyse("ligA_ligB/vacuum")

.. automodule:: BioSimSpace.FreeEnergy

.. toctree::
   :maxdepth: 1

As well as the :class:`protocol <BioSimSpace.Protocol.FreeEnergy>` used for production
simulations, it is also possible to use
:class:`FreeEnergy.Relative <BioSimSpace.FreeEnergy.Relative>` to setup and run simulations
for minimising or equilibrating structures for each lambda window. See the
:class:`FreeEnergyMinimisation <BioSimSpace.Protocol.FreeEnergyMinimisation>` and
:class:`FreeEnergyEquilibration <BioSimSpace.Protocol.FreeEnergyEquilibration>`
protocols for details. At present, these protocols are only supported when not
using :class:`SOMD <BioSimSpace.Process.Somd>` as the simulation engine.

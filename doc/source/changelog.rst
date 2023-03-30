Changelog
=========

OpenBioSim
----------

:mod:`BioSimSpace` is supported by `OpenBioSim <https://openbiosim.org>`__, a community interest
company supporting open-source development of fostering academic/industrial collaboration
within the biomolecular simulation community. Our software is hosted via the `OpenBioSim`
`GitHub <https://github.com/OpenBioSim/biosimspace>`__ organisation.

`2023.2.0 <https://github.com/openbiosim/biosimspace/compare/2023.1.2...2023.2.0>`_ - Mar 30 2023
-------------------------------------------------------------------------------------------------

* Make sure that system properties are preserved when creating a new Sire system.
* Fixed an issue with the OpenMM minimisation protocol that meant that the number of steps was ignored (`#12 <https://github.com/OpenBioSim/biosimspace/pull/12>`__).
* Use native Sire PDB downloading functionality to remove ``pypdb`` dependency.
* Fixed an issue with SMILES characters in molecule names causing issues for ``gmx grompp`` (`#14 <https://github.com/OpenBioSim/biosimspace/pull/14>`__).
* Increase default SOMD cut-off since it uses reaction field (`#15 <https://github.com/OpenBioSim/biosimspace/pull/15>`__).
* No longer downcast molecules to single residues and atoms when searching (`#19 <https://github.com/OpenBioSim/biosimspace/pull/19>`__).
* Remove velocities when combining molecules if the property isn't present for all molecules (`#21 <https://github.com/OpenBioSim/biosimspace/pull/21>`__).
* Set default-valued properties when merging molecules to avoid issues with zero values when units are stripped (`#24 <https://github.com/OpenBioSim/biosimspace/pull/24>`__).
* Remove ``watchdog`` to avoid non-deterministic parsing of AMBER output (`#27 <https://github.com/OpenBioSim/biosimspace/pull/27>`__).
* Improved handling of disulphide bonds in multi-chain PDBs sharing the same residue numbers (`#28 <https://github.com/OpenBioSim/biosimspace/pull/28>`__).
* Allow keyword arguments to be passed through to ``lomap`` in :func:`generateNetwork <BioSimSpace.Align.generateNetwork>` (`#29 <https://github.com/OpenBioSim/biosimspace/pull/29>`__).
* Add mixin classes to allow position restraints to be used with a wider range of protocols (`@xiki-tempula <https://github.com/xiki-tempula>`_) and alchemical simulations for non-production protocols (`@msuruzhon <https://github.com/msuruzhon>`_). Switch to using ``gmx energy`` to parse GROMACS energy records (`@xiki-tempula <https://github.com/xiki-tempula>`_) (`#30 <https://github.com/OpenBioSim/biosimspace/pull/30>`__).
* Switch to using native RDKit conversion throughout to avoid conversion via an intermediate file format.
* Expose Sire to OpenMM conversion functionality in :mod:`BioSimSpace.Convert <BioSimSpace.Convert>`.
* Added Python 3.10 support and now build Python 3.10 packages. This is now the default version of Python for BioSimSpace, and the version we recommend for new workflows. Note that we will drop automatic building of Python 3.8 packages later this year (likely Q3 or Q4). This will be timed to co-incide with when we add Python 3.11 support, and when (we anticipate) conda-forge will drop Python 3.8. Our aim is to only build packages for a maximum of 3 Python versions at a time.

`2023.1.2 <https://github.com/openbiosim/biosimspace/compare/2023.1.1...2023.1.2>`_ - Feb 24 2023
-------------------------------------------------------------------------------------------------

* Refactor code to use a unified :class:`WorkDir <BioSimSpace._Utils.WorkDir>` class to simplify the creation of working directories (`#2 <https://github.com/OpenBioSim/biosimspace/pull/2>`__).
* Added :meth:`isSame <BioSimSpace._SireWrappers.System.isSame>` method to compare systems using a sub-set of system and molecular properties. This improves our file caching support, allowing a user to exclude properties when comparing cached systems prior to write, e.g. ignoring coordinates and velocities, if those are the only things that differ between the systems `(#3 <https://github.com/OpenBioSim/biosimspace/pull/3>`__).
* Added the initial version of :mod:`BioSimSpace.Convert <BioSimSpace.Convert>`, which provides support for converting between native `BioSimSpace`, `Sire <http://sire.openbiosim.org>`__, and `RDKit <https://www.rdkit.org>`__ objects (`#9 <https://github.com/OpenBioSim/biosimspace/pull/9>`__).
* Fixed several formatting issues with the website documentation.

`2023.1.1 <https://github.com/openbiosim/biosimspace/compare/2023.1.0...2023.1.1>`_ - Feb 07 2023
-------------------------------------------------------------------------------------------------

* Minor fixes to website documentation.
* Fixed issues with API documentation introduced by `pydocstringformatter <https://pypi.org/project/pydocstringformatter>`__.
* Fixed globbing of GROMACS trajectory files.

`2023.1.0 <https://github.com/openbiosim/biosimspace/compare/2022.3.0...2023.1.0>`_ - Feb 03 2023
-------------------------------------------------------------------------------------------------

* Wrapped the new `sire.load <https://sire.openbiosim.org/api/sire.html#sire.load>`__ function to allow loading of URLs.
* Add basic file caching support to avoid re-writing files for the same molecular system.
* Added :data:`BioSimSpace._Config` sub-package to simplify the generation of configuration files for molecular dynamics engines and improve flexiblity. (Adapted from code written by `@msuruzhon <https://github.com/msuruzhon>`_.)
* Deprecated ``BioSimSpace.IO.glob`` since globbing is now performed automatically.
* Autoformat entire codebase using `black <https://github.com/psf/black>`__.
* Fix issues following Sire 2023 API updates.
* Update documentation for new OpenBioSim website.

Michellab
---------

Prior to January 2023, :mod:`BioSimSpace` was hosted within the `michellab <https://github.com/michellab/BioSimSpace>`__
GitHub organisation. The following releases were made during that time.

`2022.3.0 <https://github.com/openbiosim/biosimspace/compare/2022.2.1...2022.3.0>`_ - Sep 28 2022 (Pre-release)
---------------------------------------------------------------------------------------------------------------

* Improved NAMD restraint implementation for consistency with other engines.
* Make sure we wait for ``trjconv`` to finish when calling as a sub-process.
* Added wrapper for ``Sire.Units.GeneralUnit``.
* Improved interoperability of ``BioSimSpace.Trajectory`` sub-package.
* Added ``BioSimSpace.Sandpit`` for experimental features from external collaborators.
* Added functionality to check for molecules in a ``BioSimSpace.System``.
* Added functionality to extract atoms and residues by absolute index.
* Allow continuation for GROMACS equilibration simulations. (`@kexul <https://github.com/kexul>`_)
* Update BioSimSpace to work with the new Sire 2023.0.0 Python API.

`2022.2.1 <https://github.com/openbiosim/biosimspace/compare/2022.2.0...2022.2.1>`_ - Mar 30 2022
-------------------------------------------------------------------------------------------------

* Fix performance issues when ensuring unique molecule numbering when adding molecules to ``BioSimSpace._SireWrappers.System`` and ``BioSimSpace._SireWrappers.Molecules`` objects.
* Fix extraction of box vector magnitudes for triclinic boxes.

`2022.2.0 <https://github.com/openbiosim/biosimspace/compare/2022.1.0...2022.2.0>`_ - Mar 24 2022
-------------------------------------------------------------------------------------------------

* Use fast C++ wrappers for updating coordinates and velocities during SOMD simulations.
* Fix import issues caused by change in module layout for conda-forge OpenMM package.
* Don't check for structural ions when parameterising with GAFF/GAFF2.
* Fix errors in funnel correction calculation.
* Switch to using conda-forge lomap2 package, removing need to vendor lomap code.
* Use py3Dmol to visualise maximum common substructure mappings.
* Rename ``.magnitude()`` method on ``BioSimSpace.Type`` objects to ``.value()`` to avoid confusion.
* Handle trjconv frame extraction failures within ``BioSimSpace.Process.Gromacs.getSystem()``.
* Catch and handle possible GSL error during singular valued decomposition routine used for molecular alignment.

`2022.1.0 <https://github.com/openbiosim/biosimspace/compare/2020.1.0...2022.1.0>`_ - Jan 26 2022
-------------------------------------------------------------------------------------------------

* Added basic support for cleaning PDB files with `pdb4amber <https://github.com/Amber-MD/pdb4amber>`_ prior to read.
* Added basic support for exporting BioSimSpace Nodes as Common Workflow Language wrappers.
* Added support for parameterising molecules using OpenForceField.
* Added support for using SMILES strings for input to parameterisation functions.
* Added support for funnel metadynamics simulations (`@dlukauskis <https://github.com/dlukauskis>`_).
* Added support for steered molecular dynamics simulations (`@AdeleLip <https://github.com/AdeleLip>`_).
* Added support for generating perturbation networks using LOMAP (`@JenkeScheen <https://github.com/JenkeScheen>`_).
* Fixed bug affecting certain improper/dihedral terms in SOMD perturbation file writer.
* Numerous performance improvements, particularly involving the manipulation and
  combination of molecular systems.
* Native Python pickling support for wrapped Sire types (`@chryswoods <https://github.com/chryswoods>`_).
* Numerous free-energy perturbation pipeline fixes and improvements. Thanks to `@kexul <https://github.com/kexul>`_ and `@msuruzhon <https://github.com/msuruzhon>`_ for their help testing and debugging.
* Switch continuous integration to GitHub actions using conda-forge compliant build and upload to Anaconda cloud.

`2020.1.0 <https://github.com/openbiosim/biosimspace/compare/2019.3.0...2020.1.0>`_ - July 28 2020
--------------------------------------------------------------------------------------------------

* Added logo to website and update theme (`@ppxasjsm <https://github.com/ppxasjsm>`_).
* Make sure potential terms are sorted when writing to SOMD perturbation files (`@ptosco <https://github.com/ptosco>`_).
* Switch to using ipywidgets.FileUpload to eliminate non-conda dependencies.
* Added support for single-leg free energy simulations.
* Created a KCOMBU mirror to avoid network issues during install.
* Allow AMBER simulations when system wasn't loaded from file.
* Handle GROMACS simulations with non-periodic boxes.
* Run vacuum simulations on a single thread when using GROMACS to avoid domain decomposition.
* Make sure BioSimSpace is always built against the latest version of Sire during conda build.

`2019.3.0 <https://github.com/openbiosim/biosimspace/compare/2019.2.0...2019.3.0>`_ - Nov 22 2019
-------------------------------------------------------------------------------------------------

* Make FKCOMBU download during conda build resilient to server downtime.
* Added support for xtc trajectory files and custom protocols with GROMACS.
* Fixed numerous typos in Sphinx documentation.
* Added Journal of Open Source Software paper.

`2019.2.0 <https://github.com/openbiosim/biosimspace/compare/2019.1.0...2019.2.0>`_ - Sep 11 2019
-------------------------------------------------------------------------------------------------

* Switched to using `RDKit <https://www.rdkit.org/>`_ for maximum common substructure (MCS) mappings.
* Handle perturbable molecules for non free-energy protocols with SOMD and GROMACS.
* Added basic metadynamics functionality with support for distance and torsion collective variables.
* Added support for inferring formal charge of molecules.
* Numerous MCS mapping fixes and improvements. Thanks to `@maxkuhn <https://github.com/maxkuhn>`_, `@dlukauskis <https://github.com/dlukauskis>`_, and `@ptosco <https://github.com/ptosco>`_ for help testing and debugging.
* Added Dockerfile to build thirdparty packages required by the BioSimSpace notebook server.
* Exposed Sire search functionality.
* Added thin-wrappers for several additional Sire objects, e.g. Residue, Atom, and Molecules container.
* Performance improvements for searching, indexing, and extracting objects from molecular containers, e.g. System, Molecule.

`2019.1.0 <https://github.com/openbiosim/biosimspace/compare/2018.1.1...2019.1.0>`_ - May 02 2019
-------------------------------------------------------------------------------------------------

* Added support for parameterising proteins and ligands.
* Added support for solvating molecular systems.
* Molecular dynamics drivers updated to support SOMD and GROMACS.
* Support free energy perturbation simulations with SOMD and GROMACS.
* Added Azure Pipeline to automatically build, test, document, and deploy BioSimSpace.
* Created automatic Conda package pipeline.

`2018.1.1 <https://github.com/openbiosim/biosimspace/compare/2018.1.0...2018.1.1>`_ - May 02 2018
-------------------------------------------------------------------------------------------------

* Fixed conda NetCDF issue on macOS. Yay for managing `python environments <https://xkcd.com/1987>`_\ !
* Install conda `ambertools <https://anaconda.org/AmberMD/ambertools>`_ during `setup <python/setup.py>`_.
* Search for bundled version of ``sander`` when running `AMBER <http://ambermd.org>`_ simulation processes.
* Pass executable found by ``BioSimSpace.MD`` to ``BioSimSpace.Process`` constructor.
* Fixed error in RMSD calculation within ``BioSimSpace.Trajectory`` class.
* Improved example scripts and notebooks.

2018.1.0 - May 01 2018
----------------------

* Initial public release of BioSimSpace.

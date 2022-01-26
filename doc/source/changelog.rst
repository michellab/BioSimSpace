Changelog
=========

`2022.1.0 <https://github.com/michellab/BioSimSpace/compare/2020.1.0...2022.1.0>`_ - Jan 26 2022
------------------------------------------------------------------------------------------------

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

`2020.1.0 <https://github.com/michellab/BioSimSpace/compare/2019.3.0...2020.1.0>`_ - July 28 2020
-------------------------------------------------------------------------------------------------

* Added logo to website and update theme (`@ppxasjsm <https://github.com/ppxasjsm>`_).
* Make sure potential terms are sorted when writing to SOMD perturbation files (`@ptosco <https://github.com/ptosco>`_).
* Switch to using ipywidgets.FileUpload to eliminate non-conda dependencies.
* Added support for single-leg free energy simulations.
* Created a KCOMBU mirror to avoid network issues during install.
* Allow AMBER simulations when system wasn't loaded from file.
* Handle GROMACS simulations with non-periodic boxes.
* Run vacuum simulations on a single thread when using GROMACS to avoid domain decomposition.
* Make sure BioSimSpace is always built against the latest version of Sire during conda build.

`2019.3.0 <https://github.com/michellab/BioSimSpace/compare/2019.2.0...2019.3.0>`_ - Nov 22 2019
------------------------------------------------------------------------------------------------

* Make FKCOMBU download during conda build resilient to server downtime.
* Added support for xtc trajectory files and custom protocols with GROMACS.
* Fixed numerous typos in Sphinx documentation.
* Added Journal of Open Source Software paper.

`2019.2.0 <https://github.com/michellab/BioSimSpace/compare/2019.1.0...2019.2.0>`_ - Sep 11 2019
------------------------------------------------------------------------------------------------

* Switched to using `RDKit <https://www.rdkit.org/>`_ for maximum common substructure (MCS) mappings.
* Handle perturbable molecules for non free-energy protocols with SOMD and GROMACS.
* Added basic metadynamics functionality with support for distance and torsion collective variables.
* Added support for inferring formal charge of molecules.
* Numerous MCS mapping fixes and improvements. Thanks to `@maxkuhn <https://github.com/maxkuhn>`_, `@dlukauskis <https://github.com/dlukauskis>`_, and `@ptosco <https://github.com/ptosco>`_ for help testing and debugging.
* Added Dockerfile to build thirdparty packages required by the BioSimSpace notebook server.
* Exposed Sire search functionality.
* Added thin-wrappers for several additional Sire objects, e.g. Residue, Atom, and Molecules container.
* Performance improvements for searching, indexing, and extracting objects from molecular containers, e.g. System, Molecule.

`2019.1.0 <https://github.com/michellab/BioSimSpace/compare/2018.1.1...2019.1.0>`_ - May 02 2019
------------------------------------------------------------------------------------------------

* Added support for parameterising proteins and ligands.
* Added support for solvating molecular systems.
* Molecular dynamics drivers updated to support SOMD and GROMACS.
* Support free energy perturbation simulations with SOMD and GROMACS.
* Added Azure Pipeline to automatically build, test, document, and deploy BioSimSpace.
* Created automatic Conda package pipeline.

`2018.1.1 <https://github.com/michellab/BioSimSpace/compare/2018.1.0...2018.1.1>`_ - May 02 2018
------------------------------------------------------------------------------------------------

* Fixed conda NetCDF issue on macOS. Yay for managing `python environments <https://xkcd.com/1987>`_\ !
* Install conda `ambertools <https://anaconda.org/AmberMD/ambertools>`_ during `setup <python/setup.py>`_.
* Search for bundled version of ``sander`` when running `AMBER <http://ambermd.org>`_ simulation processes.
* Pass executable found by `BioSimSpace.MD <python/BioSimSpace/MD>`_ to `BioSimSpace.Process <python/BioSimSpace/Process>`_ constructor.
* Fixed error in RMSD calculation within `BioSimSpace.Trajectory <python/BioSimSpace/Trajectory>`_ class.
* Improved example scripts and notebooks.

2018.1.0 - May 01 2018
----------------------

* Initial public release of BioSimSpace.

Changelog
=========

`2019.2.0 <https://github.com/michellab/BioSimSpace/compare/2019.2.0...2019.1.0>`_ - Sep 11 2019
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

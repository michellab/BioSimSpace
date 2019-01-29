
Changelog
=========

`2018.1.1 <https://github.com/michellab/BioSimSpace/compare/2018.1.0...2018.1.1>`_ - May 02 2018
----------------------------------------------------------------------------------------------------


* Fixed conda NetCDF issue on macOS. Yay for managing `python environments <https://xkcd.com/1987>`_\ !
* Install conda `ambertools <https://anaconda.org/AmberMD/ambertools>`_ during `setup <python/setup.py>`_.
* Search for bundled version of ``sander`` when running `AMBER <http://ambermd.org>`_ simulation processes.
* Pass executable found by `BioSimSpace.MD <python/BioSimSpace/MD>`_ to `BioSimSpace.Process <python/BioSimSpace/Process>`_ constructor.
* Fixed error in RMSD calculation within `BioSimSpace.Trajectory <python/BioSimSpace/Trajectory>`_ class.
* Improved example scripts and notebooks.

2018.1.0 - May 01 2018
----------------------


* Initial public release of BioSimSpace.

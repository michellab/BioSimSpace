=================
Quick Start Guide
=================

Import :mod:`BioiSmSpace` using

>>> import BioSimSpace as BSS

Load a molecular system from a URL, via :func:`BioSimSpace.IO.readMolecules`.

>>> url = BSS.tutorialUrl()
>>> system = BSS.IO.readMolecules([f"{url}/ala.top", f"{url}/ala.crd"])

.. note ::

   :data:`BioSimSpace.tutorialUrl()` expands to the base URL that contains
   all tutorial files.

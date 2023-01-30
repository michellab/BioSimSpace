=================
Quick Start Guide
=================

Import :mod:`BioiSmSpace` using

>>> import BioSimSpace as BSS

Load a molecular system from a URL, via :func:`BioSimSpace.IO.readMolecules`:

>>> url = BSS.tutorialUrl()
>>> system = BSS.IO.readMolecules([f"{url}/ala.top", f"{url}/ala.crd"])

.. note ::

   :data:`BioSimSpace.tutorialUrl()` expands to the base URL that contains
   all tutorial files.

Write the system back to a bunch of different molecular formats:

>>> BSS.IO.saveMolecules("ala", system, ["gro87", "grotop", "pdb", "mol2"])
['/home/lester/code/openbiosim/biosimspace/ala.gro',
 '/home/lester/code/openbiosim/biosimspace/ala.top',
 '/home/lester/code/openbiosim/biosimspace/ala.pdb',
 '/home/lester/code/openbiosim/biosimspace/ala.mol2']

.. note ::

   The list of options specify the file formats, not the extensions.

Check what file formats are supported:

>>> BSS.IO.fileFormats()
['DCD', 'Gro87', 'GroTop', 'MOL2', 'PDB', 'PRM7', 'PSF', 'RST', 'RST7', 'SDF']

Find out more information about a specific format:

>>> BSS.IO.formatInfo("rst7")
'Amber coordinate/velocity text (ascii) restart files supported from Amber 7 upwards.'

Create a toluene molecule parameterised with the General AMBER Force Field (GAFF):

>>> molecule = BSS.Parameters.gaff("Cc1ccccc1").getMolecule()

Check what force fields are supported:

>>> BSS.Parameters.forceFields()
['ff03',
 'ff99',
 'ff99SB',
 'ff99SBildn',
 'ff14SB',
 'gaff',
 'gaff2',
 'openff_unconstrained-1.0.0-RC2',
 'openff_unconstrained-1.0.1',
 'openff_unconstrained-1.2.1',
 'openff_unconstrained-1.3.1-alpha.1',
 'openff_unconstrained-1.2.0',
 'openff_unconstrained-2.0.0-rc.1',
 'openff_unconstrained-1.1.0',
 'openff_unconstrained-1.1.1',
 'openff_unconstrained-1.3.1',
 'openff_unconstrained-1.0.0-RC1',
 'openff_unconstrained-2.0.0-rc.2',
 'openff_unconstrained-2.0.0',
 'openff_unconstrained-1.3.0',
 'openff_unconstrained-1.0.0']

.. note ::

   The list of supported force fields may differ depending on what version
   of the `openforcefields` package is installed.

Solvate the molecule in a truncated octahedral box of TIP3P water:

>>> box, angles = BSS.Box.truncatedOctahedron(35 * BSS.Units.Length.angstrom)
>>> solvated = BSS.Solvent.tip3p(molecule=molecule, box=box, angles=angles)

.. note ::

   By default the solvated system will be neutralised.

Run 5000 steps of minimisation using any available molecular dynamics engine
and get back the final structure:

>>> protocol = BSS.Protocol.Minimisation(steps=5000)
>>> process = BSS.MD.run(solvated, protocol)
>>> minimised_system = process.getSystem(block=True)

.. note ::

   BioSimSpace processes can also be run interactively. By passing :data:`block=True`
   ensures that we wait for it to finish before returning the minimised system.

# coding: utf-8
# # Parameterise, Solvate and minimise molecule for molecular simulations
# 
# Authors: Antonia Mey, Julien Michel, Lester Hedges

import BioSimSpace as BSS

# Setting up the input and output bits of the node
node = BSS.Gateway.Node(
    "A node to parameterise, solvate and minimise a molecule (protein or small molecule) ready for molecular simulation.")
node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edinburgh")
node.addAuthor(name="Antonia Mey", email="antonia.mey@ed.ac.uk", affiliation="University of Edinburgh")
node.setLicense("GPLv3")

node.addInput("molecule", BSS.Gateway.File(help="A molecular input file, e.g. a PDB file."))
node.addInput("forcefield", BSS.Gateway.String(help="The name of the force field to use for parameterisation.",
                                               allowed=BSS.Parameters.forceFields(), default="ff14SB"))
node.addInput("water", BSS.Gateway.String(help="The name of the water model to use for solvation.",
                                          allowed=BSS.Solvent.waterModels(), default="tip3p"))
node.addInput("shell", BSS.Gateway.Length(help="The extent of the water shell along each axis around the solute.",
                                            unit="angstrom"))
node.addInput("ion_conc", BSS.Gateway.Float(help="The ionic concentration in mol/litre.",
                                            minimum=0, maximum=1, default=0))
node.addInput("minimisation_steps",
              BSS.Gateway.Integer(help="The number of minimisation steps.", minimum=0, maximum=1000000, default=10000))
node.addInput("verbose", BSS.Gateway.Boolean(help="Be loud and noisy", default=False))

# parametrisation output control
node.addInput("parametrisation_base", BSS.Gateway.String(help="The root name of the output files."))
node.addOutput("parametrisation_out",
               BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER format."))

# minimisation output control
node.addInput("minimisation_base", BSS.Gateway.String(help="The base name for the minimised system files"))
node.addOutput("minimisation_out",
               BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER format."))
node.showControls()

isVerbose = node.getInput("verbose")

# No we actually do the work with the node
if isVerbose:
    print("Reading in molecule from file %s" % node.getInput("molecule"))
system = BSS.IO.readMolecules(node.getInput("molecule"))
molecule = system[0]

if isVerbose:
    print("Running a parametrisation with forcefiled: %s" % node.getInput("forcefield"))
molecule = BSS.Parameters.parameterise(molecule, node.getInput("forcefield")).getMolecule()

if isVerbose:
    print("Saving parameterised molecule....")
node.setOutput("parametrisation_out",
               BSS.IO.saveMolecules(node.getInput("parametrisation_base"), molecule, ["prm7", "rst7"]))

if isVerbose:
    print("Solvating system in %s water in a box of size %s with an ion concentrations of %f." % (
    node.getInput("water"), str(node.getInput("shell")), node.getInput('ion_conc')))
system_solvated = BSS.Solvent.solvate(node.getInput("water"), molecule=molecule, shell=node.getInput("shell"),
                                      ion_conc=node.getInput("ion_conc"))

if isVerbose:
    print("Running a minimsation with %d steps." % node.getInput("minimisation_steps"))
protocol = BSS.Protocol.Minimisation(steps=node.getInput("minimisation_steps"))
process = BSS.MD.run(system_solvated, protocol)

node.setOutput("minimisation_out",
               BSS.IO.saveMolecules(node.getInput("minimisation_base"), process.getSystem(block=True),
                                    ["prm7", "rst7"]))

node.validate()

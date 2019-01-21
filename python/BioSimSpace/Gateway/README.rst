
BioSimSpace.Gateway
===================

This package provides functionality to act as a bridge between BioSimSpace
and the outside world. It allows the user to define robust and portable workflow
components (nodes) that can be run from a variety of environments, e.g. within
`Jupyter <http://jupyter.org>`_\ , or from the command-line.

Requirements
------------

Workflow components require input from the user and generate output.
``BioSimSpace.Gateway`` provides a collection of ``Requirement`` types that allow
the user to document the required inputs and outputs, and to specify default
values and constraints. These are:


* ``BioSimSpace.Gateway.Boolean``
* ``BioSimSpace.Gateway.Integer``
* ``BioSimSpace.Gateway.Float``
* ``BioSimSpace.Gateway.String``
* ``BioSimSpace.Gateway.File``
* ``BioSimSpace.Gateway.FileSet``
* ``BioSimSpace.Gateway.Length``
* ``BioSimSpace.Gateway.Area``
* ``BioSimSpace.Gateway.Volume``
* ``BioSimSpace.Gateway.Charge``
* ``BioSimSpace.Gateway.Energy``
* ``BioSimSpace.Gateway.Pressure``
* ``BioSimSpace.Gateway.Temperature``
* ``BioSimSpace.Gateway.Time``

As an example:

.. code-block:: python

   import BioSimSpace as BSS

   # Create an integer requirement with an allowed range of values.
   num_steps = BSS.Gateway.Integer(help="An integer between 1 and 10.", minimum=1, maximum=10)

   # Create a string requirement that has allowed values.
   animal = BSS.Gateway.String(help="Type of animal.", allowed=["cat", "dog", "fish"], default="dog")

   # Create a requirement for a set of input files.
   files = BSS.Gateway.FileSet(help="A set of input files.")

BioSimSpace.Gateway.Node
------------------------

A node is used to collect and validate user input, document the intentions of
the workflow component, track and report errors, and validate output.

As a simple example let us write a generic energy minimisation script:

.. code-block:: python

   import BioSimSpace as BSS

   # Create the node with a description of what it does.
   node = BSS.Gateway.Node("Perform energy minimisation")

   # Add an author and set the license.
   node.addAuthor(name="Lester Hedges", email="lester.hedges@bristol.ac.uk", affiliation="University of Bristol")
   node.setLicence("GPLv3")

   # Specify the input requirements for the node.
   node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
   node.addInput("steps", BSS.Gateway.Integer(help="The number of minimisation steps.", minimum=0, maximum=100000, default=10000))

   # Specify the output of the node.
   node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system."))

   # Show the control panel so that the user to set the inputs. This will bring
   # up a GUI if running from within Jupyter.
   node.showControls()

   # Create the molecular system and minimisation protocol from the user input.
   system = BSS.IO.readMolecules(node.getInput("files"))
   protocol = BSS.Protocol.Minimisation(steps=node.getInput("steps"))

   # Find a molecular dynamics package and run the protocol.
   process = BSS.MD.run(system, protocol)

   # Wait for the process to finish and write the final molecular configuration to file,
   # binding the file names to the output requirement of the node.
   node.setOutput("minimised", BSS.IO.saveMolecules("minimised", process.getSystem(block=True), system.fileFormat()))

   # Finally, validate that the node completed successfully.
   node.validate()

"""
@package biosimspace
@author  Lester Hedges
@brief   A class for interfacing with BioSimSpace nodes.
"""

from BioSimSpace import _is_notebook

from .requirements import *

from collections import OrderedDict

import argparse
import sys

class Node():
    """A class for interfacing with BioSimSpace nodes."""

    _is_notebook = _is_notebook()

    def __init__(self, description):
        """Constructor.

           Positional arguments:

           description -- A description of the node.
        """

        if type(description) is not str:
            raise TypeError("The 'description' keyword must be of type 'str'.")

        # Set the node description string.
        self._description = description

        # Initialise dictionaries for the inputs/outputs.
        self._inputs = OrderedDict()
        self._outputs = OrderedDict()

        # The input has not yet been validated.
        self._is_validated = False

        # Create the parser.
        self._parser = argparse.ArgumentParser(description=self._description)

        # Add an option to allow the user to load a configuration from file.
        config = File(help="path to a configuration file (optional)", optional=True)
        self.addInput("config", config)

    def addInput(self, name, input):
        """Add an input requirement.

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
        """

        # Can't add requirements if the input has already been validated.
        if self._is_validated:
            return

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'.")

        if not isinstance(input, Requirement):
            raise TypeError("'input' must be of type 'Requirement'.")

        # Add the input to the dictionary.
        self._inputs[name] = input

	# Append long-form argument name if not present.
        if (len(name) > 2):
            if name[0:2] != '--':
                name = '--' + name
        else:
            name = '--' + name

        if input.isOptional():
            if input.getDefault() is not None:
                if input.isMulti() is not False:
                    self._parser.add_argument(name, type=input.getArgType(), nargs='+',
                        help=input.getHelp(), default=input.getDefault())
                else:
                    self._parser.add_argument(name, type=input.getArgType(),
                        help=input.getHelp(), default=input.getDefault())
            else:
                if input.isMulti() is not False:
                    self._parser.add_argument(name, type=input.getArgType(), nargs='+',
                        help=input.getHelp())
                else:
                    self._parser.add_argument(name, type=input.getArgType(),
                        help=input.getHelp())
        else:
            if input.isMulti() is not False:
                self._parser.add_argument(name, type=input.getArgType(), nargs='+',
                    help=input.getHelp(), required=True)
            else:
                self._parser.add_argument(name, type=input.getArgType(),
                    help=input.getHelp(), required=True)

    def addOutput(self, name, output):
        """Add an output requirement.

           Positional arguments:

           name   -- The name of the output.
           output -- The output requirement object.
        """

        # Can't add requirements if the input has already been validated.
        if self._is_validated:
            return

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'.")

        if not isinstance(output, Requirement):
            raise TypeError("'output' must be of type 'Requirement'.")

        # Add the output to the dictionary.
        self._outputs[name] = output

    def setOutput(self, name, value):
        """Set the value of an output.

           Positional arguments:

           name  -- The name of the output.
           value -- The value of the output.
        """

        try:
            self._outputs[name].setValue(value)
        except KeyError:
            raise

    def getInput(self, name):
        """Get the value of the named input.

           Positional arguments:

           name -- The name of the input requirement.
        """

        if not self._is_validated:
            self._validateInputs()
            self._is_validated = True

        if type(name) is not str:
            raise TypeError("The name must be of type 'str'")

        try:
            return self._inputs[name].getValue()
        except KeyError:
            raise

    def getInputs(self):
        """Get all of the input requirements."""

        if not self._is_validated:
            self._validateInputs()
            self._is_validated = True

        return self._inputs

    def validate(self):
        """Whether the output requirements are satisfied."""

        # Check no outputs are None.
        for name, output in self._outputs.items():
            if output.getValue() is None:
                raise SystemExit("Missing output for requirement '%s'" % name)

        # All ouputs are found.
        return True

    def _validateInputs(self):
        """Validate the parsed inputs."""

	# Parse the arguments into a dictionary.
        args = vars(self._parser.parse_args())

        # Now loop over the arguments and set the input values.
        for key, value in args.items():
            self._inputs[key].setValue(value)

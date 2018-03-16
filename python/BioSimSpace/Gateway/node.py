"""
@package biosimspace
@author  Lester Hedges
@brief   A class for interfacing with BioSimSpace nodes.
"""

from BioSimSpace import _is_notebook

from .requirements import *

import argparse

class Node():
    """A class for interfacing with BioSimSpace nodes."""

    _is_notebook = _is_notebook()

    def __init__(self, description):
        """Constructor.

           Positional arguments:

           description -- A description of the node.
        """

        if type(description) is not str:
            raise ValueError("The 'description' keyword must be of type 'str'.")

        # Set the node description string.
        self._description = description

        # Initialise a dictionary of requirements.
        self._requirements = {}

        # Create the parser.
        self._parser = argparse.ArgumentParser(description=self._description)

        # Add an option to allow the user to load a configuration from file.
        config = File(name="config", help="path to a configuration file (optional)", optional=True)
        self.addRequirement(config)

    def setRequirements(self, *args):
        """Set the nodes requirements.

           Positional arguments:

           args -- The requirements for the node.
        """

        # Loop over all of the inputs.
        for arg in args:
            self.addRequirement(arg)

	# Validate the requirements.
        self.validateRequirements()

    def addRequirement(self, requirement):
        """Add a requirement.

           Positional arguments:

           requirement -- A requirement object.
        """

        if not isinstance(requirement, Requirement):
            raise ValueError("The 'requirement' must be of type 'Requirement'.")

	# Get the name of the requirement.
        name = requirement.name()

        # Add the requirement to the dictionary.
        self._requirements[name] = requirement

	# Append long-form argument name if not present.
        if (len(name) > 2):
            if name[0:2] != '--':
                name = '--' + name
        else:
            name = '--' + name

        if requirement.isOptional():
            if requirement.default() is not None:
                if requirement.isMulti() is not False:
                    self._parser.add_argument(name, type=requirement.argType(), nargs='+',
                        help=requirement.helpText(), default=requirement.default())
                else:
                    self._parser.add_argument(name, type=requirement.argType(),
                        help=requirement.helpText(), default=requirement.default())
            else:
                if requirement.isMulti() is not False:
                    self._parser.add_argument(name, type=requirement.argType(), nargs='+',
                        help=requirement.helpText())
                else:
                    self._parser.add_argument(name, type=requirement.argType(),
                        help=requirement.helpText())
        else:
            if requirement.isMulti() is not False:
                self._parser.add_argument(name, type=requirement.argType(), nargs='+',
                    help=requirement.helpText(), required=True)
            else:
                self._parser.add_argument(name, type=requirement.argType(),
                    help=requirement.helpText(), required=True)

    def validateRequirements(self):
        """Validate the parsed requirements."""

	# Parse the arguments into a dictionary.
        args = vars(self._parser.parse_args())

        # Now loop over the arguments and set the requirement values.
        for key, value in args.items():
            self._requirements[key].setValue(value)

    def getRequirement(self, name):
        """Get the value of the named requirement.

           Positional arguments:

           name -- The name of the requirement.
        """

        if type(name) is not str:
            raise ValueError("The name must be of type 'str'")

        try:
            return self._requirements[name].value()
        except KeyError:
            raise

    def getRequirements(self):
        """Get all of the requirements."""

        x = []
        for key, value in self._requirements.items():
            x.append(value)

        return x

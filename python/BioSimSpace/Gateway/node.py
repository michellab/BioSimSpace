"""
@package biosimspace
@author  Lester Hedges
@brief   A class for interfacing with BioSimSpace nodes.
"""

from BioSimSpace import _is_notebook

import argparse

class Node():
    """A class for interfacing with BioSimSpace nodes."""

    _is_notebook = _is_notebook()

    def __init__(self, description=None, parser=None):
        """Constructor.

           Keyword arguments:

           description -- A description of the node.
	   parser      -- An argparse.ArgumentParser object.
        """

        if type(description) is not str:
            raise ValueError("The 'description' keyword must be of type 'str'.")

        # Create the parser object.
        if parser is None:
            self._parser = argparse.ArgumentParser(description=description)

        elif type(parser) is argparse.ArgumentParser:
            self._parser = parser

    def setRequirements(self, inputs=None, outputs=None):
        """Set the nodes requirements.

           Keyword arguments:

           inputs  -- A dictionary of the required inputs.
           outputs -- A dictionary of the required outputs.
        """

        # Add the input requirements.
        if inputs is not None:

            # A single argument dictionary.
            if type(inputs) is dict:
                inputs = [inputs]

            # Make sure all inputs are dicts.
            if all(isinstance(x, dict) for x in inputs):

                # Loop over all of the inputs.
                for input in inputs:

                    # Extract the keyword values.

                    # These are required...

                    try:
                        name = input['name']

                        # Append long-form argument name if not present.
                        if (len(name) > 2):
                            if name[0:2] != '--':
                                name = '--' + name
                        else:
                            name = '--' + name

                    except KeyError:
                        raise("Input requirements must have a 'name' keyword!")

                    try:
                        arg_type = input['type']
                    except KeyError:
                        raise("Input requirements must have a 'type' keyword!")

                    try:
                        doc = input['doc']
                    except KeyError:
                        raise("Input requirements must have a 'doc' keyword!")

                    # These are optional...

                    try:
                        default = input['default']
                    except:
                        default = None

                    try:
                        multi = input['multi']
                    except:
                        multi = False

                    try:
                        required = input['required']
                    except:
                        required = True

                    # Argument is never required if a default is set.
                    if required and default is not None:
                        required = False

                    # Add the argument to the parser.

                    if required is not False:
                        if default is not None:
                            if multi is not False:
                                self._parser.add_argument(name, type=arg_type, nargs='+',
                                    help=doc, default=default, required=True)
                            else:
                                self._parser.add_argument(name, type=arg_type, help=doc,
                                    default=default, required=True)
                        else:
                            if multi is not False:
                                self._parser.add_argument(name, type=arg_type, nargs='+',
                                    help=doc, required=True)
                            else:
                                self._parser.add_argument(name, type=arg_type, help=doc,
                                    required=True)
                    else:
                        if default is not None:
                            if multi is not False:
                                self._parser.add_argument(name, type=arg_type, nargs='+',
                                    help=doc, default=default)
                            else:
                                self._parser.add_argument(name, type=arg_type, help=doc,
                                    default=default)
                        else:
                            if multi is not False:
                                self._parser.add_argument(name, type=arg_type, nargs='+',
                                    help=doc)
                            else:
                                self._parser.add_argument(name, type=arg_type, help=doc)

            # Parse the arguments.
            self._args = self._parser.parse_args()

    def getInput(self, arg):
        """Get the value of a command-line argument."""

        if type(arg) is not str:
            raise ValueError("The arg must be of type 'str'")

        # Convert the arguments to a dictionary.
        arg_dict = vars(self._args)

        try:
            value = arg_dict[arg]
        except KeyError:
            print("Input argument '%s' doesn't exist!" % arg)
            value = None

        return value

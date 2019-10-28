######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for creating BioSimSpace workflow components (nodes).
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Node"]

import configargparse as _argparse
import collections as _collections
import __main__
import os as _os
import sys as _sys
import textwrap as _textwrap
import warnings as _warnings
import yaml as _yaml

from BioSimSpace import _is_notebook
from BioSimSpace import setVerbose

# Enable Jupyter widgets.
if _is_notebook:
    from IPython.display import FileLink as _FileLink

    import fileupload as _fileupload
    import ipywidgets as _widgets
    import zipfile as _zipfile

from BioSimSpace.Types._type import Type as _Type

from ._requirements import Area as _Area
from ._requirements import Boolean as _Boolean
from ._requirements import File as _File
from ._requirements import FileSet as _FileSet
from ._requirements import Float as _Float
from ._requirements import Charge as _Charge
from ._requirements import Energy as _Energy
from ._requirements import Integer as _Integer
from ._requirements import Length as _Length
from ._requirements import Pressure as _Pressure
from ._requirements import Requirement as _Requirement
from ._requirements import String as _String
from ._requirements import Temperature as _Temperature
from ._requirements import Time as _Time
from ._requirements import Volume as _Volume

# Float types (including those with units).
_float_types = [_Float, _Charge, _Energy, _Pressure, _Length, _Area, _Volume,
    _Temperature, _Time]

class Node():
    """A class for interfacing with BioSimSpace nodes.

       Nodes are used to collect and validate user input, document the
       intentions of the workflow component, track and report errors,
       and validate output. Once written, a node can be run from within
       Jupyter, from the command-line, or plugged into a workflow engine,
       such as Knime.

       Example
       -------

       A generic energy minimisation node:

       >>> import BioSimSpace as BSS
       >>> node = BSS.Gateway.Node("Perform energy minimisation")
       >>> node.addAuthor(name="Lester Hedges", email="lester.hedges@bristol.ac.uk", affiliation="University of Bristol")
       >>> node.setLicence("GPLv3")
       >>> node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
       >>> node.addInput("steps", BSS.Gateway.Integer(help="The number of minimisation steps.", minimum=0, maximum=100000, default=10000))
       >>> node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system."))
       >>> node.showControls()
       >>> system = BSS.IO.readMolecules(node.getInput("files"))
       >>> protocol = BSS.Protocol.Minimisation(steps=node.getInput("steps"))
       >>> process = BSS.MD.run(system, protocol)
       >>> node.setOutput("minimised", BSS.IO.saveMolecules("minimised", process.getSystem(block=True), system.fileFormat()))
       >>> node.validate()
    """

    # Whether the node is run from Knime.
    _is_knime = False

    # Whether the node is run from a Jupyter notebook.
    _is_notebook = _is_notebook

    def __init__(self, description, name=None):
        """Constructor.

           Parameters
           ----------

           description : str
               A description of the node.

           name : str
               The name of the node.
        """

        if type(description) is not str:
            raise TypeError("The 'description' keyword must be of type 'str'.")

        # Set the node name.
        if name is None:
            try:
                self._name = _os.path.basename(__main__.__file__)
            except:
                self._name = None
        else:
            if type(name) is not str:
                raise TypeError("The 'name' keyword must be of type 'str'.")
            self._name = name

        # Set the node description string.
        self._description = description

        # Initalise the authors.
        self._authors = None

        # Initalise the license.
        self._license = None

        # Initialise dictionaries for the inputs/outputs.
        self._inputs = _collections.OrderedDict()
        self._outputs = _collections.OrderedDict()

        # A dictionary of Jupyter widgets.
        self._widgets = _collections.OrderedDict()

        # Whether the input is valid.
        self._is_valid_input = False

        # Whether the output has been validated.
        self._is_output_validated = False

        # A list of user error messages.
        self._errors = []

        # Initialise the parser.
        self._parser = None
        self._required = None

        # Intialise the Jupyter input panel.
        self._control_panel = None

        # Running from the command-line.
        if not self._is_knime and not self._is_notebook:
            # Generate the node help description.
            description = self._generate_description()

            # Create the parser.
            self._parser = _argparse.ArgumentParser(description=description,
                formatter_class=_argparse.RawTextHelpFormatter, add_help=False,
                config_file_parser_class=_argparse.YAMLConfigFileParser,
                add_config_file_help=False)

            # Add argument groups.
            self._required = self._parser.add_argument_group("Required arguments")
            self._optional = self._parser.add_argument_group("Optional arguments")
            self._optional.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
            self._optional.add_argument("-c", "--config", is_config_file=True, help="Path to configuration file.")
            self._optional.add_argument("-v", "--verbose", type=_str2bool, nargs='?', const=True, default=False,
                                        help="Print verbose error messages.")

            # Overload the "_check_value" method for more flexible string support.
            # (Ignore whitespace and case insensitive.)
            self._parser._check_value = _check_value

    def __del__(self):
        """Destructor."""

        # Validate the node if the user hasn't already done so.
        if self._is_valid_input:
            if not self._is_output_validated:
                self.validate()

    def addInput(self, name, input):
        """Add an input requirement.

           Parameters
           ----------

           name : str
               The name of the input.

           input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
               The input requirement object.
        """

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'.")

        if not isinstance(input, _Requirement):
            raise TypeError("'input' must be of type 'Requirement'.")

        # We already have an input with this name.
        reset = False
        if name in self._inputs:
            if self._is_notebook:
                _warnings.warn("Duplicate input requirement '%s'"  % name)
                reset = True
            else:
                raise ValueError("Duplicate input requirement '%s'"  % name)

        # Add the input to the dictionary.
        self._inputs[name] = input

        # Create a Knime GUI widget.
        if self._is_knime:
            self._addInputKnime(name, input)

        # Create a Jupyter GUI widget.
        elif self._is_notebook:
            return self._addInputJupyter(name, input, reset)

        # Command-line argparse ArgumentParser.
        else:
            self._addInputCommandLine(name, input)

    def _addInputCommandLine(self, name, input):
        """Add an input requirement for the command-line.

           Parameters
           ----------

           name : str
               The name of the input.

           input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
               The input requirement object.
        """

        # Append long-form argument name if not present.
        if (len(name) > 2):
            if name[0:2] != '--':
                name = '--' + name
        else:
            name = '--' + name

        if input.isOptional():
            if input.getDefault() is not None:
                if input.isMulti() is not False:
                    self._optional.add_argument(name, type=input.getArgType(), nargs='+',
                        help=self._create_help_string(input), default=input.getDefault())
                else:
                    if input.getArgType() is bool:
                        self._optional.add_argument(name, type=_str2bool, nargs='?',
                            const=True, default=input.getDefault(), help=self._create_help_string(input))
                    else:
                        if input.getAllowedValues() is not None:
                            self._optional.add_argument(name, type=input.getArgType(),
                                help=self._create_help_string(input), default=input.getDefault(),
                                choices=input.getAllowedValues())
                        else:
                            self._optional.add_argument(name, type=input.getArgType(),
                                help=self._create_help_string(input), default=input.getDefault())
            else:
                if input.isMulti() is not False:
                    self._optional.add_argument(name, type=input.getArgType(), nargs='+',
                        help=self._create_help_string(input))
                else:
                    self._optional.add_argument(name, type=input.getArgType(),
                        help=self._create_help_string(input))
        else:
            if input.isMulti() is not False:
                self._required.add_argument(name, type=input.getArgType(), nargs='+',
                    help=self._create_help_string(input), required=True)
            else:
                if input.getAllowedValues() is not None:
                    self._required.add_argument(name, type=input.getArgType(),
                        help=self._create_help_string(input), required=True,
                        choices=input.getAllowedValues())
                else:
                    if input.getArgType() is bool:
                        self._required.add_argument(name, type=_str2bool, nargs='?',
                            const=True, help=self._create_help_string(input))
                    else:
                        self._required.add_argument(name, type=input.getArgType(),
                            help=self._create_help_string(input), required=True)

    def _addInputKnime(self, name, input):
        """Add an input requirement for Knime.

           Parameters
           ----------

           name : str
               The name of the input.

           input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
               The input requirement object.
        """
        return None

    def _addInputJupyter(self, name, input, reset=False):
        """Add an input requirement for Jupyter.

           Parameters
           ----------

           name : str
               The name of the input.

           input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
               The input requirement object.

           reset : bool
               Whether to reset the widget data.
        """

        # Create a widget button to indicate whether the requirement value
        # has been set.
        button = _widgets.Button(
            tooltip="The input requirement is unset.",
            button_style="warning",
            icon="fa-exclamation-triangle",
            layout=_widgets.Layout(flex="1 1 auto", width="auto"),
            disabled=False,
        )

        # Add a Jupyter widget for each of the supported requirement types.

        # Boolean.
        if type(input) is _Boolean:
            # Create a Jupyter toggle button.
            widget = _widgets.ToggleButton(
                value=False,
                description=name,
                tooltip=input.getHelp(),
                button_style="",
                icon="check",
                disabled=False
            )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Get the default value.
            default = input.getDefault()

            if default is not None:
                widget.value = default
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Integer.
        elif type(input) is _Integer:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = _widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False
                )

            else:
                # Get the range of the input.
                _min = input.getMin()
                _max = input.getMax()

                # Whether the integer is unbounded.
                is_unbounded = True

                if _min is not None:
                    # Set the default.
                    if default is None:
                        default = _min

                    # Bounded integer.
                    if _max is not None:
                        step=int((_max - _min) / 100)
                        # Create an int slider widget.
                        widget = _widgets.IntSlider(
                            value=default,
                            min=_min,
                            max=_max,
                            step=step,
                            description=name,
                            tooltip=input.getHelp(),
                            continuous_update=False,
                            orientation="horizontal",
                            readout=True,
                            readout_format="d",
                            disabled=False
                        )

                        # Flag that the integer is bounded.
                        is_unbounded = False

                # Unbounded integer.
                if is_unbounded:
                    # Create an integer widget.
                    widget = _widgets.IntText(
                        value=default,
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is not None:
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Float types (including those with units).
        elif type(input) in _float_types:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            # Get the magnitude of types with units.
            if isinstance(default, _Type):
                default = default.magnitude()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                    # Get the magnitude of types with units.
                    if isinstance(default, _Type):
                        default = default.magnitude()

                # Create a dropdown for the list of allowed values.
                widget = _widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False
                )

            else:
                # Get the range of the input.
                _min = input.getMin()
                _max = input.getMax()

                # Get the magnitude of types with units.
                if isinstance(_min, _Type):
                    _min = _min.magnitude()
                if isinstance(_max, _Type):
                    _max = _max.magnitude()

                # Whether the float is unbounded.
                is_unbounded = True

                if _min is not None:
                    # Set the default.
                    if default is None:
                        default = _min

                    # Bounded float.
                    if _max is not None:
                        step=(_max - _min) / 100
                        # Create a float slider widget.
                        widget = _widgets.FloatSlider(
                            value=default,
                            min=_min,
                            max=_max,
                            step=step,
                            description=name,
                            tooltip=input.getHelp(),
                            continuous_update=False,
                            orientation="horizontal",
                            readout=True,
                            readout_format=".1f",
                            disabled=False
                        )

                        # Flag that the float is bounded.
                        is_unbounded = False

                # Unbounded float.
                if is_unbounded:
                    # Create a float widget.
                    widget = _widgets.FloatText(
                        value=default,
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is not None:
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # String.
        elif type(input) is _String:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = _widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False
                )

            else:
                if default is None:
                    # Create a text widget without a default.
                    widget = _widgets.Text(
                        placeholder="Type something",
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )
                else:
                    # Create a text widget.
                    widget = _widgets.Text(
                        value=default,
                        placeholder="Type something",
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is not None:
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # File.
        elif type(input) is _File:
            # Create a fileupload widget.
            widget = _fileupload.FileUploadWidget()

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Flag that this is just a single file upload.
            widget._is_multi = False

            # Set the value to None.
            widget.value = None

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_file_upload, names="data")

            # Store the widget.
            self._widgets[name] = widget

        # File set.
        elif type(input) is _FileSet:
            # Create a fileupload widget.
            widget = _fileupload.FileUploadWidget()

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Flag that this is is a set of files.
            widget._is_multi = True

            # Set the value to None.
            widget.value = None

            # Store a reference to the node.
            widget._node = self

            # Store the requirement name.
            widget._name = name

            # Store the requirement.
            widget._input = input

            # Bind the callback function.
            widget.observe(_on_file_upload, names="data")

            # This is a new widget.
            if not name in self._widgets or reset:
                self._widgets[name] = [widget]
            else:
                self._widgets[name].append(widget)

        # Unsupported input.
        else:
            raise ValueError("Unsupported requirement type '%s'" % type(input))

    def addOutput(self, name, output):
        """Add an output requirement.

           Parameters
           ----------

           name : str
               The name of the output.

           output : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
               The output requirement object.
        """

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'.")

        if not isinstance(output, _Requirement):
            raise TypeError("'output' must be of type 'Requirement'.")

        # We already have an ouput requirement with this name.
        if name in self._outputs:
            _warnings.warn("Duplicate input requirement. Overwriting existing value!")

        # Add the output to the dictionary.
        self._outputs[name] = output

        # Update the parser description.
        if not self._is_notebook:
            self._parser.description = self._generate_description()

    def setOutput(self, name, value):
        """Set the value of an output.

           Parameters
           ----------

           name : str
               The name of the output.

           value :
               The value of the output.
        """
        try:
            self._outputs[name].setValue(value, name=name)
        except KeyError:
            raise

    def getInput(self, name):
        """Get the value of the named input.

           Parameters
           ----------

           name : str
               The name of the input requirement.

           Returns
           -------

           input :
               The value of the named input requirement.
        """

        if type(name) is not str:
            raise TypeError("The name must be of type 'str'")

        # Validate the inputs.
        self._is_valid_input = self._validateInput()

        try:
            value = self._inputs[name].getValue()
            if type(value) is list:
                return value.copy()
            else:
                return value
        except KeyError:
            raise

    def getInputs(self):
        """Get all of the input requirements.

           Returns
           -------

           inputs : { str : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>` }
               The dictionary of input requirements.
        """

        # Validate the inputs.
        self._is_valid_input = self._validateInput()

        return self._inputs.copy()

    def addError(self, error):
        """Add an error message.

           Parameters
           ----------

           error : str
               The error message.
        """

        if type(error) is not str:
            raise TypeError("The error message must be of type 'str'")
        else:
            self._errors.append(error)

    def addAuthor(self, name=None, email=None, affiliation=None):
        """Add an author for the node.

           Parameters
           ----------

           name : str
               The author's name.

           email : str
               The author's email address.

           affiliation : str
               The author's affiliation.
        """

        if name is None:
            raise ValueError("Missing required 'name' keyword argument.")

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'")

        if email is not None and type(email) is not str:
            raise TypeError("'email' must be of type 'str'")

        if affiliation is not None and type(affiliation) is not str:
            raise TypeError("'affiliation' must be of type 'str'")

        if self._authors is None:
            self._authors = [{"name" : name, "email" : email, "affiliation" : affiliation}]
        else:
            author = {"name" : name, "email" : email, "affiliation" : affiliation}
            if not author in self._authors:
                self._authors.append(author)

    def getAuthors(self):
        """Return the list of authors.

           Returns
           -------

           authors : [dict]
              A list of author dictionaries.
        """
        return self._authors.copy()

    def setLicense(self, license):
        """Set the license for the node.

           Parameters
           ----------

           license : str
               The license type.
        """
        if type(license) is not str:
            raise TypeError("The license must be of type 'str'")
        else:
            self._license = license

    def getLicense(self):
        """Return the license.

           Returns
           -------

           license : str
               The license of the node.
        """
        return self._license

    def showControls(self):
        """Show the Jupyter widget GUI to allow the user to enter input.

           Returns
           -------

           controls : ipywidgets.form
              A gui control panel for setting input requirements.
        """

        if not self._is_notebook:
            return

        # Create the layout object.
        layout = _widgets.Layout(
            display="flex",
            flex_flow="row",
            justify_content="space-between"
        )

        # Initialise the list of form items.
        requirements = []

        # Initialise the list of input requirement 'set' indicators.
        indicators = []

        # Loop over all of the widgets.
        for name, widget in self._widgets.items():

            # Credate the label string.
            string = "%s: %s" % (name, self._inputs[name].getHelp())

            # Add the unit information.
            unit = self._inputs[name].getUnit()
            if unit is not None:
                string += " (%s)" % self._inputs[name]._print_unit

            # Create the widget label.
            label = _widgets.Label(value=string)

            # This is a FileSet requirement with multiple widgets.
            if type(widget) is list:
                items = [label] + widget
                indicator = widget[0]._button
            else:
                items = [label, widget]
                indicator = widget._button

            # Create a box for the widget.
            box = _widgets.Box(items, layout=layout)

            # Add the widget box to the list.
            requirements.append(box)

            # Add the indicator to the list.
            indicators.append(indicator)

        # Create the widget form.
        form1 = _widgets.Box(requirements, layout=_widgets.Layout(
            display="flex",
            flex_flow="column",
            border="solid 2px",
            align_items="stretch",
            width="100%"
        ))

        # Create the indicator form.
        form2 = _widgets.VBox(indicators, layout=_widgets.Layout(
            display="flex",
            flex_flow="column",
            align_items="stretch",
            width="5%"
        ))

        # Combine the two forms.
        form = _widgets.Box([form1, form2], layout=layout)

        # Store the form.
        self._control_panel = form

        return form

    def _validateInput(self):
        """Validate the parsed inputs."""

        # Knime.
        if self._is_knime:
            pass

        # Jupyter.
        elif self._is_notebook:
            # Loop over the widgets and set the input values.
            for key, widget in self._widgets.items():

                # Use the widget value if it has been set, otherwise, set the value to None.
                # This ensures that the user actually sets a value.

                # This is a FileSet requirement with multiple widgets.
                if type(widget) is list:
                    value = []
                    # Loop over all of the widgets.
                    for w in widget:
                        if w._is_set:
                            value.append(w.value)
                    # If there are no values, set to None.
                    if len(value) == 0:
                        value = None
                # Single widget.
                else:
                    if widget._is_set:
                        value = widget.value
                    else:
                        value = None

                self._inputs[key].setValue(value, name=key)

        # Command-line.
        else:
            # Parse the arguments into a dictionary.
            args = vars(self._parser.parse_known_args()[0])

            # Now loop over the arguments and set the input values.
            for key, value in args.items():
                if key is "verbose":
                    setVerbose(value)
                else:
                    if key is not "config":
                        self._inputs[key].setValue(value, name=key)

    def validate(self, file_prefix="output"):
        """Whether the output requirements are satisfied.

           Parameters
           ----------

           file_prefix : str
               The prefix of the output file name.

           Returns
           -------

           output : IPython.lib.display.FileLink, str
               If running interatvely: A link to a zipfile containing the
               validated output, else the name of a YAML file containing
               the node output.
        """

        if type(file_prefix) is not str:
            raise TypeError("The 'file_prefix' keyword must be of type 'str'.")

        # Flag that we have validated output.
        self._is_output_validated = True

        # A list of File and FileSet outputs.
        file_outputs = []

        # Check no outputs are None.
        for name, output in self._outputs.items():
            if output.getValue() is None:
                self._errors.append("Missing output for requirement '%s'" % name)
            else:
                if type(output) is _File or type(output) is _FileSet:
                    file_outputs.append(output)

        # Node failed.
        if len(self._errors) > 0:
            for error in self._errors:
                print("%s" % error, file=_sys.stderr)

            if self._name is not None:
                raise SystemExit("Node '%s' failed!" % self._name)
            else:
                raise SystemExit("Node failed!")

        # Create a compressed archive containing all file output for the node.
        if self._is_notebook:
            # There are files.
            if len(file_outputs) > 0:
                # Create the archive name.
                zipname = "%s.zip" % file_prefix

                # Append the files to the archive.
                with _zipfile.ZipFile(zipname, "w") as zip:
                    # Loop over all of the file outputs.
                    for output in file_outputs:
                        if type(output) is _File:
                            file = output.getValue()
                            zip.write(file, arcname=_os.path.basename(file))
                        else:
                            for file in output.getValue():
                                zip.write(file, arcname=_os.path.basename(file))

                # Return a link to the archive.
                return _FileLink(zipname)
        else:
            # Initialise an empty dictionary to store the output data.
            data = {}

            # Create the YAML file name.
            yamlname = "%s.yaml" % file_prefix

            # Populate the dictionary.
            for name, output in self._outputs.items():
                data[name] = output.getValue()

            # Write the outputs to a YAML file.
            with open(yamlname, "w") as file:
                _yaml.dump(data, file, default_flow_style=False)

            return yamlname

    def _create_help_string(self, input):
        """Create a nicely formatted argparse help string.

           Parameters
           ----------

           input : :class:`Requirement <BioSimSpace.Gateway.Requirement>`
               The input requirement.

           Returns
           -------

           help : str
               The formatted help string.
        """

        # Initialise the help string.
        help = "\n".join(_textwrap.wrap(input.getHelp(), 56))

        # Get the unit.
        units = input.getUnit()

        # Add the units to the help string.
        if units is not None:
            help += "\n  units=%s" % (units[0] + units[1:].lower())

        # Get the default value.
        default = input.getDefault()

        # Add default value to help string.
        if default is not None:
            help += "\n  default=%s" % default

        # Get the min/max allowed values.
        minimum = input.getMin()
        maximum = input.getMax()

        # Add min/max values to help string.
        if minimum is not None:
            if maximum is not None:
                help += "\n  min=%s, max=%s" % (minimum, maximum)
            else:
                help += "\n  min=%s" % minimum
        else:
            if maximum is not None:
                help += "\n  max=%s" % maxmimum

        return help

    def _generate_description(self):
        """Generate a formatted output section for the argparse help description.

           Returns
           -------

           output : str
               A string listing the output requirements.
        """

        string = "\n".join(_textwrap.wrap(self._description, 80))

        yaml_help = ("Args that start with '--' (e.g. --arg) can also be set "
                     "in a config file (specified via -c). The config file uses "
                     "YAML syntax and must represent a YAML 'mapping' "
                     "(for details, see http://learn.getgrav.org/advanced/yaml). "
                     "If an arg is specified in more than one place, then "
                     "commandline values override config file values which "
                     "override defaults."
                    )

        string += "\n\n" + "\n".join(_textwrap.wrap(yaml_help, 80))

        # Initialise the output string.
        string += "\n\nOutput:\n"

        # Add documentation for each output.
        for name, output in self._outputs.items():
            num_whitespace = 19 - len(output.__class__.__name__) - len(name)
            s = "  %s: %s " % (name, output.__class__.__name__)
            for i, line in enumerate(_textwrap.wrap(output.getHelp(), 56)):
                if i == 0:
                    for _ in range(0, num_whitespace):
                        s += " "
                else:
                    for _ in range(0, 24):
                        s += " "
                s +=  line + "\n"
            string += s

        return string

def _on_value_change(change):
    """Helper function to flag that a widget value has been set."""
    change["owner"]._is_set = True
    change["owner"]._button.tooltip = "The input requirement is set."
    change["owner"]._button.button_style = "success"
    change["owner"]._button.icon = "fa-check"

def _on_file_upload(change):
    """Helper function to handle file uploads."""

    # Store the number of bytes.
    num_bytes = len(change["owner"].data)

    # Return if there is no data.
    if num_bytes == 0:
        return

    # Get the file name.
    filename = change["owner"].filename

    # Create the uploads directory if it doesn't already exist.
    if not _os.path.isdir("uploads"):
        _os.makedirs("uploads")

    # Append the upload directory to the file name.
    new_filename = "uploads/%s" % filename

    # Has this file already been uploaded?
    if _os.path.isfile(new_filename):

        # We'll append a number to the file name.
        index = 1
        new_filename_append = new_filename + ".%d" % index

        # Keep trying until a unique name is found.
        while _os.path.isfile(new_filename_append):
            index += 1
            new_filename_append = new_filename + ".%d" % index

        # Copy back into the new_filename variable.
        new_filename = new_filename_append

    # Write the file to disk.
    with open(new_filename, "wb") as file:
        file.write(change["owner"].data)

    # Report that the file was uploaded.
    print("Uploaded '{}' ({:.2f} kB)".format(
        filename, num_bytes / 2 **10))

    # Clear the redundant data from the widget.
    change["owner"].data = b""

    # Flag that the widget value has been set.
    change["owner"]._is_set = True
    change["owner"]._button.tooltip = "The input requirement is set."
    change["owner"]._button.button_style = "success"
    change["owner"]._button.icon = "fa-check"

    # Now update the widget value.

    # Truncate the filename string if it is more than 15 characters.
    label = (filename[:15] + "...") if len(filename) > 15 else filename

    # This is a file set widget.
    if change["owner"]._is_multi:
        # This is the first time the value has been set.
        if change["owner"].value is None:
            # Whether a match has been found.
            is_match = False

            # Store the name of the input requirement.
            name = change["owner"]._name

            # Loop over the widgets in the control panel and find the one
            # that with the matching name.
            for index, child in enumerate(change["owner"]._node._control_panel.children[0].children):
                # The widget name matches.
                if child.children[1]._name == name:
                    is_match = True
                    break

                # Increment the index.
                index += 1

            # No match!
            if not is_match:
                raise RunTimeError("Missing widget for requirement name: '%s'" % name)

            # Create a new widget.
            change["owner"]._node._addInputJupyter(name, change["owner"]._input)

            # Convert the children of the control panel to a list.
            boxes = list(change["owner"]._node._control_panel.children[0].children[index].children)

            # Append the new widget to the list.
            boxes.append(change["owner"]._node._widgets[name][-1])

            # Add the updated box back into the list of boxes.
            change["owner"]._node._control_panel.children[0].children[index].children = tuple(boxes)

    # Update the widget value and label.
    change["owner"].value = new_filename
    change["owner"].label = label

def _check_value(action, value):
    """Helper function to overload argparse's choice checker."""
    if action.choices is not None and value not in action.choices:
        args = {"value"   : value,
                "choices" : ", ".join(map(repr, action.choices))
               }
        msg = _argparse.argparse._("invalid choice: %(value)r (choose from %(choices)s)")

        # If the value is a string, then strip whitespace and try a case insensitive search.
        if type(value) is str:
            new_value = value.replace(" ", "").upper()
            choices = [x.replace(" ", "").upper() for x in action.choices]

            # Check whether we now have a match.
            if new_value not in choices:
                raise _argparse.ArgumentError(action, msg % args)
        else:
            raise _argparse.ArgumentError(action, msg % args)

def _str2bool(v):
    """Convert an argument string to a boolean value."""
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise _argparse.ArgumentTypeError("Boolean value expected.")

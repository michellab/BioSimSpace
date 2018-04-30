######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from BioSimSpace import _is_notebook

from ._requirements import Boolean as _Boolean
from ._requirements import File as _File
from ._requirements import FileSet as _FileSet
from ._requirements import Float as _Float
from ._requirements import Integer as _Integer
from ._requirements import Requirement as _Requirement
from ._requirements import String as _String

import argparse as _argparse
import collections as _collections
import io as _io
import __main__
import os as _os
import sys as _sys
import warnings as _warnings

# Enable Jupyter widgets.
if _is_notebook():
    from IPython.display import FileLink as _FileLink

    import fileupload as _fileupload
    import ipywidgets as _widgets
    import zipfile as _zipfile

__all__ = ["Node"]

class Node():
    """A class for interfacing with BioSimSpace nodes."""

    # Whether the node is run from Knime.
    _is_knime = False

    # Whether the node is run from a Jupyter notebook.
    _is_notebook = _is_notebook()

    def __init__(self, description, name=None):
        """Constructor.

           Positional arguments:

           description -- A description of the node.

           Keyword arguments:

           name        -- The name of the node.
        """

        if type(description) is not str:
            raise TypeError("The 'description' keyword must be of type 'str'.")

        # Set the node name.
        if name is None:
            try:
                self._name = _os.path.basename(__main__.__file__)
            except:
                self._name = None

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
        self._required_args = None

        # Intialise the Jupyter input panel.
        self._control_panel = None

        # Running from the command-line.
        if not self._is_knime and not self._is_notebook:
            # Create the parser.
            self._parser = _argparse.ArgumentParser(description=self._description,
                formatter_class=_argparse.ArgumentDefaultsHelpFormatter)

            # TODO: Add an option to allow the user to load a configuration from file.
            # config = File(help="path to a configuration file", optional=True)
            # self.addInput("config", config)

    def __del__(self):
        """Destructor."""

        # Validate the node if the user hasn't already done so.
        if self._is_valid_input:
            if not self._is_output_validated:
                self.validate()

    def addInput(self, name, input):
        """Add an input requirement.

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
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

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
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
                    self._parser.add_argument(name, type=input.getArgType(), nargs='+',
                        help=input.getHelp(), default=input.getDefault())
                else:
                    if input.getArgType() is bool:
                        self._parser.add_argument(name, type=_str2bool, nargs='?',
                            const=True, default=input.getDefault(), help=input.getHelp())
                    else:
                        if input.getAllowedValues() is not None:
                            self._parser.add_argument(name, type=input.getArgType(), help=input.getHelp(),
                                default=input.getDefault(), choices=input.getAllowedValues())
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
            # Create the required arguments group.
            if self._required_args is None:
                self._required_args = self._parser.add_argument_group("required arguments")

            if input.isMulti() is not False:
                self._required_args.add_argument(name, type=input.getArgType(), nargs='+',
                    help=input.getHelp(), required=True)
            else:
                if input.getAllowedValues() is not None:
                    self._required_args.add_argument(name, type=input.getArgType(),
                        help=input.getHelp(), required=True, choices=input.getAllowedValues())
                else:
                    if input.getArgType() is bool:
                        self._required_args.add_argument(name, type=_str2bool, nargs='?',
                            const=True, help=input.getHelp())
                    else:
                        self._required_args.add_argument(name, type=input.getArgType(),
                            help=input.getHelp(), required=True)

    def _addInputKnime(self, name, input):
        """Add an input requirement for Knime.

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
        """
        return None

    def _addInputJupyter(self, name, input, reset=False):
        """Add an input requirement for Jupyter.

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
           reset -- Whether to reset the widget data.
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
                min_ = input.getMin()
                max_ = input.getMax()

                # Whether the integer is unbounded.
                is_unbounded = True

                if min_ is not None:
                    # Set the default.
                    if default is None:
                        default = min_

                    # Bounded integer.
                    if max_ is not None:
                        # Create an int slider widget.
                        widget = _widgets.IntSlider(
                            value=default,
                            min=min_,
                            max=max_,
                            step=1,
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

        # Float.
        elif type(input) is _Float:
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
                min_ = input.getMin()
                max_ = input.getMax()

                # Whether the float is unbounded.
                is_unbounded = True

                if min_ is not None:
                    # Set the default.
                    if default is None:
                        default = min_

                    # Bounded float.
                    if max_ is not None:
                        # Create a float slider widget.
                        widget = _widgets.FloatSlider(
                            value=default,
                            min=min_,
                            max=max_,
                            step=0.1,
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

           Positional arguments:

           name   -- The name of the output.
           output -- The output requirement object.
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
        """Get all of the input requirements."""

        # Validate the inputs.
        self._is_valid_input = self._validateInput()

        return self._inputs.copy()

    def addError(self, error):
        """Add an error message."""

        if type(error) is not str:
            raise TypeError("The error message must be of type 'str'")
        else:
            self._errors.append(error)

    def addAuthor(self, name=None, email=None, affiliation=None):
        """Add an author for the node.

           Keyword arguments:

           name        -- The author's name.
           email       -- The author's email address.
           affiliation -- The author's affiliation.
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
        """Return the list of authors."""
        return self._authors.copy()

    def setLicense(self, license):
        """Set the license for the node."""

        if type(license) is not str:
            raise TypeError("The license must be of type 'str'")
        else:
            self._license = license

    def getLicense(self):
        """Return the license."""
        return self._license

    def showControls(self):
        """Show the Jupyter widget GUI to allow the user to enter input."""

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
            # Create the widget label.
            label = _widgets.Label(value="%s: %s" % (name, self._inputs[name].getHelp()))

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

                self._inputs[key].setValue(value)

        # Command-line.
        else:
            # Parse the arguments into a dictionary.
            args = vars(self._parser.parse_args())

            # Now loop over the arguments and set the input values.
            for key, value in args.items():
                self._inputs[key].setValue(value)

    def validate(self):
        """Whether the output requirements are satisfied."""

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
                if self._name is None:
                    zipname = "output.zip"
                else:
                    zipname = "%s.zip" % self._name

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
            return True

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
            is_found = False

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

def _str2bool(v):
    """Convert an argument string to a boolean value."""
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise _argparse.ArgumentTypeError("Boolean value expected.")

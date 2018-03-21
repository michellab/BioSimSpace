"""
@package biosimspace
@author  Lester Hedges
@brief   A class for interfacing with BioSimSpace nodes.
"""

from BioSimSpace import _is_notebook

from Sire import try_import

# Enable Jupyter widgets.
if _is_notebook():
    try:
        widgets = try_import("ipywidgets")
    except ImportError:
        raise ImportError("Ipywidgets is not installed. Please install ipywidgets in order to use BioSimSpace.")
    try:
        fileupload = try_import("fileupload")
    except ImportError:
        raise ImportError("Fileupload is not installed. Please install fileupload in order to use BioSimSpace.")

from .requirements import *

from collections import OrderedDict
from os import makedirs
from os.path import basename
from warnings import warn

import argparse
import io
import __main__ as main
import sys

class Node():
    """A class for interfacing with BioSimSpace nodes."""

    # Whether the node is run from Knime.
    _is_knime = False

    # Whether the node is run from a Jupyter notebook.
    _is_notebook = _is_notebook()

    def __init__(self, description):
        """Constructor.

           Positional arguments:

           description -- A description of the node.
        """

        if type(description) is not str:
            raise TypeError("The 'description' keyword must be of type 'str'.")

        # Set the node name.
        try:
            self._name = basename(main.__file__)
        except:
            self._name = None

        # Set the node description string.
        self._description = description

        # Initalise the authors.
        self._authors = None

        # Initalise the license.
        self._license = None

        # Initialise dictionaries for the inputs/outputs.
        self._inputs = OrderedDict()
        self._outputs = OrderedDict()

        # A dictionary of Jupyter widgets.
        self._widgets = OrderedDict()

        # Whether the output has been validated.
        self._is_output_validated = False

        # A list of user error messages.
        self._errors = []

        # Initialise the parser.
        self._parser = None

        # Running from the command-line.
        if not self._is_knime and not self._is_notebook:
            # Create the parser.
            self._parser = argparse.ArgumentParser(description=self._description)

            # Add an option to allow the user to load a configuration from file.
            config = File(help="path to a configuration file (optional)", optional=True)
            self.addInput("config", config)

    def __del__(self):
        """Destructor."""

        # Validate the node if the user hasn't already done so.
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

        if not isinstance(input, Requirement):
            raise TypeError("'input' must be of type 'Requirement'.")

        # We already have an input with this name.
        if name in self._inputs:
            warn("Duplicate input requirement '%s'"  % name)

        # Add the input to the dictionary.
        self._inputs[name] = input

        # Create a Knime GUI widget.
        if self._is_knime:
            self._addInputKnime(name, input)

        # Create a Jupyter GUI widget.
        elif self._is_notebook:
            return self._addInputJupyter(name, input)

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

    def _addInputKnime(self, name, input):
        """Add an input requirement for Knime.

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
        """
        return None

    def _addInputJupyter(self, name, input):
        """Add an input requirement for Jupyter.

           Positional arguments:

           name  -- The name of the input.
           input -- The input requirement object.
        """

        # Add a Jupyter widget for each of the supported requirement types.

        # Boolean.
        if type(input) is Boolean:
            # Create a Jupyter toggle button.
            widget = widgets.ToggleButton(
                value=False,
                description=name,
                tooltip=input.getHelp(),
                button_style='',
                icon='check',
                disabled=False
            )

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is None:
                widget._is_set = False

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Integer.
        elif type(input) is Integer:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = widgets.Dropdown(
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
                        widget = widgets.IntSlider(
                            value=default,
                            min=min_,
                            max=max_,
                            step=1,
                            description=name,
                            tooltip=input.getHelp(),
                            continuous_update=False,
                            orientation='horizontal',
                            readout=True,
                            readout_format='d',
                            disabled=False
                        )

                        # Flag that the integer is bounded.
                        is_unbounded = False

                # Unbounded integer.
                if is_unbounded:
                    # Create an integer widget.
                    widget = widgets.IntText(
                        value=default,
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is None:
                widget._is_set = False

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Float.
        elif type(input) is Float:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = widgets.Dropdown(
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
                        widget = widgets.FloatSlider(
                            value=default,
                            min=min_,
                            max=max_,
                            step=0.1,
                            description=name,
                            tooltip=input.getHelp(),
                            continuous_update=False,
                            orientation='horizontal',
                            readout=True,
                            readout_format='.1f',
                            disabled=False
                        )

                        # Flag that the float is bounded.
                        is_unbounded = False

                # Unbounded float.
                if is_unbounded:
                    # Create a float widget.
                    widget = widgets.IntText(
                        value=default,
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is None:
                widget._is_set = False

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # String.
        elif type(input) is String:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False
                )

            else:
                if default is None:
                    # Create a text widget without a default.
                    widget = widgets.Text(
                        placeholder='Type something',
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )
                else:
                    # Create a text widget.
                    widget = widgets.Text(
                        value=default,
                        placeholder='Type something',
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False
                    )

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is None:
                widget._is_set = False

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # File.
        elif type(input) is File:
            # Create a fileupload widget.
            widget = fileupload.FileUploadWidget()

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            widget._is_set = False

            # Flag that this is just a single file upload.
            widget._is_multi = False

            # Set the value to None.
            widget.value = None

            # Bind the callback function.
            widget.observe(_on_file_upload, names="data")

            # Store the widget.
            self._widgets[name] = widget

        # File set.
        elif type(input) is FileSet:
            # Create a fileupload widget.
            widget = fileupload.FileUploadWidget()

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            widget._is_set = False

            # Flag that this is is a set of files.
            widget._is_multi = True

            # Set the value to None.
            widget.value = None

            # Bind the callback function.
            widget.observe(_on_file_upload, names="data")

            # Store the widget.
            self._widgets[name] = widget

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

        if not isinstance(output, Requirement):
            raise TypeError("'output' must be of type 'Requirement'.")

        # We already have an ouput requirement with this name.
        if name in self._outputs:
            warn("Duplicate input requirement. Overwriting existing value!")

        # Add the output to the dictionary.
        self._outputs[name] = output

    def setInput(self, name, value):
        """Directly set an input requirement.

           Positional arguments:

           name  -- The name of the output.
           value -- The value of the output.
        """

        try:
            self._inputs[name].setValue(value)
        except KeyError:
            raise

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
        self._validateInput()

        try:
            return self._inputs[name].getValue()
        except KeyError:
            raise

    def getInputs(self):
        """Get all of the input requirements."""

        # Validate the inputs.
        self._validateInput()

        return self._inputs

    def addError(self, error):
        """Add an error message."""

        if type(error) is not str:
            raise TypeError("The error message must be of type 'str'")
        else:
            self._errors.append(error)

    def addAuthor(self, author):
        """Set the author of the node."""

        if type(author) is not str:
            raise TypeError("The author must be of type 'str'")
        else:
            self._authors.append(error)

    def getAuthors(self):
        """Return the list of authors."""

        if len(self._authors) == 0:
            return None
        else:
            return self._authors

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

        # Update the widgets if values have been set.
        self._update_widgets()

        # Create the layout object.
        layout = widgets.Layout(
            display='flex',
            flex_flow='row',
            justify_content='space-between'
        )

        # Initialise the list of form items.
        form_items = []

        # Loop over all of the widgets.
        for name, widget in self._widgets.items():
            # Create the widget label.
            label = widgets.Label(value="%s: %s" % (name, self._inputs[name].getHelp()))

            # Create a box for the widget.
            box = widgets.Box([label, widget], layout=layout)

            # Add the box to the list of items.
            form_items.append(box)

        # Create the form.
        form = widgets.Box(form_items, layout=widgets.Layout(
            display='flex',
            flex_flow='column',
            border='solid 2px',
            align_items='stretch',
            width='100%'
        ))

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

        # Check no outputs are None.
        for name, output in self._outputs.items():
            if output.getValue() is None:
                self._errors.append("Missing output for requirement '%s'" % name)

        # Node failed.
        if len(self._errors) > 0:
            for error in self._errors:
                print("%s" % error, file=sys.stderr)

            if self._name is not None:
                raise SystemExit("Node '%s' failed!" % self._name)
            else:
                raise SystemExit("Node failed!")

        # Node completed successfully.
        return True

    def _update_widgets(self):
        """Update widget defaults if values have been set."""

        for name, widget in self._widgets.items():
            value = self._inputs[name].getValue()
            if value is not None:
                widget.value = value

def _on_value_change(change):
    """Helper function to flag that a widget value has been set."""
    change['owner']._is_set = True

def _on_file_upload(change):
    """Helper function to handle file uploads."""

    # Decode the data stream.
    decoded = io.StringIO(change['owner'].data.decode('utf-8'))

    # Get the file name.
    filename = change['owner'].filename

    # Print that the file was uploaded.
    print('Uploaded `{}` ({:.2f} kB)'.format(
        filename, len(decoded.read()) / 2 **10))

    # Create the uploads directory if it doesn't already exist.
    if not path.isdir("uploads"):
        makedirs("uploads")

    # Append the upload directory to the file name.
    new_filename = "uploads/%s" % filename

    # Write the file to disk.
    with open(new_filename, 'w') as file:
        file.write(change['owner'].data.decode('utf-8'))

    # Flag that the widget value has been set.
    change['owner']._is_set = True

    # Now update the widget value.

    # Truncate the filename string if it is more than 15 characters.
    label = (filename[:15] + '...') if len(filename) > 15 else filename

    # This is a file set widget.
    if change['owner']._is_multi:
        # No previous value set.
        if change['owner'].value is None:
            change['owner'].value = [new_filename]
            change['owner'].label = "[%s]" % label
        else:
            change['owner'].value.append(new_filename)
            change['owner'].label = \
                change['owner'].label.replace(']', '') + ", %s]" % label
    # This is a file widget.
    else:
        change['owner'].value = new_filename
        change['owner'].label = label

def _str2bool(v):
    """Convert an argument string to a boolean value."""
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

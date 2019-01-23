import BioSimSpace as BSS

from collections import OrderedDict

import os
import pytest

# Make sure AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_minimise():
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_equilibrate():
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_heat():
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(0, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(300, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_cool():
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(300, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(0, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_production():
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_args():
    """Test setting an manipulation of command-line args."""

    # Create a default minimisation protocol. This doesn't matter since
    # we're going to clear the default arguments anyway.
    protocol = BSS.Protocol.Minimisation()

    # Create the process object.
    process = create_process(protocol)

    # Clear the existing arguments.
    process.clearArgs()

    # Now add some arguments. Firstly one-by-one, using a mixture of
    # arguments and flags.
    process.setArg("-a", "A")   # Regular argument.
    process.setArg("-b", "B")   # Regular argument.
    process.setArg("-c", True)  # Boolean flag.
    process.setArg("-d", "D")   # Regular argument.
    process.setArg("-e", True)  # Boolean flag.
    process.setArg("-f", 6)     # Argument value is an integer.

    # Get the arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 6

    # Make sure the string is correct.
    assert len(arg_string_list) == 10
    assert arg_string == "-a A -b B -c -d D -e -f 6"

    # Turn off one of the flags.
    process.setArg("-c", False)

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the same number of arguments.
    assert len(args) == 6

    # Make sure the new string is correct. (The "False" flag should be missing.)
    assert len(arg_string_list) == 9
    assert arg_string == "-a A -b B -d D -e -f 6"

    # Create a new dictionary of extra arguments. This could be a regular
    # dictionary, but we use an OrderedDict for testing purposes.
    extra_args = OrderedDict()

    # Populate the arguments.
    extra_args["-g"] = True
    extra_args["-h"] = "H"
    extra_args["-i"] = False
    extra_args["-k"] = "K"

    # Add the arguments.
    process.addArgs(extra_args)

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 10

    # Make sure the new string is correct.
    assert len(arg_string_list) == 14
    assert arg_string == "-a A -b B -d D -e -f 6 -g -h H -k K"

    # Now we'll delete an argument.
    process.deleteArg("-d")

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 9

    # Make sure the new string is correct.
    assert len(arg_string_list) == 12
    assert arg_string == "-a A -b B -e -f 6 -g -h H -k K"

    # Now test insertion of additional arguments.
    process.insertArg("-x", "X", 0)    # Insert at beginning.
    process.insertArg("-y", True, 4)   # Insert at middle.
    process.insertArg("-z", "Z", 11)   # Insert at end.

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 12

    # Make sure the new string is correct.
    assert len(arg_string_list) == 17
    assert arg_string == "-x X -a A -b B -y -e -f 6 -g -h H -k K -z Z"

def create_process(protocol):
    """Create an Amber process for a given prototol."""

    # Glob the input files.
    files = BSS.IO.glob("test/io/amber/ala/*")

    # Load the molecular system.
    system = BSS.IO.readMolecules(files)

    # Initialise the AMBER process.
    return BSS.Process.Amber(system, protocol, name="test")

def run_process(protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the AMBER process.
    process = create_process(protocol)

    # Start the AMBER simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Return the process exit code.
    return not process.isError()

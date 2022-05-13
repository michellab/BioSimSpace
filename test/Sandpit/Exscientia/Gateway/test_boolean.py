from BioSimSpace.Gateway import Boolean

import pytest

# Tests for Boolean requirements.

def test_no_arguments():
    """Test that calling constructor with no arguments will raise a ValueError."""

    with pytest.raises(ValueError):
        b = Boolean()

@pytest.mark.parametrize("value", [True, False])
def test_value(value):
    """Test whether object is initialised correctly and value is set."""

    # Create boolean with name and help text.
    b = Boolean(help="Help!")

    # Set the value of the boolean.
    b.setValue(value)

    # Assert that all arguments are parsed correctly and the value is set.
    assert b.getValue() == value
    assert b.getHelp() == "Help!"
    assert b.getDefault() is None
    assert b.isOptional() is False
    assert b.getMin() == None
    assert b.getMax() == None
    assert b.getAllowedValues() is None
    assert b.isMulti() == False
    assert b.getArgType() == bool

@pytest.mark.parametrize('default', [True, False])
def test_default(default):
    """Test that default value is set correctly."""

    b = Boolean(help="Help!", default=default)

    assert b.getHelp() == "Help!"
    assert b.getDefault() == default

    # Make sure that the argument is now optional.
    assert b.isOptional() is True

@pytest.mark.parametrize("value", [5, 'True', 5.7, 'string', [2, 3, 4]])
def test_bad_value(value):
    """Test that TypeError is raised when "value" is not a boolean."""

    b = Boolean(help="Help!")
    with pytest.raises(TypeError):
        b.setValue(value)

@pytest.mark.parametrize('default', [5, 'False', 5.7, 'string', [2, 3, 4]])
def test_bad_defaul(default):
    """Test that TypeError is raised when the default is the wrong type."""

    with pytest.raises(TypeError):
        b = Boolean(help="Help!", default=default)

@pytest.mark.parametrize("help", [5, 5.7, True, [5, 7, -4]])
def test_bad_help(help):
    """Test that TypeError is raised when 'help' is not a string."""

    with pytest.raises(TypeError):
        b = Boolean(help=help)

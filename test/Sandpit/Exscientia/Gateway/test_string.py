from BioSimSpace.Gateway import String

import pytest

# Tests for String requirements.

def test_no_arguments():
    """Test that calling constructor with no arguments will raise a ValueError."""

    with pytest.raises(ValueError):
        s = String()

@pytest.mark.parametrize("value", ["foo", "bar", "baz"])
def test_value(value):
    """Test whether object is initialised correctly and value is set."""

    s = String(help="Help!")
    s.setValue(value)

    assert s.getValue() == value
    assert s.getHelp() == "Help!"
    assert s.isOptional() == False
    assert s.getDefault() == None
    assert s.getMin() == None
    assert s.getMax() == None
    assert s.getAllowedValues() is None
    assert s.isMulti() == False
    assert s.getArgType() == str

@pytest.mark.parametrize("value", [1.5, 7, True])
def test_bad_value(value):
    """Check that TypeError is raised when "value" is the wrong type."""

    s = String(help='no')
    with pytest.raises(TypeError):
        s.setValue(value)

@pytest.mark.parametrize("default", ["foo", "bar", "baz"])
def test_default(default):
    """Check that the default value is set correctly."""

    s = String(help="Help!", default=default)

    assert s.getHelp() == "Help!"
    assert s.isOptional() == True
    assert s.getDefault() == default
    assert s.getMin() == None
    assert s.getMax() == None
    assert s.getAllowedValues() is None
    assert s.isMulti() == False
    assert s.getArgType() == str

    # Make sure that the argument is now optional.
    assert s.isOptional() is True

@pytest.mark.parametrize("value", ["foo", "bar"])
def test_valid_allowed(value):
    """Test that all is okay if the value is in the allowed set."""

    allowed = ["foo", "bar", "baz", "spam", "ham", "eggs"]
    s = String(help="Help!", allowed=allowed)
    s.setValue(value)
    assert s.getValue()== value
    assert s.getHelp() == "Help!"
    assert s.isOptional() == False
    assert s.getDefault() == None
    assert s.getMin() == None
    assert s.getMax() == None
    assert s.getAllowedValues() == allowed
    assert s.isMulti() == False
    assert s.getArgType() == str

@pytest.mark.parametrize("value", ["lorem", "ipsum"])
def test_invalid_allowed(value):
    """Check that ValueError is raised if the value is not in the allowed set."""

    allowed = ["foo", "bar", "baz", "spam", "ham", "eggs"]
    s = String(help="Help!", allowed=allowed)
    with pytest.raises(ValueError):
        s.setValue(value)

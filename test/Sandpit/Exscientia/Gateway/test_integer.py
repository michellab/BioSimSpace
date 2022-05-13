from BioSimSpace.Gateway import Integer

import pytest

# Tests for Integer requirements.

def test_no_arguments():
    """Test that calling constructor with no arguments will raise a ValueError."""

    with pytest.raises(ValueError):
        i = Integer()

@pytest.mark.parametrize("value", [-10, 0, 5])
def test_value(value):
    """Test whether object is initialised correctly and value is set."""

    i = Integer(help="Help!")
    i.setValue(value)

    assert i.getValue() == value
    assert i.getHelp() == "Help!"
    assert i.isOptional() == False
    assert i.getDefault() == None
    assert i.getMin() == None
    assert i.getMax() == None
    assert i.getAllowedValues() is None
    assert i.isMulti() == False
    assert i.getArgType() == int

@pytest.mark.parametrize("value", [1.5, '1.5', True, 'True'])
def test_bad_value(value):
    """Check that TypeError is raised when "value" is the wrong type."""

    i = Integer(help='no')
    with pytest.raises(TypeError):
        i.setValue(value)

@pytest.mark.parametrize("default", [-10, 0, 5])
def test_default(default):
    """Check that the default value is set correctly."""

    i = Integer(help="Help!", default=default)

    assert i.getHelp() == "Help!"
    assert i.isOptional() == True
    assert i.getDefault() == default
    assert i.getMin() == None
    assert i.getMax() == None
    assert i.getAllowedValues() is None
    assert i.isMulti() == False
    assert i.getArgType() == int

    # Make sure that the argument is now optional.
    assert i.isOptional() is True

@pytest.mark.parametrize("value", [-1, 0, 6])
def test_valid_min_max(value):
    """Check that all is okay when value is within min/max."""

    i = Integer(help= "Help!", minimum=-2, maximum=7)
    i.setValue(value)
    assert i.getValue()==value
    assert i.getHelp() == "Help!"
    assert i.isOptional() == False
    assert i.getDefault() == None
    assert i.getMin() == -2
    assert i.getMax() == 7
    assert i.getAllowedValues() is None
    assert i.isMulti() == False
    assert i.getArgType() == int

@pytest.mark.parametrize("value", [-7, 9])
def test_invalid_min_max(value):
    """Check that ValueError is raised if value is outside of min/max range."""

    i = Integer(help= "Help!", minimum=-2, maximum=7)
    with pytest.raises(ValueError):
        i.setValue(value)

@pytest.mark.parametrize("default", [-1, 0, 6])
def test_valid_min_max_default_valid(default):
    """Check that default is handled correctly when min/max are set."""

    i = Integer(help="Help!", minimum=-2, maximum=7, default=default)
    assert i.getValue()==None
    assert i.getHelp() == "Help!"
    assert i.isOptional() == True
    assert i.getDefault() == default
    assert i.getMin() == -2
    assert i.getMax() == 7
    assert i.getAllowedValues() is None
    assert i.isMulti() == False
    assert i.getArgType() == int

def test_invalid_min_max():
    """Check that ValueError is raised when min and max are incompatible."""

    with pytest.raises(ValueError):
        i = Integer(help="Help!", minimum=10, maximum=0)

@pytest.mark.parametrize("default", [-7, 9])
def test_invalid_min_max_default(default):
    """Check that ValueError is raised when default value is outside of min/max range."""

    with pytest.raises(ValueError):
        i = Integer(help="Help!", minimum=-2, maximum=7, default=default)

@pytest.mark.parametrize("value", [-7, 9])
def test_valid_allowed(value):
    """Test that all is okay if the value is in the allowed set."""

    allowed = [-4, -5, -7, 20, 9, 7]
    i = Integer(help="Help!", allowed=allowed)
    i.setValue(value)
    assert i.getValue()== value
    assert i.getHelp() == "Help!"
    assert i.isOptional() == False
    assert i.getDefault() == None
    assert i.getMin() == None
    assert i.getMax() == None
    assert i.getAllowedValues() == allowed
    assert i.isMulti() == False
    assert i.getArgType() == int

@pytest.mark.parametrize("value", [-7, 9])
def test_invalid_allowed(value):
    """Check that ValueError is raised if the value is not in the allowed set."""

    allowed = [-4, -5, 20, 7]
    i = Integer(help="Help!", allowed=allowed)
    with pytest.raises(ValueError):
        i.setValue(value)

def test_mixed_constraints():
    """Check that ValueError is raised if min/max and allowed are set."""

    with pytest.raises(ValueError):
        i = Integer(help="Help!", minimum=0, maximum=10, allowed=[1, 2, 3])

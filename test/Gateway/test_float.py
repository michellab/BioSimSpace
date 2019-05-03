from BioSimSpace.Gateway import Float

import pytest

# Tests for Float requirements.

def test_no_arguments():
    """Test that calling constructor with no arguments will raise a ValueError."""

    with pytest.raises(ValueError):
        f = Float()

@pytest.mark.parametrize("value", [-10.0, 0.0, 5.0])
def test_value(value):
    """Test whether object is initialised correctly and value is set."""

    # Create Float object with minimium required arguments.
    f = Float(help="Help!")

    # Set the value of the float.
    f.setValue(value)

    # Assert that all arguments are parsed correctly and the value is set.
    assert f.getValue() == value
    assert f.getHelp() == "Help!"
    assert f.isOptional() == False
    assert f.getDefault() == None
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() is None
    assert f.isMulti() == False
    assert f.getArgType() == float

@pytest.mark.parametrize("default", [-10.0, 0.0, 5.0])
def test_default(default):
    """Check that default value is set correctly."""

    f = Float(help="Help!", default=default)

    assert f.getHelp() == "Help!"
    assert f.isOptional() == True
    assert f.getDefault() == default
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() is None
    assert f.isMulti() == False
    assert f.getArgType() == float

    # Make sure that the argument is now optional.
    assert f.isOptional() is True

@pytest.mark.parametrize("value", ['1.5', True])
def test_bad_value(value):
    """Test that TypeError is raised when "value" is not a float or int."""

    f = Float(help="Help!")
    with pytest.raises(TypeError):
        f.setValue(value)

@pytest.mark.parametrize("value", [-1.0, 0.0, 6.0])
def test_valid_min_max(value):
    """Check that all is okay when the value is within min/max."""

    f = Float(help= "Help!", minimum=-2.0, maximum=7.9)
    f.setValue(value)
    assert f.getValue()==value
    assert f.getHelp() == "Help!"
    assert f.isOptional() == False
    assert f.getDefault() == None
    assert f.getMin() == -2.0
    assert f.getMax() == 7.9
    assert f.getAllowedValues() is None
    assert f.isMulti() == False
    assert f.getArgType() == float

@pytest.mark.parametrize("value", [-7.0, 9.0])
def test_invalid_min_max(value):
    """Check that ValueError is raised if value is outside of min/max range."""

    f = Float(help= "Help!", minimum=-2.0, maximum=7.9)
    with pytest.raises(ValueError):
        f.setValue(value)

@pytest.mark.parametrize("default", [-1.0, 0.0, 6.0])
def test_valid_min_max_default_valid(default):
    """Check that default is handled correctly when min/max are set."""

    f = Float(help="Help!", minimum=-2.0, maximum=7.9, default=default)
    assert f.getValue()==None
    assert f.getHelp() == "Help!"
    assert f.isOptional() == True
    assert f.getDefault() == default
    assert f.getMin() == -2.0
    assert f.getMax() == 7.9
    assert f.getAllowedValues() is None
    assert f.isMulti() == False
    assert f.getArgType() == float

@pytest.mark.parametrize("default", [-7.0, 9.0])

def test_invalid_min_max_default(default):
    """Check that ValueError is raised when default value is outside of min/max range."""

    with pytest.raises(ValueError):
        f = Float(help="Help!", minimum=-2.0, maximum=7.0, default=default)

def test_invalid_min_max():
    """Check that ValueError is raised when min and max are incompatible."""

    with pytest.raises(ValueError):
        f = Float(help="Help!", minimum=10.0, maximum=0.0)

@pytest.mark.parametrize("value", [-7.0, 9.0])
def test_valid_allowed(value):
    """Test that all is okay if the value is in the allowed set."""

    allowed = [-4.5, -5.1, -7.0, 20.0, 9.0, 7.0]
    f = Float(help="Help!", allowed=allowed)
    f.setValue(value)
    assert f.getValue()== value
    assert f.getHelp() == "Help!"
    assert f.isOptional() == False
    assert f.getDefault() == None
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() == allowed
    assert f.isMulti() == False
    assert f.getArgType() == float

@pytest.mark.parametrize("value", [-7.0, 9.0])
def test_invalid_allowed(value):
    """Check that ValueError is raised if the value is not in the allowed set."""

    allowed = [-4.0, -5.0, 20.0, 7.0]
    f = Float(help="Help!", allowed=allowed)
    with pytest.raises(ValueError):
        f.setValue(value)

def test_mixed_constraints():
    """Check that ValueError is raised if min/max and allowed are set."""

    with pytest.raises(ValueError):
        f = Float(help="Help!", minimum=0.0, maximum=10.0, allowed=[1.0, 2.0, 3.0])

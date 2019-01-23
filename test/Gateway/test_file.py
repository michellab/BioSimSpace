from BioSimSpace.Gateway import File

import pytest

# Tests for File requirements.

def test_no_arguments():
    """Test that calling constructor with no arguments will raise a ValueError."""

    with pytest.raises(ValueError):
        f = File()

@pytest.mark.parametrize("value", ["test/io/amber/ala/ala.crd", "test/io/amber/ala/ala.top"])
def test_value(value):
    """Test whether object is initialised correctly and value is set."""

    f = File(help="Help!")
    f.setValue(value)

    assert f.getValue() == value
    assert f.getHelp() == "Help!"
    assert f.isOptional() == False
    assert f.getDefault() == None
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() is None
    assert f.isMulti() == False
    assert f.getArgType() == str

@pytest.mark.parametrize("value", [1.5, 7, True])
def test_bad_value(value):
    """Check that TypeError is raised when "value" is the wrong type."""

    f = File(help="Help!")
    with pytest.raises(TypeError):
        f.setValue(value)

def test_missing_file():
    """Check that IOError is raised when the file doesn't exist."""

    with pytest.raises(IOError):
        f = File(help="Help!")
        f.setValue("missing.txt")

@pytest.mark.parametrize("optional", [True, False])
def test_optional(optional):
    """Check that the optional flag is set correctly."""

    f = File(help="Help!", optional=optional)

    assert f.getHelp() == "Help!"
    assert f.isOptional() == optional
    assert f.getDefault() == None
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() is None
    assert f.isMulti() == False
    assert f.getArgType() == str

from BioSimSpace.Gateway import FileSet

import pytest

# Tests for FileSet requirements.

def test_no_arguments():
    """Test that calling constructor with no arguments will raise a ValueError."""

    with pytest.raises(ValueError):
        f = FileSet()

@pytest.mark.parametrize("value", [["test/io/amber/ala/ala.crd", "test/io/amber/ala/ala.top"],
    ["test/io/namd/alanin/alanin.pdb", "test/io/namd/alanin/alanin.psf", "test/io/namd/alanin/alanin.params"]])
def test_value(value):
    """Test whether object is initialised correctly and value is set."""

    f = FileSet(help="Help!")
    f.setValue(value)

    assert f.getValue() == value
    assert f.getHelp() == "Help!"
    assert f.isOptional() == False
    assert f.getDefault() == None
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() is None
    assert f.isMulti() == True
    assert f.getArgType() == str

@pytest.mark.parametrize("value", [[1.5, 7, True], [33.0, False]])
def test_bad_value(value):
    """Check that TypeError is raised when "value" is the wrong type."""

    f = FileSet(help="Help!")
    with pytest.raises(TypeError):
        f.setValue(value)

def test_missing_files():
    """Check that IOError is raised when a file doesn't exist."""

    # All files missing.
    with pytest.raises(IOError):
        f = FileSet(help="Help!")
        f.setValue(["missing1.txt", "missing2.txt"])

    # One file missing.
    with pytest.raises(IOError):
        f = FileSet(help="Help!")
        f.setValue(["test/io/amber/ala/ala.crd", "missing2.txt"])

@pytest.mark.parametrize("optional", [True, False])
def test_optional(optional):
    """Check that the optional flag is set correctly."""

    f = FileSet(help="Help!", optional=optional)

    assert f.getHelp() == "Help!"
    assert f.isOptional() == optional
    assert f.getDefault() == None
    assert f.getMin() == None
    assert f.getMax() == None
    assert f.getAllowedValues() is None
    assert f.isMulti() == True
    assert f.getArgType() == str

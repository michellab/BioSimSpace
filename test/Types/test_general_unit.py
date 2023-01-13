import BioSimSpace.Types as Types
import BioSimSpace.Units as Units

import pytest


@pytest.mark.parametrize(
    "string, dimensions",
    [
        ("kilo Cal oriEs per Mole / angstrom **2", (0, 0, 0, 1, -1, 0, -2)),
        ("k Cal_per  _mOl / nm^2", (0, 0, 0, 1, -1, 0, -2)),
        ("kj p  eR  moles / pico METERs2", (0, 0, 0, 1, -1, 0, -2)),
        ("coul oMbs / secs * ATm os phereS", (0, 1, -1, 1, 0, 0, -3)),
        ("pm**3 * rads * de grEE", (2, 0, 3, 0, 0, 0, 0)),
    ],
)
def test_supported_units(string, dimensions):
    """Test that we can create GeneralUnit objects with the correct dimensions
    by evaluating strings as unit based algebraic expressions.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that the dimensions match.
    assert general_unit.dimensions() == dimensions


@pytest.mark.parametrize(
    "string, matching_type",
    [
        ("radian * degree**2 / radian^2", Types.Angle),
        ("angstrom**3 / nanometer", Types.Area),
        ("coulombs * angstrom**-2 * nanometer**2", Types.Charge),
        ("kcal_per_mol / angstrom**2 * nanometer**2", Types.Energy),
        ("angstrom**3 * nanometer^-1 / picometer", Types.Length),
        ("bar * kJ_per_mol**2 / (kcal_per_mol * kJ_per_mol)", Types.Pressure),
        ("coulomb * kelvin^-3 * celsius**2 * kelvin^2 / e_charge", Types.Temperature),
        ("nanoseconds^3 * kelvin^-3 * celsius**3 / milliseconds**2", Types.Time),
        ("angstroms cubed * atm^-3 * bar**3", Types.Volume),
    ],
)
def test_type_conversion(string, matching_type):
    """Test that GeneralUnit objects can be converted to a type with matching
    dimensions.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that the types match.
    assert type(general_unit) is matching_type


@pytest.mark.parametrize(
    "string, default_unit",
    [
        ("degree", Units.Angle.radian),
        ("meters2", Units.Area.angstrom2),
        ("coulombs", Units.Charge.electron_charge),
        ("kJ_per_mol", Units.Energy.kcal_per_mol),
        ("nanometer", Units.Length.angstrom),
        ("bar", Units.Pressure.atm),
        ("fahrenheit", Units.Temperature.kelvin),
        ("days", Units.Time.nanosecond),
        ("picometers**3", Units.Volume.angstrom3),
    ],
)
def test_default_conversion(string, default_unit):
    """Test that GeneralUnit objects are always converted to the default
    unit for that type.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that units match.
    assert general_unit.unit() == default_unit.unit()


@pytest.mark.parametrize(
    "unit_type",
    [
        Units.Angle.radian,
        Units.Area.angstrom2,
        Units.Charge.electron_charge,
        Units.Energy.kcal_per_mol,
        Units.Length.angstrom,
        Units.Pressure.atm,
        Units.Temperature.kelvin,
        Units.Time.nanosecond,
        Units.Volume.angstrom3,
    ],
)
def test_pos_pow(unit_type):
    """Test that unit-based types can be raised to positive powers."""

    # Store the dimensions associated with the original type.
    old_dimensions = unit_type.dimensions()

    # Square the unit-based type.
    unit_type = unit_type**2

    # Store the new dimensions.
    new_dimensions = unit_type.dimensions()

    # Each dimension entry should be twice the old value.
    for d0, d1 in zip(old_dimensions, new_dimensions):
        assert d1 == 2 * d0


@pytest.mark.parametrize(
    "unit_type",
    [
        Units.Angle.radian,
        Units.Area.angstrom2,
        Units.Charge.electron_charge,
        Units.Energy.kcal_per_mol,
        Units.Length.angstrom,
        Units.Pressure.atm,
        Units.Temperature.kelvin,
        Units.Time.nanosecond,
        Units.Volume.angstrom3,
    ],
)
def test_neg_pow(unit_type):
    """Test that unit-based types can be raised to negative powers."""

    # Store the dimensions associated with the original type.
    old_dimensions = unit_type.dimensions()

    # Invert the unit-based type.
    unit_type = unit_type**-1

    # Store the new dimensions.
    new_dimensions = unit_type.dimensions()

    # Each dimension entry should be the inverse of the old value.
    for d0, d1 in zip(old_dimensions, new_dimensions):
        assert d1 == -d0


@pytest.mark.parametrize(
    "string",
    [
        "degree",
        "meters2",
        "coulombs",
        "kJ_per_mol",
        "nanometer",
        "bar",
        "fahrenheit",
        "days",
        "picometers**3",
    ],
)
def test_dimensionless(string):
    """Test that GeneralUnit objects convert to dimensionless float values
    when divided by themself.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Check that we get back a float when divided by itself.
    assert isinstance(general_unit / general_unit, float)


def test_dimensionless_value():
    """Check that conversion to a dimensionless unit preserves the value
    of the unit conversion.
    """

    value = (Units.Energy.kcal_per_mol / Units.Length.angstrom**2) / (
        Units.Energy.kj_per_mol / Units.Length.nanometer**2
    )

    assert value == pytest.approx(418.4)

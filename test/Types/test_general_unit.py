import BioSimSpace.Types as Types
import BioSimSpace.Units as Units

import pytest

@pytest.mark.parametrize("string, dimensions",
        [("kilo Cal oriEs per Mole / angstrom **2", (0, 0, 0, 1, -1, 0, -2)),
         ("k Cal_per  _mOl / nm^2", (0, 0, 0, 1, -1, 0, -2)),
         ("kj p  eR  moles / pico METERs2", (0, 0, 0, 1, -1, 0, -2)),
         ("coul oMbs / secs * ATm os phereS", (0, 1, -1, 1, 0, 0, -3)),
         ("pm**3 * rads * de grEE", (2, 0, 3, 0, 0, 0, 0)),
        ])
def test_supported_units(string, dimensions):
    """Test that we can create GeneralUnit objects with the correct dimensions
       by evaluating strings as unit based algebraic expressions.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that the dimensions match.
    assert general_unit.dimensions() == dimensions

@pytest.mark.parametrize("string, matching_type",
        [("radian * degree**2 / radian^2", Types.Angle),
         ("angstrom**3 / nanometer", Types.Area),
         ("coulombs * angstrom**-2 * nanometer**2", Types.Charge),
         ("kcal_per_mol / angstrom**2 * nanometer**2", Types.Energy),
         ("angstrom**3 * nanometer^-1 / picometer", Types.Length),
         ("bar * kJ_per_mol**2 / (kcal_per_mol * kJ_per_mol)", Types.Pressure),
         ("coulomb * kelvin^-3 * celsius**2 * kelvin^2 / e_charge", Types.Temperature),
         ("nanoseconds^3 * kelvin^-3 * celsius**3 / milliseconds**2", Types.Time),
         ("angstroms cubed * atm^-3 * bar**3", Types.Volume)
        ])
def test_type_conversion(string, matching_type):
    """Test that GeneralUnit objects can be converted to a type with matching
       dimensions.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that the types match.
    assert type(general_unit) is matching_type

@pytest.mark.parametrize("string, default_unit",
        [("degree", Units.Angle.radian),
         ("meters2", Units.Area.angstrom2),
         ("coulombs", Units.Charge.electron_charge),
         ("kJ_per_mol", Units.Energy.kcal_per_mol),
         ("nanometer", Units.Length.angstrom),
         ("bar", Units.Pressure.atm),
         ("fahrenheit", Units.Temperature.kelvin),
         ("days", Units.Time.nanosecond),
         ("picometers**3", Units.Volume.angstrom3)
        ])
def test_default_conversion(string, default_unit):
    """Test that GeneralUnit objects are always converted to the default
       unit for that type.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that units match.
    assert general_unit.unit() == default_unit.unit()

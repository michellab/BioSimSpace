import BioSimSpace.Types as Types

import pytest
import random

# Create a list of the Type classes.
types = []
for var in dir(Types):
    if var[0] == "_":
        continue
    elif var[0].upper() == var[0]:
        # Exclude non unit based types.
        if var != "Vector" and var != "Coordinate":
            types.append(getattr(Types, var))

def insert_space(string, index):
    return string[:index] + " " + string[index:]

def lower_case(string, index):
    x = len(string)
    if index == 0:
        return string[0].lower() + string[1:]
    elif index == x - 1:
        return string[:x-1] + string[x-1].lower()
    else:
        return string[:index] + string[index].lower() + string[index+1:]

@pytest.mark.parametrize("Type", types)
def test_supported_units(Type):
    """Test that it's possible to instatiate objects for all supported units."""

    # Store a list of the supported units.
    units = list(Type._supported_units.keys())

    # Loop over all units.
    for unit in units:

        # Attempt to instantiate an object.
        my_type = Type(1.0, unit)

        # Swap the case of random characters.
        for x in range(0, 5):
            unit = lower_case(unit, random.randint(0, len(unit)-1))

        # Add random whitespace.
        for x in range(0, 5):
            unit = insert_space(unit, random.randint(0, len(unit)-1))

        # Attempt to instantiate an object.
        my_type = Type(1.0, unit)

@pytest.mark.parametrize("Type", types)
def test_abbreviations(Type):
    """Test that all unit abbreviations work."""

    # Store a list of the unit abbreviations.
    abbreviations = list(Type._abbreviations.keys())

    # Loop over all abbreviations.
    for abbrev in abbreviations:

        # Attempt to instantiate an object.
        my_type = Type(1.0, abbrev)

        # Swap the case of random characters.
        for x in range(0, 5):
            abbrev = lower_case(abbrev, random.randint(0, len(abbrev)-1))

        # Add random whitespace.
        for x in range(0, 5):
            unit = insert_space(abbrev, random.randint(0, len(abbrev)-1))

        # Attempt to instantiate an object.
        my_type = Type(1.0, abbrev)

@pytest.mark.parametrize("Type", types)
def test_round_trip(Type):
    """Test that round-trip unit conversion works."""

    # Store a list of the supported units.
    units = list(Type._supported_units.keys())

    # Loop over all units.
    for index, unit in enumerate(units):

        # Create the initial object.
        my_type = Type(1.0, unit)

        # Convert to the next unit.
        for x in range(1, len(units)+1):
            my_type = my_type._convert_to(units[index-x])

        # Make sure the unit and magnitude are correct.
        assert my_type.magnitude() == pytest.approx(1.0)
        assert unit == my_type.unit()

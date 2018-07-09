# BioSimSpace.Units

This package provides a set of common physical units required by BioSimSpace.

## `Length`

* `meter`
* `centimeter`
* `millimeter`
* `nanometer`
* `angstrom`
* `picometer`

## `Area`

* `meter2`
* `nanometer2`
* `angstrom2`
* `picometer2`

## `Volume`

* `meter3`
* `nanometer3`
* `angstrom3`
* `picometer3`

## `Energy`

* `kcal_per_mol`
* `kj_per_mol`
* `kt`

## `Pressure`

* `atm`
* `bar`

## `Temperature`

* `celsius`
* `fahrenheit`
* `kelvin`

## `Time`

* `day`
* `hour`
* `minute`
* `second`
* `millisecond`
* `nanosecond`
* `picosecond`
* `femtosecond`

Units are intrinsically linked to [BioSimSpace.Types](../Types), i.e. units can
be used to generate an object of the type associated with the physical unit.

For example:

```python
import BioSimSpace as BSS

# This will return a temperature object.
temp = 300 * BSS.Units.Temperature.kelvin

# This will return a length object.
length = 10 * BSS.Units.Length.angstrom
```

######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""Functionality for equilibration protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Equilibration"]

import math as _math
import warnings as _warnings

from .. import Types as _Types
from .. import Units as _Units

from ._protocol import Protocol as _Protocol


class Equilibration(_Protocol):
    """A class for storing equilibration protocols."""

    # Supported restraint keywords.
    _restraints = ["backbone", "heavy", "all", "none"]

    def __init__(
        self,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(0.2, "nanoseconds"),
        temperature_start=_Types.Temperature(300, "kelvin"),
        temperature_end=_Types.Temperature(300, "kelvin"),
        temperature=None,
        pressure=None,
        report_interval=100,
        restart_interval=500,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
    ):
        """
        Constructor.

        Parameters
        ----------

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature_start : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.

        temperature_end : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The equilibration temperature. This takes precedence of over
            the other temperatures, i.e. to run at fixed temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. If this argument is omitted then the simulation
            is run using the NVT ensemble.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory
            frames are saved. (In integration steps.)

        restraint : str, [int]
            The type of restraint to perform. This should be one of the
            following options:
                "backbone"
                     Protein backbone atoms. The matching is done by a name
                     template, so is unreliable on conversion between
                     molecular file formats.
                "heavy"
                     All non-hydrogen atoms that aren't part of water
                     molecules or free ions.
                "all"
                     All atoms that aren't part of water molecules or free
                     ions.
            Alternatively, the user can pass a list of atom indices for
            more fine-grained control. If None, then no restraints are used.

        force_constant : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
            The force constant for the restraint potential. If a 'float' is
            passed, then default units of 'kcal_per_mol / angstrom**2' will
            be used.
        """

        # Call the base class constructor.
        super().__init__()

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the running time.
        self.setRunTime(runtime)

        # Constant temperature equilibration.
        if temperature is not None:
            self.setStartTemperature(temperature)
            self.setEndTemperature(temperature)
            self._is_const_temp = True

        # Heating / cooling simulation.
        else:
            self._is_const_temp = False

            # Set the start temperature.
            self.setStartTemperature(temperature_start)

            # Set the final temperature.
            self.setEndTemperature(temperature_end)

            # Constant temperature simulation.
            if self._temperature_start == self._temperature_end:
                self._is_const_temp = True

        # Constant pressure simulation.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        # Set the report interval.
        self.setReportInterval(report_interval)

        # Set the restart interval.
        self.setRestartInterval(restart_interval)

        # Set the restraint.
        if restraint is not None:
            self.setRestraint(restraint)
        else:
            self._restraint = None

        # Set the force constant.
        self.setForceConstant(force_constant)

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return (
                "<BioSimSpace.Protocol.Equilibration: timestep=%s, runtime=%s, "
                "temperature_start=%s, temperature_end=%s, pressure=%s, "
                "report_interval=%d, restart_interval=%d,restraint=%r, "
                "force_constant=%3.2f kcal_per_mol/angstrom**2>"
            ) % (
                self._timestep,
                self._runtime,
                self._temperature_start,
                self._temperature_end,
                self._pressure,
                self._report_interval,
                self._restart_interval,
                self._restraint,
                self._force_constant.value(),
            )

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return (
                "BioSimSpace.Protocol.Equilibration(timestep=%s, runtime=%s, "
                "temperature_start=%s, temperature_end=%s, pressure=%s, "
                "report_interval=%d, restart_interval=%d, restraint=%r, force_constant=%3.2f)"
            ) % (
                self._timestep,
                self._runtime,
                self._temperature_start,
                self._temperature_end,
                self._pressure,
                self._report_interval,
                self._restart_interval,
                self._restraint,
                self._force_constant.value(),
            )

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
        """
        Set the time step.

        Parameters
        ----------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        if isinstance(timestep, _Types.Time):
            self._timestep = timestep
        else:
            raise TypeError("'timestep' must be of type 'BioSimSpace.Types.Time'")

    def getRunTime(self):
        """
        Return the running time.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        return self._runtime

    def setRunTime(self, runtime):
        """
        Set the running time.

        Parameters
        ----------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        if isinstance(runtime, _Types.Time):
            self._runtime = runtime
        else:
            raise TypeError("'runtime' must be of type 'BioSimSpace.Types.Time'")

    def getStartTemperature(self):
        """
        Return the starting temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.
        """
        return self._temperature_start

    def setStartTemperature(self, temperature):
        """
        Set the starting temperature.

        Parameters
        ----------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.
        """

        if isinstance(temperature, _Types.Temperature):
            if _math.isclose(temperature.kelvin().value(), 0, rel_tol=1e-6):
                temperature._value = 0.01
            self._temperature_start = temperature
        else:
            raise TypeError(
                "'temperature_start' must be of type 'BioSimSpace.Types.Temperature'"
            )

    def getEndTemperature(self):
        """
        Return the final temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.
        """
        return self._temperature_end

    def setEndTemperature(self, temperature):
        """
        Set the final temperature.

        Parameters
        ----------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.
        """
        if isinstance(temperature, _Types.Temperature):
            if _math.isclose(temperature.kelvin().value(), 0, rel_tol=1e-6):
                temperature._value = 0.01
            self._temperature_end = temperature
        else:
            raise TypeError(
                "'temperature_end' must be of type 'BioSimSpace.Types.Temperature'"
            )

    def getPressure(self):
        """
        Return the pressure.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """
        Set the pressure.

        Parameters
        ----------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        if isinstance(pressure, _Types.Pressure):
            self._pressure = pressure
        else:
            raise TypeError("'pressure' must be of type 'BioSimSpace.Types.Pressure'")

    def getReportInterval(self):
        """
        Return the interval between reporting statistics. (In integration steps.).

        Returns
        -------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        return self._report_interval

    def setReportInterval(self, report_interval):
        """
        Set the interval at which statistics are reported. (In integration steps.).

        Parameters
        ----------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        if not type(report_interval) is int:
            raise TypeError("'report_interval' must be of type 'int'")

        if report_interval <= 0:
            _warnings.warn("'report_interval' must be positive. Using default (100).")
            report_interval = 100

        self._report_interval = report_interval

    def getRestartInterval(self):
        """
        Return the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Returns
        -------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        return self._restart_interval

    def setRestartInterval(self, restart_interval):
        """
        Set the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Parameters
        ----------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval

    def getRestraint(self):
        """
        Return the type of restraint.

        Returns
        -------

        restraint : str, [int]
            The type of restraint.
        """
        return self._restraint

    def setRestraint(self, restraint):
        """
        Set the type of restraint.

        Parameters
        ----------

        restraint : str, [int]
            The type of restraint to perform. This should be one of the
            following options:
                "backbone"
                     Protein backbone atoms. The matching is done by a name
                     template, so is unreliable on conversion between
                     molecular file formats.
                "heavy"
                     All non-hydrogen atoms that aren't part of water
                     molecules or free ions.
                "all"
                     All atoms that aren't part of water molecules or free
                     ions.
            Alternatively, the user can pass a list of atom indices for
            more fine-grained control.
        """

        if isinstance(restraint, str):
            # Convert to lower case and strip whitespace.
            restraint = restraint.lower().replace(" ", "")
            if restraint not in self._restraints:
                raise ValueError(f"'restraint' must be one of: {self._restraints}")
            # Set to NoneType if equal to "none", since this makes checking
            # whether a restraint is set elsewhere much easier.
            if restraint == "none":
                restraint = None

        elif isinstance(restraint, (list, tuple)):
            if not all(type(x) is int for x in restraint):
                raise ValueError("'restraint' must be a list of 'int' types!")
            # Create a set to sort and ensure no duplicates, then convert back to a list.
            restraint = list(set(restraint))
            restraint.sort()

        else:
            raise TypeError(
                "'restraint' must be of type 'str', or a list of 'int' types."
            )

        self._restraint = restraint

    def getForceConstant(self):
        """
        Return the force constant for the restraint.

        Returns
        -------

        force_constant :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the restraint, in units of
            kcal_per_mol/angstrom**2.
        """
        return self._force_constant

    def setForceConstant(self, force_constant):
        """
        Set the type force constant for the restraint.

        Parameters
        ----------

        force_constant : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
        """

        # Convert int to float.
        if type(force_constant) is int:
            force_constant = float(force_constant)

        if isinstance(force_constant, float):
            # Use default units.
            force_constant *= _Units.Energy.kcal_per_mol / _Units.Area.angstrom2

        elif isinstance(force_constant, _Types._GeneralUnit):
            # Validate the dimensions.
            if force_constant.dimensions() != (0, 0, 0, 1, -1, 0, -2):
                raise ValueError(
                    "'force_constant' has invalid dimensions! "
                    f"Expected dimensions are 'M Q-1 T-2', found '{force_constant.unit()}'"
                )

        else:
            raise TypeError(
                "'force_constant' must be of type 'BioSimSpace.Types._GeneralUnit', or 'float'."
            )

        self._force_constant = force_constant

    def isConstantTemp(self):
        """
        Return whether the protocol has a constant temperature.

        Returns
        -------

        is_const_temp : bool
            Whether the temperature is fixed.
        """
        return self._temperature_start == self._temperature_end

    @classmethod
    def restraints(cls):
        """
        Return a list of the supported restraint keywords.

        Returns
        -------

        restraints : [str]
            A list of the supported restraint keywords.
        """
        return cls._restraints.copy()

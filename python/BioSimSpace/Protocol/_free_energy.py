######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
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

"""
Functionality for free energy protocols.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["FreeEnergy"]

from BioSimSpace import Types as _Types

from ._free_energy_mixin import _FreeEnergyMixin
from ._production import Production as _Producton


class FreeEnergy(_Producton, _FreeEnergyMixin):
    """A class for storing free energy production protocols."""

    def __init__(self,
                 lam=0.0,
                 lam_vals=None,
                 min_lam=0.0,
                 max_lam=1.0,
                 num_lam=11,
                 timestep=_Types.Time(2, "femtosecond"),
                 runtime=_Types.Time(4, "nanosecond"),
                 temperature=_Types.Temperature(300, "kelvin"),
                 pressure=_Types.Pressure(1, "atmosphere"),
                 report_interval=200000,
                 restart_interval=20000,
                 perturbation_type="full"
                ):
        """Constructor.

           Parameters
           ----------

           lam : float
               The perturbation parameter: [0.0, 1.0]

           lam_vals : [float]
               The list of lambda parameters.

           min_lam : float
               The minimum lambda value.

           max_lam : float
               The maximum lambda value.

           num_lam : int
               The number of lambda values.

           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration timestep.

           runtime : :class:`Time <BioSimSpace.Types.Time>`
               The running time.

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The temperature.

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure. Pass pressure=None to use the NVT ensemble.

           report_interval : int
               The frequency at which statistics are recorded. (In integration steps.)

           restart_interval : int
               The frequency at which restart configurations and trajectory

           perturbation_type : str
               The type of perturbation to perform. Options are:
                "full" : A full perturbation of all terms (default option).
                "discharge_soft" : Perturb all discharging soft atom charge terms (i.e. value->0.0).
                "vanish_soft" : Perturb all vanishing soft atom LJ terms (i.e. value->0.0).
                "flip" : Perturb all hard atom terms as well as bonds/angles.
                "grow_soft" : Perturb all growing soft atom LJ terms (i.e. 0.0->value).
                "charge_soft" : Perturb all charging soft atom LJ terms (i.e. 0.0->value).

                Currently perturubation_type != "full" is only supported by
                BioSimSpace.Process.Somd.
        """

        # Call the base class constructors.
        _Producton.__init__(self,
                            timestep=timestep,
                            runtime=runtime,
                            temperature=temperature,
                            pressure=pressure,
                            report_interval=report_interval,
                            restart_interval=restart_interval)

        _FreeEnergyMixin.__init__(self,
                                  lam=lam,
                                  lam_vals=lam_vals,
                                  min_lam=min_lam,
                                  max_lam=max_lam,
                                  num_lam=num_lam,
                                  perturbation_type=perturbation_type)

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("<BioSimSpace.Protocol.FreeEnergy: lam=%5.4f, lam_vals=%r, timestep=%s, "
                    "runtime=%s, temperature=%s, pressure=%s, report_interval=%d, restart_interval=%d>"
                   ) % (self._lambda, self._lambda_vals, self._timestep, self._runtime,
                        self._temperature, self._pressure, self._report_interval, self._restart_interval)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("BioSimSpace.Protocol.FreeEnergy(lam=%5.4f, lam_vals=%r, timestep=%s, "
                    "runtime=%s, temperature=%s, pressure=%s, report_interval=%d, restart_interval=%d)"
                   ) % (self._lambda, self._lambda_vals, self._timestep, self._runtime,
                        self._temperature, self._pressure, self._report_interval, self._restart_interval)

    def getTimeStep(self):
        """Return the time step.

           Returns
           -------

           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
        """Set the time step.

           Parameters
           ----------

           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration time step.
        """
        if type(timestep) is _Types.Time:
            self._timestep = timestep
        else:
            raise TypeError("'timestep' must be of type 'BioSimSpace.Types.Time'")

    def getRunTime(self):
        """Return the running time.

           Returns
           -------

           runtime : :class:`Time <BioSimSpace.Types.Time>`
               The simulation run time.
        """
        return self._runtime

    def setRunTime(self, runtime):
        """Set the running time.

           Parameters
           ----------

           runtime : :class:`Time <BioSimSpace.Types.Time>`
               The simulation run time.
        """
        if type(runtime) is _Types.Time:
            self._runtime = runtime
        else:
            raise TypeError("'runtime' must be of type 'BioSimSpace.Types.Time'")

    def getTemperature(self):
        """Return temperature.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The simulation temperature.
        """
        return self._temperature

    def setTemperature(self, temperature):
        """Set the temperature.

           Parameters
           ----------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The simulation temperature.
        """
        if type(temperature) is _Types.Temperature:
            self._temperature = temperature
        else:
            raise TypeError("'temperature' must be of type 'BioSimSpace.Types.Temperature'")

    def getPressure(self):
        """Return the pressure.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """Set the pressure.

           Parameters
           ----------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure.
        """
        if type(pressure) is _Types.Pressure:
            self._pressure = pressure
        else:
            raise TypeError("'pressure' must be of type 'BioSimSpace.Types.Pressure'")

    def getReportInterval(self):
        """Return the interval between reporting statistics. (In integration steps.)

           Returns
           -------

           report_interval : int
               The number of integration steps between reporting statistics.
        """
        return self._report_interval

    def setReportInterval(self, report_interval):
        """Set the interval at which statistics are reported. (In integration steps.)

           Parameters
           ----------

           report_interval : int
               The number of integration steps between reporting statistics.
        """
        if type(report_interval) is not int:
            raise TypeError("'report_interval' must be of type 'int'")

        if report_interval <= 0:
            _warnings.warn("'report_interval' must be positive. Using default (100).")
            report_interval = 100

        self._report_interval = report_interval

    def getRestartInterval(self):
        """Return the interval between saving restart confiugrations, and/or
           trajectory frames. (In integration steps.)

           Returns
           -------

           restart_interval : int
               The number of integration steps between saving restart
               configurations and/or trajectory frames.
        """
        return self._restart_interval

    def setRestartInterval(self, restart_interval):
        """Set the interval between saving restart confiugrations, and/or
           trajectory frames. (In integration steps.)

           Parameters
           ----------

           restart_interval : int
               The number of integration steps between saving restart
               configurations and/or trajectory frames.
        """
        if type(restart_interval) is not int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval

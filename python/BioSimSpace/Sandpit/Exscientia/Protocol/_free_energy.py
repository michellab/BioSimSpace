######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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

import math as _math
import warnings as _warnings

from .. import Types as _Types

from ._free_energy_mixin import _FreeEnergyMixin
from ._production import Production as _Production


class FreeEnergy(_Production, _FreeEnergyMixin):
    """A class for storing free energy production protocols."""

    def __init__(
        self,
        lam=0.0,
        lam_vals=None,
        min_lam=0.0,
        max_lam=1.0,
        num_lam=11,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(4, "nanosecond"),
        temperature=_Types.Temperature(300, "kelvin"),
        pressure=_Types.Pressure(1, "atmosphere"),
        report_interval=200,
        restart_interval=1000,
        first_step=0,
        restart=False,
        perturbation_type="full",
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

        first_step : int
            The initial time step (for restart simulations).

        restart : bool
            Whether this is a continuation of a previous simulation.

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
        _Production.__init__(
            self,
            timestep=timestep,
            runtime=runtime,
            temperature=temperature,
            pressure=pressure,
            report_interval=report_interval,
            restart_interval=restart_interval,
            first_step=first_step,
            restart=restart,
        )

        _FreeEnergyMixin.__init__(
            self,
            lam=lam,
            lam_vals=lam_vals,
            min_lam=min_lam,
            max_lam=max_lam,
            num_lam=num_lam,
            perturbation_type=perturbation_type,
        )

    def _get_parm(self):
        """Return a string representation of the parameters."""

        return ", ".join(
            [_Production._get_parm(self), _FreeEnergyMixin._get_parm(self)]
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol.FreeEnergy: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "BioSimSpace.Protocol.Custom"
        else:
            return f"BioSimSpace.Protocol.FreeEnergy({self._get_parm()})"

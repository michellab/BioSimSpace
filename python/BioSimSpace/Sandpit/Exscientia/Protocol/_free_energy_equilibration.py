__all__ = ["FreeEnergyEquilibration"]

from BioSimSpace import Types as _Types
from BioSimSpace import Units as _Units

from ._free_energy_mixin import _FreeEnergyMixin
from ._equilibration import Equilibration as _Equilibration


class FreeEnergyEquilibration(_Equilibration, _FreeEnergyMixin):
    """A class for storing free energy equilibration protocols."""

    def __init__(self,
                 lam=0.0,
                 lam_vals=None,
                 min_lam=0.0,
                 max_lam=1.0,
                 num_lam=11,
                 timestep=_Types.Time(2, "femtosecond"),
                 runtime=_Types.Time(0.2, "nanoseconds"),
                 temperature_start=_Types.Temperature(300, "kelvin"),
                 temperature_end=_Types.Temperature(300, "kelvin"),
                 temperature=None,
                 pressure=None,
                 report_interval=200,
                 restart_interval=1000,
                 restraint=None,
                 force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
                 restart=False,
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
        _Equilibration.__init__(self,
                                timestep=timestep,
                                runtime=runtime,
                                temperature_start=temperature_start,
                                temperature_end=temperature_end,
                                temperature=temperature,
                                pressure=pressure,
                                report_interval=report_interval,
                                restart_interval=restart_interval,
                                restraint=restraint,
                                force_constant=force_constant,
                                restart=restart)

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
            return ("<BioSimSpace.Protocol.FreeEnergyEquilibration: lam=%5.4f, lam_vals=%r, timestep=%s, "
                    "runtime=%s, temperature_start=%s, temperature_end=%s, pressure=%s, report_interval=%d, "
                    "restart_interval=%d, restart_interval=%d,restraint=%r, restart=%r, force_constant=%3.2f>"
                    ) % (self._lambda, self._lambda_vals, self._timestep, self._runtime,
                         self._temperature_start, self._temperature_end, self._pressure, self._report_interval,
                         self._restart_interval, self._restraint, self._restart, self._force_constant.value())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("<BioSimSpace.Protocol.FreeEnergyEquilibration(lam=%5.4f, lam_vals=%r, timestep=%s, "
                    "runtime=%s, temperature_start=%s, temperature_end=%s, pressure=%s, report_interval=%d, "
                    "restart_interval=%d, restart_interval=%d,restraint=%r, restart=%r, force_constant=%3.2f)>"
                    ) % (self._lambda, self._lambda_vals, self._timestep, self._runtime,
                         self._temperature_start, self._temperature_end, self._pressure, self._report_interval,
                         self._restart_interval, self._restraint, self._restart, self._force_constant.value())

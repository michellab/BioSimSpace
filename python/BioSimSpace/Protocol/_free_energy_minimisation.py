__all__ = ["FreeEnergyMinimisation"]

from ._free_energy_mixin import _FreeEnergyMixin
from ._minimisation import Minimisation as _Minimisation


class FreeEnergyMinimisation(_Minimisation, _FreeEnergyMixin):
    """A class for storing free energy minimisation protocols."""

    def __init__(self,
                 lam=0.0,
                 lam_vals=None,
                 min_lam=0.0,
                 max_lam=1.0,
                 num_lam=11,
                 steps=10000,
                 restraint=None,
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

           steps : int
               The maximum number of steps to perform.

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
        _Minimisation.__init__(self, steps=steps, restraint=restraint)

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
            return "<BioSimSpace.Protocol.FreeEnergyMinimisation: steps=%d, restraint=%r, lam=%5.4f, lam_vals=%r>" % (
                self._steps, self._restraint, self._lambda, self._lambda_vals)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return "<BioSimSpace.Protocol.FreeEnergyMinimisation: steps=%d, restraint=%r, lam=%5.4f, lam_vals=%r>" % (
                self._steps, self._restraint, self._lambda, self._lambda_vals)

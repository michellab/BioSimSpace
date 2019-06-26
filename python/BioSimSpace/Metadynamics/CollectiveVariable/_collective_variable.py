######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
Functionality for handling collective variables for metadynamics simulations.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["CollectiveVariable"]

class CollectiveVariable():
    """A base class for holding collective variables."""

    def __init__(self, pbc=True):
        """Constructor.

           Parameters
           ----------

           pbc : bool
              Whether to use periodic boundaries conditions.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is CollectiveVariable:
            raise Exception("<CollectiveVariable> must be subclassed.")

        # Whether this is a newly created object.
        self._is_new_object = True

        self.setPeriodicBoundaries(pbc)

    def setPeriodicBoundaries(self, pbc):
        """Set whether to use periodic_boundaries when calculating the
           collective variable.

           Parameters
           ----------

           pbc : bool
               Whether to use periodic boundaries conditions.
        """
        if type(pbc) is not bool:
            raise TypeError("'pbc' must be of type 'bool'")
        self._pbc = pbc

    def getPeriodicBoundaries(self):
        """Return whether to take account of periodic boundary conditions
           when computing the collective variable.

           Returns
           -------

           pbc : bool
               Whether to use periodic boundaries conditions.
        """
        return self._pbc

    @staticmethod
    def _setBound(bound, bound_type, name, force_constant=100.0, exponent=2.0, epsilon=1.0):
        """Internal function to validate input for a bound on a collective
           variable.

           Parameters
           ----------

           bound : dict, :class:`Type <BioSimSpace.Types._type.Type>`
               The value of the bound, along with parameters defining the
               bias potential.

           bound_type : :class:`Type <BioSimSpace.Types._types.Type>`
               The type of the bound, e.g. Length, Angle, etc.

           name : str
               The name of the bound ('lower_bound' or 'upper_bound').

           force_constant : float
               The force constant (k) for the bias potential.

           exponent : float
               The exponent (e) for the bias potential.

           epsilon : float
               The rescaling factor (s) for the bias potential.

           Returns
           -------

           bound : dict
               The validate bound dictionary.
        """

        if type(bound) is bound_type:
            bound = { "value"          : bound,
                      "force_constant" : force_constant,
                      "exponent"       : exponent,
                      "epsilon"        : epsilon }

        elif type(bound) is dict:
            pass

        else:
            raise TypeError("%r must be of type 'dict', or 'BioSimSpace.Types.%s'" % (name, bound_type))

        keys = bound.keys()

        if "value" not in keys:
            raise ValueError("Missing 'value' key for %r" % name)
        else:
            if type(bound["value"]) is not bound_type:
                raise ValueError("'value' must be of type 'BioSimSpace.Types.%s'" % bound_type.__qualname__)

        if "force_constant" not in keys:
            bound["force_constant"] = force_constant
        else:
            try:
                force_constant = float(bound["force_constant"])
                bound["force_constant"] = force_constant
            except:
                raise ValueError("'force_constant' must be of type 'float'")

        if "exponent" not in keys:
            bound["exponent"] = exponent
        else:
            try:
                exponent = float(bound["exponent"])
                bound["exponent"] = exponent
            except:
                raise ValueError("'exponent' must be of type 'float'")

        if "epsilon" not in keys:
            bound["epsilon"] = epsilon
        else:
            try:
                epsilon = float(bound["epsilon"])
                bound["epsilon"] = epsilon
            except:
                raise ValueError("'epsilon' must be of type 'float'")

        return bound

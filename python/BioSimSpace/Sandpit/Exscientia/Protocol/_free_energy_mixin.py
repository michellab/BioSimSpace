__all__ = ["_FreeEnergyMixin"]

from ._protocol import Protocol as _Protocol

class _FreeEnergyMixin(_Protocol):
    """A mixin for storing free energy protocols."""

    def __init__(self,
                 lam=0.0,
                 lam_vals=None,
                 min_lam=0.0,
                 max_lam=1.0,
                 num_lam=11,
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

        # Call the base class constructor.
        super().__init__()

        # Validate and set the lambda values.
        self.setLambdaValues(lam, lam_vals, min_lam, max_lam, num_lam)

        # Set the perturbation type. Default is "full", i.e. onestep protocol.
        self.setPerturbationType(perturbation_type)

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("<BioSimSpace.Protocol._FreeEnergyMixin: lam=%5.4f, lam_vals=%r>"
                   ) % (self._lambda, self._lambda_vals)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("<BioSimSpace.Protocol._FreeEnergyMixin: lam=%5.4f, lam_vals=%r>"
                    ) % (self._lambda, self._lambda_vals)

    def getPerturbationType(self):
        """Get the perturbation type.

           Returns
           -------

           perturbation_type : str
               The perturbation type.
        """
        return self._perturbation_type

    def setPerturbationType(self, perturbation_type):
        """Set the perturbation type.

           Parameters
           ----------

           perturbation_type : str
               The type of perturbation to perform. Options are:
                "full" : A full perturbation of all terms (default option).
                "discharge_soft" : Perturb all discharging soft atom charge terms (i.e. value->0.0).
                "vanish_soft" : Perturb all vanishing soft atom LJ terms (i.e. value->0.0).
                "flip" : Perturb all hard atom terms as well as bonds/angles.
                "grow_soft" : Perturb all growing soft atom LJ terms (i.e. 0.0->value).
                "charge_soft" : Perturb all charging soft atom LJ terms (i.e. 0.0->value).
        """
        if type(perturbation_type) is not str:
            raise TypeError("'perturbation_type' must be of type 'str'")

        # Convert to lower case and strip whitespace.
        perturbation_type = perturbation_type.lower().replace(" ", "")

        allowed_perturbation_types = ["full",
                                      "discharge_soft",
                                      "vanish_soft",
                                      "flip",
                                      "grow_soft",
                                      "charge_soft"]

        if perturbation_type not in allowed_perturbation_types:
            raise ValueError(f"'perturbation_type' must be one of: {allowed_perturbation_types}")

        self._perturbation_type = perturbation_type

    def getLambda(self):
        """Get the value of the perturbation parameter.

           Returns
           -------

           lam : float
               The value of the perturbation parameter.
        """
        return self._lambda

    def getLambdaValues(self):
        """Get the list of lambda values.

           Returns
           -------

           lam_vals : [float]
               The list of lambda values.
        """
        return self._lambda_vals

    def setLambdaValues(self, lam, lam_vals=None, min_lam=None, max_lam=None, num_lam=None):
        """Set the list of lambda values.

           Parameters
           ----------

           lam : float
               The perturbation parameter: [0.0, 1.0]

           lam_vals : [float]
               A list of lambda values.

           min_lam : float
               The minimum lambda value.

           max_lam : float
               The maximum lambda value.

           num_lam : int
               The number of lambda values.
        """

        # Convert int to float.
        if type(lam) is int:
            lam = float(lam)

        # Validate the lambda parameter.
        if type(lam) is not float:
            raise TypeError("'lam' must be of type 'float'.")

        self._lambda = lam

        # A list of lambda values takes precedence.
        if lam_vals is not None:
            # Convert tuple to list.
            if type(lam_vals) is tuple:
                lam_vals = list(lam_vals)

            # Make sure list (or tuple) has been passed.
            if type(lam_vals) is not list:
                raise TypeError("'lam_vals' must be of type 'list'.")

            # Make sure all lambda values are of type 'float'.
            if not all(isinstance(x, float) for x in lam_vals):
                raise TypeError("'lam_vals' must contain values of type 'float'.")

            # Make sure all values are in range [0.0, 1.0]
            for lam in lam_vals:
                if lam < 0:
                    raise ValueError("'lam_vals' must contain values in range [0.0, 1.0]")
                elif lam > 1:
                    raise ValueError("'lam_vals' must contain values in range [0.0, 1.0]")

            # Sort the values.
            lam_vals.sort()

            # Make sure there are no duplicates.
            if len(set(lam_vals)) < len(lam_vals):
                raise ValueError("'lam_vals' cannot contain duplicate values!")

            # Make sure the lambda value is in the list.
            if not lam in lam_vals:
                raise ValueError("'lam' is not a member of the 'lam_vals' list!")

            # Set the values.
            self._lambda_vals = lam_vals

        else:
            # Convert int to float.

            if type(min_lam) is int:
                min_lam = float(min_lam)

            if type(max_lam) is int:
                max_lam = float(max_lam)

            # Validate type.

            if type(min_lam) is not float:
                raise TypeError("'min_lam' must be of type 'float'.")

            if type(max_lam) is not float:
                raise TypeError("'max_lam' must be of type 'float'.")

            if type(num_lam) is not int:
                raise TypeError("'num_lam' must be of type 'int'.")

            # Validate values.

            if min_lam < 0 or min_lam > 1:
                raise ValueError("'min_lam' must be in range [0.0, 1.0]")

            if max_lam < 0 or max_lam > 1:
                raise ValueError("'max_lam' must be in range [0.0, 1.0]")

            if num_lam <=2:
                raise ValueError("'num_lam' must be greater than 2!")

            if min_lam >= max_lam:
                raise ValueError("'min_lam' must be less than 'max_lam'!")

            # Set values.

            self._lambda_vals = []
            step = (max_lam - min_lam) / (num_lam - 1)

            for x in range(0, num_lam):
                self._lambda_vals.append(round(min_lam + (x * step), 5))

            # Make sure the lambda value is in the list.
            if not self._lambda in self._lambda_vals:
                raise ValueError("'lam' (%5.4f) is not a member of the 'lam_vals' list: %s" \
                    % (self._lambda, self._lambda_vals))

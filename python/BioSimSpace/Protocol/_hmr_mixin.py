__all__ = ["_HmrMixin"]

from .. import Types as _Types

from ._protocol import Protocol as _Protocol


class _HmrMixin(_Protocol):
    """A mixin for storing HMR protocols."""

    def __init__(self,
                 hmr="auto",
                 hmr_factor="auto",
                 hmr_water="auto",
                 timestep=_Types.Time(2, "femtosecond")
                ):
        """Constructor.

           Parameters
           ----------

           hmr : "auto" or bool
               Whether HMR should be applied.
               If auto, this is based on the timestep and if a factor is chosen.

           hmr_factor : "auto" or float
               The factor used to repartition.
               "auto" indicates the recommended factor for the engine will be used.

           hmr_water : bool
               Whether the water molecules should also be repartitioned.
               "auto" indicates that recommended for the engine will be used.
           
           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration timestep.

        """

        # Call the base class constructor.
        super().__init__()

        # Validate the inputs
        self.setHmr(hmr, hmr_factor, hmr_water, timestep)

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("<BioSimSpace.Protocol._HmrMixin: hmr=%r, hmr_factor=%1.1f, hmr_water=%r>"
                   ) % (self._hmr, self._hmr_factor, self._hmr_water)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return ("<BioSimSpace.Protocol._HmrMixin: hmr=%r, hmr_factor=%1.1f, hmr_water=%r>"
                    ) % (self._hmr, self._hmr_factor, self._hmr_water)

    def getHmr(self):
        """Set the HMR flag.

           Returns
           -------

           hmr : bool
               Whether HMR should be applied.
        """
        return self._hmr

    def getHmrFactor(self):
        """Get the factor used for HMR.
           "auto" indicates the recommended factor for the engine will be used.

           Returns
           -------

           hmr_factor : "auto" or float
               The factor used to repartition.
        """
        return self._hmr_factor

    def getHmrWater(self):
        """Get whether HMR is applied to water molecules also.
          "auto" indicates that recommended for the engine will be used.

           Returns
           -------

           hmr_water : "auto" or bool
               Whether HMR should be applied to waters.
        """
        return self._hmr_water

    def setHmr(self, hmr, hmr_factor, hmr_water, timestep):
        """Set the HMR repartitioning variables.

           Parameters
           ----------

           hmr : "auto" or bool
               Whether HMR should be applied.

           hmr_factor : "auto" or float
               The factor used to repartition.
               "auto" indicates the recommended factor for the engine will be used.

           hmr_water : bool
               Whether the water molecules should also be repartitioned.
           
           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration timestep.

        """
        
        # check if the timestep is a BSS time.
        # This should be the same timestep as defined in the main protocol.
        if not isinstance(timestep, _Types.Time):
            raise TypeError("'timestep' must be of type 'BioSimSpace.Types.Time'")
          
        # check if the values are "auto" or bool for the rest
        if hmr.lower() != "auto":
            if not isinstance(hmr, bool):
                raise TypeError("'hmr' must be 'auto' or of type 'bool'.")
        elif hmr.lower() == "auto":
            hmr = "auto"

        if not isinstance(hmr_factor, str):
            if isinstance(hmr_factor, int):
                hmr_factor = float(hmr_factor)
            if not isinstance(hmr_factor, float):
                raise TypeError("'hmr_factor' must be 'auto' or of type 'int' or 'float'.")
            if hmr_factor > 5:
                raise ValueError("'hmr_factor' must not be greater than 5.")
        else:
            if hmr_factor.lower() == "auto":
                hmr_factor = "auto"
            else:
                raise TypeError("'hmr_factor' must be 'auto' or of type 'int' or 'float'.")

        if hmr_water.lower() != "auto":
            if not isinstance(hmr_water, bool):
                raise TypeError("'hmr_water' must be 'auto' or of type 'bool'.")
        elif hmr_water.lower() == "auto":
            hmr_water = "auto"

        # if hmr is auto, set to True when there is a timestep over 4fs
        if hmr == "auto":
            if timestep.femtoseconds().value() < 4:
                hmr = False
            elif timestep.femtoseconds().value() >= 4:
                hmr = True
        
        # however even if the timestep is less than 4fs,
        # if a hmr factor is chosen, also set hmr to true so this can be repartitioned.
        if hmr_factor != "auto":
            hmr = True

        # auto factors and water are later chosen when the engine is chosen for the protocol.

        # setting the values
        self._hmr = hmr # this will be either true or false
        self._hmr_factor = hmr_factor
        self._hmr_water = hmr_water

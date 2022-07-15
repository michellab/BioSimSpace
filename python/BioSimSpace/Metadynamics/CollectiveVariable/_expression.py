"""
Functionality for combining collective variables using a custom expression
"""

__author__ = "Adele Hardie"
__email__ = "adele.lip@gmail.com"

__all__ = ["Expression"]

from math import ceil as _ceil

from ._collective_variable import CollectiveVariable as _CollectiveVariable
from ...Types import Length as _Length
from ...Types import Angle as _Angle

# Store the collective variable base type.
#_colvar_type = _CollectiveVariable._collective_variable.CollectiveVariable

class Expression(_CollectiveVariable):
    """A class for combining collective variables using a custom expression."""

    def __init__(self, collective_variable, expression, variables=None, 
            arg=None, periodic=None, numerical_derivatives=False,
            hill_width=_Length(0.1, "nanometer"), lower_bound=None, upper_bound=None,
            grid=None):
        """Constructor.

            Parameters
            ----------
            collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                                [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
                The collective variable (or variables) used in the custom expression.
        
            expression : str
                The function to be evaluated. The CVs specified are referred 
                to either as x, y, and z, or as specified in 'var'.
        
            variables : [str]
                The names to give each of the CVs in the function. Can be 
                left as 'None' for up to 3 CVs, which will be called x, y,
                and z.

            arg : [str]
                specify which components will of variables will be used for 
                the calculation. E.g. 'distance' in 3D space has x, y, and 
                z components. To calculate the difference between the x
                components of two distances, use 'arg=['x.x', 'y.x']'. If
                'None', they will be set to be the same as 'variables'

            periodic : str
                if the output of the function is periodic, the periodicity
                should be specified here

            numerical_derivatives : bool
                whether to calculate the derivative of the function

            hill_width : :class:`Length <BioSimSpace.Types.Length>`
                The width of the Gaussian hill used to sample this variable.

            lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
                A lower bound on the value of the collective variable.

            upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
                An upper bound on the value of the collective variable.

            grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
                The grid on which the collective variable will be sampled.
                This can help speed up long metadynamics simulations where
                the number of Gaussian kernels can become prohibitive.        
        """

        # Call the base class constructor.
        super().__init__()
        
        # Set the collective variable
        self.setCollectiveVariable(collective_variable)

        # Set the types associated with his collective variable
        self.__setTypes()

        # Set the function
        self.setExpression(expression)

        # Set the variable names
        self.setVariableNames(variables)

        # set the arg
        self.setArg(arg)

        # set periodicity
        self.setPeriodic(periodic)

        # set derivatives
        self.setDerivatives(numerical_derivatives)

        # Set the "settable" parameters.
        self.setHillWidth(hill_width)

        # Set defaults for optional values.
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None

        # Set the optional parameters.
        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)
        if grid is not None:
            self.setGrid(grid)
    
    def __str__(self):
        """Return a human readable string representation of the object."""
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.Expression: "
        string += "collective variable(s)=%s " % ",".join(str(cv) for cv in self._collective_variable)
        string += "function=%s " % self._expression
        if self._arg is not None:
            string += "arguments=%s " % self._arg
        if self._periodic is not None:
            string += "periodic=%s " % self._periodic
        string += "derivatives=%s " % self._derivatives
        string += ", hill_width=%s" % self._hill_width
        if self._lower_bound is not None:
            string += ", lower_bound=%s" % self._lower_bound
        if self._upper_bound is not None:
            string += ", upper_bound=%s" % self._upper_bound
        if self._grid is not None:
            string += ", grid=%s" % self._grid
        string += ">"
        return string
    
    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()
    
    def setCollectiveVariable(self, collective_variable):
        """Set the collective variable (or variables).

        Parameters
        ----------
        collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                              [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
            The collective variable (or variables) for the simulation.
        """

        # A single collective variable.
        if isinstance(collective_variable, _CollectiveVariable):
            self._collective_variable = [collective_variable]

        else:
            if isinstance(collective_variable, (list, tuple)):
                if not all(isinstance(x, _CollectiveVariable) for x in collective_variable):
                    raise TypeError("'collective_variable' must all be of type "
                                    "'BioSimSpace.Metadynamics.CollectiveVariable'")
            else:
                raise TypeError("'collective_variable' must be of type "
                                "'BioSimSpace.Metadynamics.CollectiveVariable' "
                                "or a list of 'BioSimSpace.Metadynamics.CollectiveVariable' types.")

            self._collective_variable = collective_variable
    
    def getCollectiveVariable(self):
        """Return the collective variable(s)"""
        return self._collective_variable.copy()
    
    def setExpression(self, expression):
        """Set the custom expression to be evaluated"""
        if isinstance(expression, str):
            self._expression = expression
        else:
            raise TypeError("'expression' must be of type 'str'")

    def getExpression(self):
        """Return the expression"""
        return self._expression
    
    def setVariableNames(self, variables):
        """Set the collective variable names as they appear in the function expression
        
        Parameters
        ----------
        var : [str]
            Collective variable names
        """
        if isinstance(variables, str):
            variables = [variables]
        elif isinstance(variables, (list, tuple)):
            if not all(isinstance(x, str) for x in variables):
                raise TypeError("'variable' must be all of type str")

        # Check that the names match the collective variables:
        if variables is None:
            if len(self._collective_variable) > 3:
                raise ValueError("'var' must be defined when more than 3 CVs are used")
            else:
                variables = ["x", "y", "z"][:len(self.getCollectiveVariable())]       
        elif len(variables) != len(self._collective_variable):
            raise ValueError("Length of 'var' must match length of 'collective_variable'")

        # Check that the variables appear in the function:
        expression = self.getExpression()
        for var in variables:
            if var not in expression:
                raise ValueError("Variable with name '%s' does not appear in expression '%s'"%(var,expression))

        # Set variables
        self._variables = variables
        
    def getVariableNames(self):
        """Return list of variable names that map the CVs to the custom expression"""
        return self._variables
    
    def setArg(self, arg):
        """Set the arguments. If using different components, e.g. the
        difference between the x components of two distance vectors, use
        ['x.x', 'y.x']
        
        Parameters
        ----------
        arg : [str]
            List of arguments"""
        if arg is None:
            self._arg = self.getVariableNames()
        elif isinstance(arg, str):
            self._arg = [arg]
        elif isinstance(arg, (list,tuple)):
            if not all(isinstance(x, str) for x in arg):
                raise TypeError("'arg must all be of type 'str'")
            else:
                self._arg = arg
        else:
            raise TypeError("'arg' must be of type 'str' or a list of 'str'")

    def getArg(self):
        """Get the arguments"""
        return self._arg

    def setPeriodic(self, periodic):
        """Set the periodicity of the function if it is periodic. If
        'None', the function is not periodic
        
        Parameters
        ----------
        periodic : str
            The periodicity of the function"""
        if isinstance(periodic, (str)) or periodic is None:
            self._periodic = periodic
        else:
            raise TypeError("'periodic' must be 'None' or of type 'str'")

    def getPeriodic(self):
        """Get the periodicity of the function"""
        return self._periodic

    def setDerivatives(self, derivatives):
        """Set whether the derivative of the function has
        to be computed
        
        Parameters
        ----------
        derivatives : bool
            Whether to compute derivatives of the function"""
        if isinstance(derivatives, bool):
            self._derivatives = derivatives
        else:
            raise TypeError("'derivatives' must be of type 'bool'")
    
    def getDerivatives(self):
        """get whether the derivatives of the function are to be computed"""
        return self._derivatives


    def setHillWidth(self, hill_width):
        """Set the width of the Gaussian hills used to bias this collective
           variable.

           hill_width : :class:`Length <BioSimSpace.Types.Length>`
               The width of the Gaussian hill.
        """
        if not isinstance(hill_width, _Length):
            raise TypeError("'hill_width' must be of type 'BioSimSpace.Types.Length'")

        if hill_width.value() < 0:
            raise ValueError("'hill_width' must have a value of > 0")

        # Convert to the internal unit.
        self._hill_width = hill_width.nanometers()

    def getHillWidth(self):
        """Return the width of the Gaussian hill used to bias this collective
           variable.

           Returns
           -------

           hill_width : :class:`Length <BioSimSpace.Types.Length>`
               The width of the Gaussian hill.
        """
        return self._hill_width

    def __setTypes(self):
        """Set the types based on assigned CVs"""
        types = []
        for cv in self.getCollectiveVariable():
            if cv._types[0] not in types:
                types.append(cv._types[0])
        
        self._types = types
    
    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._lower_bound is not None:
            if type(self._lower_bound.getValue()) not in self._types:
                raise TypeError("'lower_bound' must be of type 'BioSimSpace.Types.Length'")
        if self._upper_bound is not None:
            if type(self._upper_bound.getValue()) not in self._types:
                raise TypeError("'upper_bound' must be of type 'BioSimSpace.Types.Length'")
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound.getValue() >= self._upper_bound.getValue():
                raise TypeError("'lower_bound' must less than 'upper_bound'")

        if self._grid is not None:
            if type(self._grid.getMinimum()) not in self._types:
                raise TypeError("'grid' minimum must be of type 'BioSimSpace.Types.Length'")
            if type(self._grid.getMaximum()) not in self._types:
                raise TypeError("Grid 'maximum' must be of type 'BioSimSpace.Types.Length'")
            if self._lower_bound is not None and self._grid.getMinimum() > self._lower_bound.getValue():
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if self._upper_bound is not None and self._grid.getMaximum() < self._upper_bound.getValue():
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")

            # If the number of bins isn't specified, estimate it out from the hill width.
            if self._grid.getBins() is None:
                grid_range = (self._grid.getMaximum() - self._grid.getMinimum()).value()
                num_bins = _ceil(5.0 * (grid_range / self._hill_width.value()))
                self._grid.setBins(num_bins)
        


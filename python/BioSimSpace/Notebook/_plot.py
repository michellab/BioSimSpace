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
Tools for plotting data.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["plot", "plotContour"]

from warnings import warn as _warn
from os import environ as _environ

from BioSimSpace import _is_interactive, _is_notebook
from BioSimSpace.Types._type import Type as _Type

# Check to see if DISPLAY is set.
if "DISPLAY" in _environ:
    _display = _environ.get("DISPLAY")
else:
    _display = None
del _environ

if _display is not None:
    _has_display = True
    try:
        import matplotlib.pyplot as _plt
        _has_matplotlib = True
    except ImportError:
        _has_matplotlib = False
else:
    if _is_notebook:
        try:
            import matplotlib.pyplot as _plt
            _has_matplotlib = True
        except ImportError:
            _has_matplotlib = False
    else:
        _has_matplotlib = False
        _has_display = False
        #_warn("The DISPLAY environment variable is unset. Plotting functionality disabled!")

del _display

if _has_matplotlib:
    # Define font sizes.
    _SMALL_SIZE = 14
    _MEDIUM_SIZE = 16
    _BIGGER_SIZE = 18

    # Set font sizes.
    _plt.rc('font', size=_SMALL_SIZE)          # controls default text sizes
    _plt.rc('axes', titlesize=_SMALL_SIZE)     # fontsize of the axes title
    _plt.rc('axes', labelsize=_MEDIUM_SIZE)    # fontsize of the x and y labels
    _plt.rc('xtick', labelsize=_SMALL_SIZE)    # fontsize of the tick labels
    _plt.rc('ytick', labelsize=_SMALL_SIZE)    # fontsize of the tick labels
    _plt.rc('legend', fontsize=_SMALL_SIZE)    # legend fontsize
    _plt.rc('figure', titlesize=_BIGGER_SIZE)  # fontsize of the figure title

def plot(x=None, y=None, xerr=None, yerr=None, xlabel=None, ylabel=None, logx=False, logy=False):
    """A simple function to create x/y plots with matplotlib.

       Parameters
       ----------

       x : list
           A list of x data values.

       y : list
           A list of y data values.

       xerr: list
           A list of error values for the x data.

       yerr: list
           A list of error values for the y data.

       xlabel : str
           The x axis label string.

       ylabel : str
           The y axis label string.

       logx : bool
           Whether the x axis is logarithmic.

       logy : bool
           Whether the y axis is logarithmic.
    """

    # Make sure were running interactively.
    if not _is_interactive:
        _warn("You can only use BioSimSpace.Notebook.plot when running interactively.")
        return None

    # Matplotlib failed to import.
    if not _has_matplotlib and _has_display:
        _warn("BioSimSpace.Notebook.plot is disabled as matplotlib failed "
            "to load. Please check your matplotlib installation.")
        return None

    # Convert tuple to a list.
    if type(x) is tuple:
        x = list(x)
    if type(y) is tuple:
        y = list(y)
    if type(xerr) is tuple:
        xerr = list(xerr)
    if type(yerr) is tuple:
        yerr = list(yerr)

    # Whether we need to convert the x and y data to floats.
    is_unit_x = False
    is_unit_y = False

    if x is None:
        if y is None:
            raise ValueError("'y' data must be defined!")

        # No x data, use array index as value.
        x = [x for x in range(0, len(y))]

    else:
        # No y data, we assume that the user wants to plot the x
        # data as a series.
        if y is None:
            y = x
            x = [x for x in range(0, len(y))]

    # The x argument must be a list of data records.
    if type(x) is not list:
        raise TypeError("'x' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(x[0])
        if not all(isinstance(xx, _type) for xx in x):
            raise TypeError("All 'x' data values must be of same type")

        # Convert int to float.
        if _type is int:
            x = [float(xx) for xx in x]
            _type = float

            try:
                xerr = [float(xx) for xx in xerr]
            except:
                pass

        # Make sure any associated error has the same unit.
        if xerr is not None:
            if not all(isinstance(xx, _type) for xx in xerr):
                raise TypeError("All 'xerr' values must be of same type as x data")

        # Does this type have units?
        if isinstance(x[0], _Type):
            is_unit_x = True

    # The y argument must be a list of data records.
    if type(y) is not list:
        raise TypeError("'y' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(y[0])
        if not all(isinstance(yy, _type) for yy in y):
            raise TypeError("All 'y' data values must be of same type")

        # Convert int to float.
        if _type is int:
            y = [float(yy) for yy in y]
            _type = float

            try:
                yerr = [float(yy) for yy in yerr]
            except:
                pass

        # Make sure any associated error has the same unit.
        if yerr is not None:
            if not all(isinstance(yy, _type) for yy in yerr):
                raise TypeError("All 'yerr' values must be of same type as y data")

        # Does this type have units?
        if isinstance(y[0], _Type):
            is_unit_y = True

    # Lists must contain the same number of records.
    # Truncate the longer list to the length of the shortest.
    if len(x) != len(y):
        _warn("Mismatch in list sizes: len(x) = %d, len(y) = %d"
            % (len(x), len(y)))

        len_x = len(x)
        len_y = len(y)

        if len_x < len_y:
            y = y[:len_x]
        else:
            x = x[:len_y]

        if xerr is not None:
            xerr = xerr[:len(x)]
        if yerr is not None:
            yerr = yerr[:len(y)]

    if xlabel is not None:
        if type(xlabel) is not str:
            raise TypeError("'xlabel' must be of type 'str'")
    else:
        if isinstance(x[0], _Type):
            xlabel = x[0].__class__.__qualname__ + " (" + x[0]._print_format[x[0].unit()] + ")"

    if ylabel is not None:
        if type(ylabel) is not str:
            raise TypeError("'ylabel' must be of type 'str'")
    else:
        if isinstance(y[0], _Type):
            ylabel = y[0].__class__.__qualname__ + " (" + y[0]._print_format[y[0].unit()] + ")"

    # Convert the x and y values to floats.
    if is_unit_x:
        x = [x.magnitude() for x in x]
        if xerr is not None:
            xerr = [x.magnitude() for x in xerr]
    if is_unit_y:
        y = [y.magnitude() for y in y]
        if yerr is not None:
            yerr = [y.magnitude() for y in yerr]

    # Set the figure size.
    _plt.figure(figsize=(8, 6))

    # Create the plot.
    if xerr is None and yerr is None:
        _plt.plot(x, y, "-bo")
    else:
        if xerr is None:
            _plt.errorbar(x, y, yerr=yerr, fmt="-bo")
        else:
            if yerr is None:
                _plt.errorbar(x, y, xerr=xerr, fmt="-bo")
            else:
                _plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="-bo")

    # Add axis labels.
    if xlabel is not None:
        _plt.xlabel(xlabel)
    if ylabel is not None:
        _plt.ylabel(ylabel)

    # Scale the axes.
    if logx:
        _plt.xscale("log")
    if logy:
        _plt.yscale("log")

    # Turn on grid.
    _plt.grid()

    return _plt.show()

def plotContour(x, y, z, xlabel=None, ylabel=None, zlabel=None):
    """A simple function to create two-dimensional contour plots with matplotlib.

       Parameters
       ----------

       x : list
           A list of x data values.

       y : list
           A list of y data values.

       z : list
           A list of z data values.

       xlabel : str
           The x axis label string.

       ylabel : str
           The y axis label string.

       zlabel : str
           The z axis label string.
    """

    import numpy as _np
    import scipy.interpolate as _interp

    from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable

    # Make sure were running interactively.
    if not _is_interactive:
        _warn("You can only use BioSimSpace.Notebook.plot when running interactively.")
        return None

    # Matplotlib failed to import.
    if not _has_matplotlib and _has_display:
        _warn("BioSimSpace.Notebook.plot is disabled as matplotlib failed "
            "to load. Please check your matplotlib installation.")
        return None

    # Convert tuple to a list.
    if type(x) is tuple:
        x = list(x)
    if type(y) is tuple:
        y = list(y)

    # Whether we need to convert the x, y, and z data to floats.
    is_unit_x = False
    is_unit_y = False
    is_unit_z = False

    # The x argument must be a list of data records.
    if type(x) is not list:
        raise TypeError("'x' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(x[0])
        if not all(isinstance(xx, _type) for xx in x):
            raise TypeError("All 'x' data values must be of same type")

        # Convert int to float.
        if _type is int:
            x = [float(xx) for xx in x]
            _type = float

        # Does this type have units?
        if isinstance(x[0], _Type):
            is_unit_x = True

    # The y argument must be a list of data records.
    if type(y) is not list:
        raise TypeError("'y' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(y[0])
        if not all(isinstance(yy, _type) for yy in y):
            raise TypeError("All 'y' data values must be of same type")

        # Convert int to float.
        if _type is int:
            y = [float(yy) for yy in y]
            _type = float

        # Does this type have units?
        if isinstance(y[0], _Type):
            is_unit_y = True

    if type(z) is not list:
        raise TypeError("'z' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(z[0])
        if not all(isinstance(zz, _type) for zz in z):
            raise TypeError("All 'z' data values must be of same type")

        # Convert int to float.
        if _type is int:
            z = [float(zz) for zz in z]
            _type = float

        # Does this type have units?
        if isinstance(z[0], _Type):
            is_unit_z = True

    # Lists must contain the same number of records.
    # Truncate the longer list to the length of the shortest.
    if len(x) != len(y) or \
       len(x) != len(z) or \
       len(y) != len(z):
        _warn("Mismatch in list sizes: len(x) = %d, len(y) = %d, len(z) = %d"
            % (len(x), len(y), len(z)))

        lens = [len(x), len(y), len(z)]
        min_len = min(lens)

        x = x[:min_len]
        y = y[:min_len]
        z = z[:min_len]

    if xlabel is not None:
        if type(xlabel) is not str:
            raise TypeError("'xlabel' must be of type 'str'")
    else:
        if isinstance(x[0], _Type):
            xlabel = x[0].__class__.__qualname__ + " (" + x[0]._print_format[x[0].unit()] + ")"

    if ylabel is not None:
        if type(ylabel) is not str:
            raise TypeError("'ylabel' must be of type 'str'")
    else:
        if isinstance(y[0], _Type):
            ylabel = y[0].__class__.__qualname__ + " (" + y[0]._print_format[y[0].unit()] + ")"

    if zlabel is not None:
        if type(zlabel) is not str:
            raise TypeError("'zlabel' must be of type 'str'")
    else:
        if isinstance(z[0], _Type):
            zlabel = z[0].__class__.__qualname__ + " (" + z[0]._print_format[z[0].unit()] + ")"

    # Convert the x and y values to floats.
    if is_unit_x:
        x = [x.magnitude() for x in x]
    if is_unit_y:
        y = [y.magnitude() for y in y]
    if is_unit_z:
        z = [z.magnitude() for z in z]

    # Convert to two-dimensional arrays. We don't assume the data is on a grid,
    # so we interpolate the z values.
    try:
        X, Y, = _np.meshgrid(_np.linspace(_np.min(x), _np.max(y), 1000),
                             _np.linspace(_np.min(y), _np.max(y), 1000))
        Z = _interp.griddata((x, y), z, (X, Y), method="linear")
    except:
        raise ValueError("Unable to interpolate x, y, and z data to a grid.")

    # Set the figure size.
    _plt.figure(figsize=(8, 8))

    # Create the contour plot.
    cp = _plt.contourf(X, Y, Z)

    # Add axis labels.
    if xlabel is not None:
        _plt.xlabel(xlabel)
    if ylabel is not None:
        _plt.ylabel(ylabel)

    # Get the current axes.
    ax = _plt.gca()

    # Make sure the axes are equal.
    ax.set_aspect("equal", adjustable="box")

    # Make sure the colour bar matches size of the axes.
    divider = _make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    # Add a colour bar and label it.
    cbar = _plt.colorbar(cp, cax=cax)
    if zlabel is not None:
        cbar.set_label(zlabel)

    return _plt.show()

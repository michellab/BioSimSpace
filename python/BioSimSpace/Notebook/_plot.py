######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from BioSimSpace import _is_interactive

from warnings import warn as _warn

try:
    import matplotlib.pyplot as _plt
    _has_matplotlib = True
except ImportError:
    _has_matplotlib = False

__all__ = ["plot"]

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

def plot(x=None, y=None, xlabel=None, ylabel=None, logx=False, logy=False):
    """A simple function to create x/y plots with matplotlib.

       Keyword arguments:

       x      -- A list of x data values.
       y      -- A list of y data values.
       xlabel -- The x axis label string.
       ylabel -- The y axis label string.
       logx   -- Whether the x axis is logarithmic.
       logy   -- Whether the y axis is logarithmic.
    """

    # Make sure were running interactively.
    if not _is_interactive():
        _warn("You can only use BioSimSpace.Notebook.plot when running interactively.")
        return None

    # Matplotlib failed to import.
    if not _has_matplotlib:
        _warn("BioSimSpace.Notebook.plot is disabled as matplotlib failed "
            "to load. Please check your matplotlib installation.")
        return None

    if x is None:
        if y is None:
            raise ValueError("'y' data must be defined!")
        else:
            if type(y) is not list:
                raise TypeError("'y' must be of type 'list'")
            # No x data, use array index as value.
            x = [x for x in range(0, len(y))]

    if type(x) is not list:
        raise TypeError("'x' must be of type 'list'")

    if len(x) != len(y):
        raise ValueError("Mismatch in list sizes: len(x) = %d, len(y) = %d"
            % (len(x), len(y)))

    if xlabel is not None:
        if type(xlabel) is not str:
            raise TypeError("'xlabel' must be of type 'str'")

    if ylabel is not None:
        if type(ylabel) is not str:
            raise TypeError("'ylabel' must be of type 'str'")

    # Set the figure size.
    _plt.figure(figsize=(8, 6))

    # Create the plot.
    _plt.plot(x, y, "-bo")

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

######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
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
Functionality for absolute free-energy simulations. Currently the Absolute
class is identical to the Relative class, but may diverge in future.
"""

__all__ = ["Absolute", "getData"]

from ._relative import Relative as _Relative


class Absolute(_Relative):
    """
    Class for configuring and running absolute free-energy perturbation simulations.
    This is currently an exact copy of the relative class but may diverge in future.
    """


def getData(name="data", file_link=False, work_dir=None):
    """Return a link to a zip file containing the data files required for
    post-simulation analysis.

    Parameters
    ----------

    name : str
        The name of the zip file.

    file_link : bool
        Whether to return a FileLink when working in Jupyter.

    work_dir : str
        The working directory for the simulation.

    Returns
    -------

    output : str, IPython.display.FileLink
        A path, or file link, to an archive of the process input.
    """
    return Absolute.getData(name=name, file_link=file_link, work_dir=work_dir)

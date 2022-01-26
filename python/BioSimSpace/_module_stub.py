######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
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
This file provides a class that can be used to stub any
module that BioSimSpace fails to import. The class will
raise a ModuleNotFound exception with a clear instruction
to the user if any attempt it made to use a module that 
has not been installed.
"""

__author__ = "Christopher Woods"
__email__ = "chryswoods@hey.com"

__all__ = ["module_stub", "try_import"]


class _ModuleStub:
    def __init__(self, name: str, install_command: str):
        self._name = name

        if install_command is None:
            self._install_command = f"conda install {name}"
        else:
            self._install_command = install_command

    def __getattr__(self, key):
        import BioSimSpace

        message = f"Cannot continue as the module {self._name} " \
                  "has not been installed. To continue, you " \
                  "should install the module using the command " \
                  f"{self._install_command}."

        if BioSimSpace._isVerbose():
            print(message)

        raise ModuleNotFoundError(message)


def module_stub(name: str, install_command: str = None):
    """Return a ModuleStub that will raise a ModuleNotFoundError
       if it is used in any way.

       Parameters
       ----------

       name: str
           The name of the module being stubbed

       install_command: str (optional)
           The command used to install the module. If
           this is not supplied, then it is assumed
           to be 'conda install {name}'

       Returns
       -------

       module: _ModuleStub
           The stubbed module
    """ 
    return _ModuleStub(name=name,
                       install_command=install_command)

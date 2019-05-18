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
.. currentmodule:: BioSimSpace.Node

Functions
=========

.. autosummary::
    :toctree: generated/

    help
    list
    run

Examples
========

Print the list of available nodes.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Node.list())

Get help on the "minimisation" node.

.. code-block:: python

   import BioSimSpace as BSS

   BSS.Node.help("minimisation")

Run the minimisation node.

.. code-block:: python

   import BioSimSpace as BSS

   # Generate a dictionary of input arguments.
   input = {"steps" : 1000, "files" : ["amber/ala/ala.top", "amber/ala/ala.crd"]}

   # Run the node and get the output as a dictionary.
   output = BSS.Node.run("minimisation", input)
"""

from glob import glob as _glob

import os as _os
import subprocess as _subprocess
import yaml as _yaml

import Sire.Base as _SireBase

# Set the node directory.
_node_dir = _os.path.dirname(__file__) + "/_nodes"

def list():
    """Return a list of the available nodes."""

    # Glob all Python scripts in the _nodes directory.
    nodes = _glob("%s/*.py" % _node_dir)

    # Strip the extension.
    nodes = [_os.path.basename(x).split(".py")[0] for x in nodes]

    return nodes

def help(name):
    """Print the help message for the named node.

       Parameters
       ----------

       name : str
           The name of the node.
    """

    if type(name) is not str:
        raise TypeError("'name' must be of type 'str'.")

    # Apped the node directory name.
    full_name = _node_dir + "/" + name

    # Make sure the node exists.
    if not _os.path.isfile(full_name):
        if not _os.path.isfile(full_name + ".py"):
            raise ValueError("Cannot find node: '%s'. " % name
                           + "Run 'Node.list()' to see available nodes!")
        else:
            full_name += ".py"

    # Create the command.
    command = "%s/python %s --help" % (_SireBase.getBinDir(), full_name)

    # Run the node as a subprocess.
    proc = _subprocess.run(command, shell=True, stdout=_subprocess.PIPE)

    # Print the standard output, decoded as UTF-8.
    print(proc.stdout.decode("utf-8"))

def run(name, args={}):
    """Run a node.

       Parameters
       ----------

       name : str
           The name of the node.

       args : dict
           A dictionary of arguments to be passed to the node.

       Returns
       -------

       output : dict
           A dictionary containing the output of the node.
    """

    # Validate the input.

    if type(name) is not str:
        raise TypeError("'name' must be of type 'str'.")

    if type(args) is not dict:
        raise TypeError("'args' must be of type 'dict'.")

    # Apped the node directory name.
    full_name = _node_dir + "/" + name

    # Make sure the node exists.
    if not _os.path.isfile(full_name):
        if not _os.path.isfile(full_name + ".py"):
            raise ValueError("Cannot find node: '%s'. " % name
                           + "Run 'Node.list()' to see available nodes!")
        else:
            full_name += ".py"

    # Write a YAML configuration file for the BioSimSpace node.
    if len(args) > 0:
        with open("input.yaml", "w") as file:
            _yaml.dump(args, file, default_flow_style=False)

        # Create the command.
        command = "%s/python %s --config input.yaml" % (_SireBase.getBinDir(), full_name)

    # No arguments.
    else:
        command = "%s/python %s" % (_SireBase.getBinDir(), full_name)

    # Run the node as a subprocess.
    proc = _subprocess.run(command, shell=True, stderr=_subprocess.PIPE)

    if proc.returncode == 0:
        # Read the output YAML file into a dictionary.
        with open("output.yaml", "r") as file:
            output = _yaml.safe_load(file)

        # Delete the redundant YAML files.
        _os.remove("input.yaml")
        _os.remove("output.yaml")

        return output

    else:
        # Print the standard error, decoded as UTF-8.
        print(proc.stderr.decode("utf-8"))

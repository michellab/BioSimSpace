
BioSimSpace.Process
===================

This package provides functionality for running different molecular
simulation processes.

The package provides a base class, ``Process``\ , that defines common properties
and methods for all processes. Derived classes, such as ``Namd``\ , define
functionality for running process with a particular software package. At
present we provide support for `AMBER <http://ambermd.org>`_\ ,
`GROMACS <http://www.gromacs.org/>`_\ , `NAMD <http://www.ks.uiuc.edu/Research/namd>`_\ ,
and `SOMD <https://siremol.org/tutorials/somd>`_. (GROMACS support is
currently limited.)

Object instantiation
--------------------

All process classes must take at least two arguments to their constructor:


* 
  ``system``\ : A Sire molecular system.

* 
  ``protocol``\ : A `\ ``BioSimSpace.Protocol`` <../Protocol>`_ object defining the
  protocol for the simulation process, e.g. an equilibration protocol.

For example, to initialise an object to run a default minimistion protocol
using AMBER:

.. code-block:: python

   import BioSimSpace as BSS

   # Create a molecular system.
   system = BSS.IO.readMolecules(["ala.crd", "ala.top"])

   # Create a default minimisation protocol.
   protocol = BSS.Protocol.Minimisation()

   # Initialise the AMBER process.
   process = BSS.Process.Amber(system, protocol)

*Working directory*
^^^^^^^^^^^^^^^^^^^^^^^

By default, each process is run in a temporary workspace. To specify
the working directory the user can pass an appropriate keyword argument:

.. code-block:: python

   # Initialise the AMBER process using a custom working directory.
   process = BSS.Process.Amber(system, work_dir="/my/custom/path")

The directory will be created if it doesn't already exist (assuming write
privileges on the path).

To change the location of temporary directory, e.g. if ``/tmp`` doesn't have
enough space, simply set the ``TMPDIR`` environment variable before using
BioSimSpace.

*Executable*
^^^^^^^^^^^^^^^^

BioSimSpace will search your ``PATH`` to find an appropriate executable to run
the process. An ``IOError`` will be raised if the executable is missing.
Alternatively, the location of the executable can be specified when creating
the process object, e.g.

.. code-block:: python

   # Initialise the AMBER process and specify the executable path.
   process = BSS.Process.Amber(system, exe="/home/amber/bin/sander")

*Configuration parameters*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once initialised, the process object will have set up all of the appropriate
input and configuration files needed to run the desired simulation protocol
using AMBER.

To get a list of the auto-generated input files:

.. code-block:: python

   # Return a list of names for the input files.
   files = process.inputFiles()

To get a list of the configuration file options:

.. code-block:: python

   # Return the contents of the configuration file as a list of strings.
   config = process.getConfig()

BioSimSpace uses a set of well chosen configuration parameters as defaults for
each simulation protocol. However, we provide lots of flexibility for overriding
these defaults. For example:

.. code-block:: python

   # Add a single additional configuration string.
   param = "some-parameter = some-value"
   process.addToConfig(param)

   # Add a list of additional configuration parameters.
   params = ["some-parameter = some-value", "some-other-parameter = some-other-value"]
   process.addToConfig(params)

   # Add some additional parameters from a file.
   process.addToConfig("params.txt")

   # Overwrite the entire configuration using a new set of parameters.
   process.setConfig(params)         # Using a list of parameter strings.
   process.setConfig("params.txt")   # Using a parameter file.

   # Write the current configuration parameters to a file.
   process.writeConfig("params.txt")

*Command-line arguments*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If necessary, BioSimSpace also configures all of the command-line arguments
needed to run the process. To get the arguments:

.. code-block:: python

   # Get the arguments as an OrderedDict, i.e. ([arg, value], ...)
   arg_dict = process.getArgs()

   # Get the arguments as a list of strings.
   arg_strings = process.getArgStringList()

   # Get the arguments as a single single string.
   arg_string = process.getArgString()

As with configuration parameters, we provide a flexible means of configuring
the command-line arguments.

.. code-block:: python

   # Add a single additional command-line argument. (This overwrites any existing argument with the same name.)
   process.setArg("-inf", "mdinfo.txt")    # Regular argument, i.e. arg / value.
   process.setArg("-O", True)              # Boolean flag.

   # Add a dictionary of arguments. (A regular 'dict' is allowed, although argument ordering is lost.)
   args = OrderedDict([('-inf', 'mdinfo.txt'), ('-O', True)])
   process.addArgs(args)

   # Insert an additional argument at a specfic position.
   process.insertArg("-inf", "mdinfo.txt", 3)

   # Delete an argument.
   process.deleteArg("-inf")

   # Disable/enable a boolean argument.
   process.setArg("-O", False)
   process.setArg("-O", True)

   # Overwrite all command-line arguments. (A regular 'dict' is allowed, although argument ordering is lost.)
   args = OrderedDict([('-inf', 'mdinfo.txt'), ('-O', True)])
   process.setArgs(args)

   # Clear the command-line arguments.
   process.clearArgs()

Running a process
-----------------

Once you are happy with the way a process is configured it can be started using:

.. code-block:: python

   # Start the process in the background and return to the main thread.
   process.start()

If you are using BioSimSpace from a regular python script then the process will
block the main thread if you try to get any data from it while it is running, e.g.:

.. code-block:: python

   # This will wait for the process to finish running before returning the final system.
   system = process.getSystem()

If you are using BioSimSpace interactively, e.g. using a `Jupyter <http://jupyter.org>`_
notebook, you can carry on with your work and query the running process in
real time. Each BioSimSpace.Process object collects output from ``stdout`` and ``stderr``
and monitors log files for updates to thermodynamic measures, such as energy
and pressure. Some examples of how to interactively query a running process
are given below.

.. code-block:: python

   # Check whether the process is still running.
   process.isRunning()

   # Get the estimated number of minutes until completion. This is supported for all
   # programs that provide regular timing statistics.
   process.eta()

   # Get the current runtime of the process (in minutes).
   runtime = process.runTime()

   # Print the last 10 lines from stdout.
   process.stdout()

   # Print the last 20 lines from stderr.
   process.stderr(20)

   # Get the whole of stdout and stderr a list of strings.
   stdout = process.getStdout()
   stderr = process.getStderr()

   # Get the current number of steps, the run time (in nanoseconds) and energy (in kcal/mol).
   # Many other record types are supported. The options available depend on the nature of the
   # program and simulation protocol.
   step = process.getStep()
   time = process.getTime()
   energy = process.getTotalEnergy()

   # Since the output is recorded periodically, we can also get a time series of records.
   step = process.getStep(time_series=True)
   time = process.getTime(time_series=True)
   energy = process.getTotalEnergy(time_series=True)

It is also possible to get the latest molecular system from the running process:

.. code-block:: python

   # Return the latest molecular configuration as a Sire.System.
   system = process.getSystem()

This could then be saved to file:

.. code-block:: python

   BSS.IO.saveMolecules("configuration", system, system.fileFormat())

To safely kill a running process:

.. code-block:: python

   process.kill()

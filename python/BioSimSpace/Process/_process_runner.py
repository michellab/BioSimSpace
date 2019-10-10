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
Functionality for running multiple processes.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["ProcessRunner"]

import os as _os
import tempfile as _tempfile

from BioSimSpace._SireWrappers import System as _System

from ._process import Process as _Process

class ProcessRunner():
    """A class for managing and running multiple simulation processes, e.g.
       a free energy simulation at multiple lambda values."""

    def __init__(self, processes, name="runner", work_dir=None, nest_dirs=True):
        """Constructor.

           Parameters
           ----------

           processes : [:class:`Process <BioSimSpace.Process>`]
               A list of process objects.

           name : str
               The name of the of processes.

           work_dir : str
               The working directory for the processes.

           nest_dirs : bool
               Whether to nest the working directory of the processes inside
               the process runner top-level directory.
        """

        # Check that the list of processes is valid.
        if not all(isinstance(process, _Process) for process in processes):
            raise TypeError("'processes' must be a list of 'BioSimSpace.Process' types.")

        # Make sure all of the processes aren't running.
        if not all(process.isRunning() == False for process in processes):
            raise ValueError("'processes' must not contain any running 'BioSimSpace.Process' objects!")

        # Check that the working directory is valid.
        if work_dir is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        # Check that the nest_dirs flag is valid.
        if type(nest_dirs) is not bool:
            raise TypeError("'nest_dirs' must be of type 'bool'")
        self._nest_dirs = nest_dirs

        # Set the list of processes.
        self._processes = processes

        # Set the name
        if name is None:
            self._name = None
        else:
            self.setName(name)

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir, exist_ok=True)

        # Nest all of the process working directories inside the runner directory.
        if nest_dirs:
            self._processes = self._nest_directories(self._processes)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Process.%s: nProcesses=%d, nRunning=%d, nQueued=%d, nError=%d, name='%s', work_dir='%s'>" \
            % (self.__class__.__name__, self.nProcesses(), self.nRunning(), self.nQueued(),
               self.nError(), self._name, self._work_dir)

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Process.%s: nProcesses=%d, nRunning=%d, nQueued=%d, nError=%d, name='%s', work_dir='%s'>" \
            % (self.__class__.__name__, self.nProcesses(), self.nRunning(), self.nQueued(),
               self.nError(), self._name, self._work_dir)

    def processes(self):
        """Return the list of processes.

           Returns
           -------

           processes : [:class:`Process<BioSimSpace.Process>`]
               The list of processes.
        """
        return self._processes

    def workDir(self):
        """Return the working directory.

           Returns
           -------

           work_dir : str
               The working directory.
        """
        return self._work_dir

    def getName(self):
        """Return the process runner name.

           Returns
           -------

           name : str
               The name of the process.

        """
        return self._name

    def setName(self, name):
        """Set the process runner name.

           Parameters
           ----------

           name : str
               The process runner name.
        """
        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'")
        else:
            self._name = name

    def addProcess(self, process):
        """Add a process to the runner.

           Parameters
           ----------

           process : :class:`Process <BioSimSpace.Process>`, \
                     [:class:`Process <BioSimSpace.Process>`]
               The process/processes to add.
        """

        # Convert to a list.
        if type(process) is not list:
            processes = [ process ]

        # Check that the list of processes is valid.
        if not all(isinstance(process, _Process) for process in processes):
            raise TypeError("'processes' must be a list of 'BioSimSpace.Process' types.")

        # Make sure all of the processes aren't running.
        if not all(process.isRunning() == False for process in processes):
            raise ValueError("'processes' must not contain any running 'BioSimSpace.Process' objects!")

        # Nest the directories inside the process runner's working directory.
        if self._nest_dirs:
            # Extend the list of procesess.
            self._processes.extend(self._nest_directories(processes, len(self._processes)))

    def removeProcess(self, index):
        """Remove a process from the runner.

           Parameters
           ----------

           index : int
               The index of the process.
        """

        try:
            # Pop the chosen process from the list.
            process = self._processes.pop(index)

            # Kill the process.
            process.kill()

        except IndexError:
            raise("'index' is out of range: [0-%d]" % len(self._processes))

    def nProcesses(self):
        """Return the number of processes.

           Returns
           -------

           n_processes : int
               The number of processes managed by the runner.
        """
        return len(self._processes)

    def nRunning(self):
        """Return the number of running processes.

           Returns
           -------

           n_running : int
               The number of processes that are running.
        """

        n = 0

        for p in self._processes:
            if p.isRunning():
                n += 1

        return n

    def nQueued(self):
        """Return the number of queued processes.

           Returns
           -------

           n_queued : int
               The number of processes that are queued.
        """

        n = 0

        for p in self._processes:
            if p.isQueued():
                n += 1

        return n

    def nError(self):
        """Return the number of errored processes.

           Returns
           -------

           n_error : int
               The number of processes that are in an error state.
        """

        n = 0

        for p in self._processes:
            if p.isError():
                n += 1

        return n

    def running(self):
        """Return the indices of the running processes.

           Returns
           -------

           idx_running : [ int ]
               A list containing the indices of the running processes.
        """

        indices = []

        for idx, p in enumerate(self._processes):
            if p.isRunning():
                indices.append(idx)

        return indices

    def queued(self):
        """Return the indices of the queued processes.

           Returns
           -------

           idx_queued : [int]
               A list containing the indices of the queued processes.
        """

        indices = []

        for idx, p in enumerate(self._processes):
            if p.isQueued():
                indices.append(idx)

        return indices

    def errored(self):
        """Return the indices of the errored processes.

           Returns
           -------

           idx_errored : [int]
               A list containing the indices of the errored processes.
        """

        indices = []

        for idx, p in enumerate(self._processes):
            if p.isError():
                indices.append(idx)

        return indices

    def isRunning(self):
        """Return whether each process is running.

           Returns
           -------

           is_running : [ bool ]
               A list indicating whether each process is running.
        """

        bool_list = []

        for p in self._processes:
            if p.isRunning():
                bool_list.append(True)
            else:
                bool_list.append(False)

        return bool_list

    def isQueued(self):
        """Return whether each process is queued.

           Returns
           -------

           is_queued : [ bool ]
               A list indicating whether each process is queued.
        """

        bool_list = []

        for p in self._processes:
            if p.isQueued():
                bool_list.append(True)
            else:
                bool_list.append(False)

        return bool_list

    def isError(self):
        """Return whether each process is in an error state.

           Returns
           -------

           is_error : [bool]
               A list indicating whether each process is in an error state.
        """

        bool_list = []

        for p in self._processes:
            if p.isError():
                bool_list.append(True)
            else:
                bool_list.append(False)

        return bool_list

    def start(self, index):
        """Start a specific process. The same can be achieved using:
               runner.processes()[index].start()

           Parameters
           ----------

           index : int
               The index of the process.
        """

        try:
            self._processes[index].start()

        except IndexError:
            raise("'index' is out of range: [0-%d]" % len(self._processes))

    def startAll(self):
        """Start all of the processes."""

        for p in self._processes:
            # Initialise the error state.
            is_error = True

            # Zero the tally of failed processes.
            num_failed = 0

            # Retry failed processes up to a maximum of 5 times.
            while is_error:
                # Start the process and wait for it to finish.
                p.start()
                p.wait()

                # Check the error state.
                is_error = p.isError()

                # Increment the number of failures.
                if is_error:
                    num_failed += 1

                    # Maximum retries reached, move to the next process.
                    if num_failed == 5:
                        break

    def kill(self, index):
        """Kill a specific process. The same can be achieved using:
               runner.processes()[index].kill()

           Parameters
           ----------

           index : int
               The index of the process.
        """

        try:
            self._processes[index].kill()

        except IndexError:
            raise("'index' is out of range: [0-%d]" % len(self._processes))

    def killAll(self):
        """Kill all of the processes."""

        for p in self._processes:
            p.kill()

    def restartFailed(self):
        """Restart any jobs that are in an error state."""

        for p in self._processes:
            if p.isError():
                p.start()

    def runTime(self):
        """Return the run time for each process.

           Returns
           -------

           run_time : [ BioSimSpace.Types.Time ]
               A list containing the run time of each process.
        """

        run_time = []

        for p in self._processes:
            run_time.append(p.runTime())

        return run_time

    def _nest_directories(self, processes):
        """Helper function to nest processes inside the runner's working
           directory.

           Parameters
           ----------

           processes : [:class:`Process <BioSimSpace.Process>`]
               A list of process objects.

           Returns
           -------

           new_processes : [:class:`Process <BioSimSpace.Process>`]
               A list of procesess with updated working directories.
        """

        # Create the list of new processes.
        new_processes = []

        # Loop over each process.
        for process in processes:
            # Create the new working directory name.
            new_dir = "%s/%s" % (self._work_dir, process._work_dir)

            # Create a new process object using the nested directory.
            if process._package_name == "SOMD":
                new_processes.append(type(process)(_System(process._system), process._protocol,
                    process._exe, process._name, process._platform, new_dir, process._seed, process._property_map))
            else:
                new_processes.append(type(process)(_System(process._system), process._protocol,
                    process._exe, process._name, new_dir, process._seed, process._property_map))

        return new_processes

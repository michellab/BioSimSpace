
BioSimSpace: Unit Tests
=======================

This directory contains all of the unit tests for BioSimSpace.

The following has been adapted from a testing workshop available
`here <http://chryswoods.com/python_and_data/testing>`_.

Directory layout
----------------

The main directory is ``test``. Inside this are subdirectories for each of the
BioSimSpace sub-modules, which include files containing all of the unit tests.

The test suite is intended to be run using `pytest <https://docs.pytest.org/en/latest/contents.html>`_.
When run, ``pytest`` searches for tests in all directories and files below the current
directory, collects the tests together, then runs them. Pytest uses name matching
to locate the tests. Valid names start or end with *test*\ , e.g.:

.. code-block:: python

   # Files:
   test_file.py       file_test.py

   # Functions:
   def test_func():   def func_test():

We use the convention of ``test_*`` when naming files and functions.

Running tests
-------------

To run the full test suite, simply type:

.. code-block:: bash

   pytest

This assumes that you have installed `Sire <https://github.com/michellab/Sire>`_
and have added its ``bin`` directory, ``$HOME/sire.app/bin``\ , to your path.

We also assume that you have installed `BioSimSpace <https://github.com/michellab/BioSimSpace>`_
into an existing `Sire <https://github.com/michellab/Sire>`_ package, making
it is accessible to ``pytest``. If not, you can prefix the ``pytest`` commands with
an appropriate ``PYTHONPATH`` environment variable. For example:

.. code-block:: bash

   PYTHONPATH=/home/Code/BioSimSpace/python pytest

To run tests for a specific sub-module, e.g.:

.. code-block:: bash

   pytest test/Process

To only run the unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest test/Process/test_namd.py

To run a specific unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest test/Process/test_namd.py::test_minimise

To get more detailed information about each test, run pytests using the
*verbose* flag, e.g.:

.. code-block::

   pytest -v

More details regarding how to invoke ``pytest`` can be found `here <https://docs.pytest.org/en/latest/usage.html>`_.

Writing tests
-------------

Basics
^^^^^^

Try to keep individual unit tests short and clear. Aim to test one thing, and
test it well. Where possible, try to minimise the use of ``assert`` statements
within a unit test. Since the test will return on the first failed assertion,
additional contextual information may be lost.

Floating point comparisons
^^^^^^^^^^^^^^^^^^^^^^^^^^

Make use of the `\ ``approx`` <https://docs.pytest.org/en/latest/builtin.html#comparing-floating-point-numbers>`_
function from the ``pytest`` package for performing floating point comparisons, e.g:

.. code-block:: python

   from pytest import approx

   assert 0.1 + 0.2 == approx(0.3)

By default, the ``approx`` function compares the result using a relative tolerance
of 1e-6. This can be changed by passing a keyword argument to the function, e.g:

.. code-block:: python

   assert 2 + 3 == approx(7, rel=2)

Skipping tests
^^^^^^^^^^^^^^

If you are using `\ *test-driven development* <https://en.wikipedia.org/wiki/Test-driven_development>`_
it might be desirable to write your tests before implementing the functionality,
i.e. you are asserting what the *output* of a function should be, not how it should
be *implemented*. In this case, you can make use of the ``pytest`` *skip* decorator
to flag that a unit test should be skipped, e.g.:

.. code-block:: python

   @pytest.mark.skip(reason="Not yet implemented.")
   def test_new_feature():
       # A unit test for an, as yet, unimplemented feature.
       ...

Parametrizing tests
^^^^^^^^^^^^^^^^^^^

Often it is desirable to run a test for a range of different input parameters.
This can be achieved using the ``parametrize`` decorator, e.g.:

.. code-block:: python

   import pytest
   from operator import mul

   @pytest.mark.parametrize("x", [1, 2])
   @pytest.mark.parametrize("y", [3, 4])
   def test_mul(x, y):
       """ Test the mul function. """
       assert mul(x, y) == mul(y, x)

Here the function test_mul is parametrized with two parameters, ``x`` and ``y``.
By marking the test in this manner it will be executed using all possible
parameter pairs ``(x, y)``\ , i.e. ``(1, 3), (1, 4), (2, 3), (2, 4)``.

Alternatively:

.. code-block:: python

   import pytest
   from operator import sub
   @pytest.mark.parametrize("x, y, expected",
                           [(1, 2, -1),
                            (7, 3,  4),
                            (21, 58, -37)])
   def test_sub(x, y, expected):
       """ Test the sub function. """
       assert sub(x, y) == -sub(y, x) == expected

Here we are passing a list containing different parameter sets, with the names
of the parameters matched against the arguments of the test function.

Testing exceptions
^^^^^^^^^^^^^^^^^^

Pytest provides a way of testing your code for known exceptions. For example,
suppose we had a function that raises an ``IndexError``\ :

.. code-block:: python

   def indexError():
       """ A function that raises an IndexError. """
       a = []
       a[3]

We could then write a test to validate that the error is thrown as expected:

.. code-block:: python

   def test_indexError():
       with pytest.raises(IndexError):
           indexError()

Custom attributes
^^^^^^^^^^^^^^^^^

It's possible to mark test functions with any attribute you like. For example:

.. code-block:: python

   @pytest.mark.slow
   def test_slow_function():
       """ A unit test that takes a really long time. """
       ...

Here we have marked the test function with the attribute ``slow`` in order to
indicate that it takes a while to run. From the command line it is possible
to run or skip tests with a particular mark.

.. code-block:: python

   pytest mypkg -m "slow"        # only run the slow tests
   pytest mypkg -m "not slow"    # skip the slow tests

The custom attribute can just be a label, as in this case, or could be your
own function decorator.

Committing
----------

If you create new tests, please make sure that they pass locally before
pushing your commits to the remote.

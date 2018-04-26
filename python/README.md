# Python

BioSimSpace uses the Python programming language. Our aim is to provide a simple
and robust API where unnecessary implementation details are hidden from the user.

"_Hold on a second, this code isn't very Pythonic!_"

Indeed it is not, but this is a design choice. BioSimSpace is intended primarily
to be used by novices, who may be unfamiliar with Python, or programming in
general. We want to make it as easy as possible for these users to get up and
running with molecular simulation. BioSimSpace also needs to be robust and
portable, hence we need to use encapsulation to shield the user from unintended
consequences.

With this in mind, we use the following coding conventions:

## Naming

We follow a C++ style naming convention.

* Packages: CamelCase
* Classes: CamelCase
* Methods: camelCase
* Functions: camelCase
* Variables: snake_case

For example, to instantiate a minimisation protocol from the `Protocol` package:

```python
import BioSimSpace as BSS
protocol = BSS.Protocol.Minimisation()
```

## Modules

BioSimSpace is a collection of packages, e.g. `BioSimSpace.Gateway` and
`BioSimSpace.Protocol`. Within each package is a set of modules that
implement the required functionality. Rather than directly exposing all of
the modules we choose to hide implementation details from the user. Instead
we use the package `__init__.py` to selectively import the required
classes and functions.

* Module files containing implementation details are prefixed with an underscore,
i.e. `_process.py`

* Each module file contains an `__all__` variable that lists the specific items
that should be imported.

* The package `__init__.py` can be used to safely expose the required
functionality to the user with:

```python
from module import *
```

This results in a clean API and documentation, with all extraneous information,
e.g. external modules, hidden from the user. This is important when working
interactively, since [IPython](https://ipython.org) and (Jupyter)[https://jupyter.org]
do not respect the `__all__` variable when auto-completing, meaning that the
user will see a full list of the available names when hitting tab. When
following the conventions above, the user will only be able to access the
exposed names. This greatly improves the clarity of the package, allowing
a new user to quickly determine the available functionality. Any user wishing
expose further implementation detail can of course type an underscore to
show the hidden names when searching.

(An alternative method is to prefix all external imports local to a module
with an underscore, e.g. `import sys as _sys`. However, this quickly becomes
unwieldy and also loses the flexibility of the above approach.)

## Encapsulation

BioSimSpace aims to provide a means of writing robust and portable workflow
components (nodes). To this end, we choose to use an object oriented approach
where data is encapsulated, with getters used to retrieve data from an object.

To avoid unintended consequences, getters that return mutable data types, e.g.
lists and dictionaries, should return a copy of the data. This prevents the
user unintentionally modifying the private data contained in the object. Setters
should be used to explicitly modify member data.

For example:

```python
# A class that holds a list of numbers.

class MyClass():
    # A private class member variable containing a list of numbers.
    _list = [1, 2, 3, 4, 5]

    def getList(self):
        return self._list

# Create an instance of the class.
c = MyClass()
n = c.getList()
print(n)
[1, 2, 3, 4, 5]

# Update n.
n.append(6)

# The private member data has been modified!
print(c.getList()
[1, 2, 3, 4, 5, 6]
```

Instead use:

```python
class MyClass():
    # A private class member variable containing a list of numbers.
    _list = [1, 2, 3, 4, 5]

    def getList(self):
        return self._list.copy()

# Create an instance of the class.
c = MyClass()
n = c.getList()
print(n)
[1, 2, 3, 4, 5]

# Update n.
n.append(6)

# The private member data is untouched.
print(c.getList()
[1, 2, 3, 4, 5]
```

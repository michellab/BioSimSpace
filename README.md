# BioSimSpace

Code and resources for the EPSRC [BioSimSpace](https://biosimspace.org) project.

To install into a [Sire](https://github.com/michellab/Sire) package, use:

```bash
cd python
python setup.py install
```

Note that this assumes that Sire's `bin` directory, `$HOME/sire.app/bin`,
has been added to your path. Alternatively, run:

```bash
$HOME/sire.app/bin/python setup.py install
```
At present we recommend using the `devel` branch of [Sire](https://github.com/michellab/Sire)
to ensure that you have access to the latest features and bug-fixes.

## Documentation

Each package has its own README page:

- [BioSimspace.Gateway](python/BioSimSpace/Gateway)
- [BioSimspace.IO](python/BioSimSpace/IO)
- [BioSimspace.MD](python/BioSimSpace/MD)
- [BioSimspace.Notebook](python/BioSimSpace/Notebook)
- [BioSimspace.Process](python/BioSimSpace/Process)
- [BioSimspace.Protocol](python/BioSimSpace/Protocol)
- [BioSimspace.Trajectory](python/BioSimSpace/Trajectory)

Full API documentation and examples can be found at [biosimspace.org](https://biosimspace.org).

## Dependencies

BioSimSpace makes use of several external python packages. These should be
automatically installed into your Sire app by the [setup.py](python/setup.py)
script.

## Developing

BioSimSpace is written in Python but uses an object-oriented style with C++
naming conventions. More details can be found [here](python).

Please create a feature branch for development work. Branches related to a
particular feature should be prefixed with `feature`, i.e. for the
BioSimSpace.Process module:

```bash
feature-process
```

When you are happy with your feature, create a pull request so that the feature
can be merged into devel. Once the merge is successful, please delete the
redundant feature branch.

## Tests

A test suite can be found [here](https://github.com/michellab/BioSimSpaceUnitTests).

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
At present you will need to use the `devel` branch of [Sire](https://github.com/michellab/Sire)
to make sure that you have access to the latest features and bug-fixes.

## Documentation

Each sub-package has its own README page:

- [BioSimspace.Driver](python/BioSimSpace/Driver)
- [BioSimspace.Gateway](python/BioSimSpace/Gateway)
- [BioSimspace.Notebook](python/BioSimSpace/Notebook)
- [BioSimspace.Process](python/BioSimSpace/Process)
- [BioSimspace.Protocol](python/BioSimSpace/Protocol)
- [BioSimspace.Trajectory](python/BioSimSpace/Trajectory)

## Dependencies

BioSimSpace makes use of several external python packages. When using BioSimSpace
these should be automatically installed into your Sire package using
`try_import` from `Sire.Base`. In some instances it may appear that the
import has failed even though the package was installed successfully. If this
occurs, simply re-run your script (the module should now be imported without
error.)

## Developing

Please create a feature branch for development work. We use the following
formatting conventions:

* Branches related to a particular feature should be prefixed with `feature`,
i.e. for the BioSimSpace.Process module:

```bash
feature-process
```

When you are happy with your feature create a pull request so that the feature
can be merged into devel. Once the merge is successful, please delete the
redundant feature branch.

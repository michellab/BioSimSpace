# BioSimSpace

Code and resources for the EPSRC [BioSimSpace](https://biosimspace.org) project.

To install into a [Sire](https://github.com/michellab/Sire) package, use:

```bash
python setup.py install
```

Note that this assumes that Sire's `bin` directory, `$HOME/sire.app/bin`,
has been added to your path. Alternatively, run:

```bash
$HOME/sire.app/bin/python setup.py install
```
We recommend using the `devel` branch of [Sire](https://github.com/michellab/Sire)
to make sure that you have access to the latest features and bug-fixes.

## Developing

Please create a feature branch for development work. We use the following
formatting conventions:

* Branches related to a particular feature should be prefixed with `feature`,
i.e. for the BioSimSpace.Sample module:

```bash
feature-sample
```

* If you wish to use sub-branches for a feature, they can be identified by
using a backslash, e.g.:

```bash
feature-sample/namd
feature-sample/amber
```

When you are happy with your feature, merge sub-branches into the main feature
branch then create a pull request so that the feature can be merged into devel.
Once the merge is successful, please delete the redundant feature branch.

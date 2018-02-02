# Demo

This directory contains example scripts that illustrate various ways of using
BioSimSpace. If you have already installed BioSimSpace into a Sire package then
you can run any of the scripts directly, i.e.

```bash
$HOME/sire.app/bin/python amber.py
```

Alternatively, set the `PYTHONPATH` environment so that the python interpreter
can locate BioSimSpace:

```bash
PYTHONPATH=../python $HOME/sire.app/bin/python amber.py
```

The interactive demos require [Jupyter](http://jupyter.org),
[NGLview](https://github.com/arose/nglview), and [py3Dmol](https://pypi.python.org/pypi/py3Dmol).
These aren't currently bundled with the default Sire installation, so you will
need to install them into your package as follows:

```bash
$HOME/sire.app/bin/pip install jupyter nglview py3Dmol
```

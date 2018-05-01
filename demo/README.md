# Demo

This directory contains example scripts that illustrate various ways of using
BioSimSpace. If you have already installed BioSimSpace into a Sire package then
you can run any of the scripts directly, i.e.

```bash
# Run an example AMBER workflow.
$HOME/sire.app/bin/python amber.py

# Run an example NAMD workflow.
$HOME/sire.app/bin/python namd.py

# Run a generic minimisation node using AMBER input files.
$HOME/sire.app/bin/python minimisation.py --steps=1000 --files amber/ala/*

# Run a generic minimisation node using NAMD input files.
$HOME/sire.app/bin/python minimisation.py --steps=1000 --files namd/ala*/*
```

Alternatively, set the `PYTHONPATH` environment so that the python interpreter
can locate BioSimSpace:

```bash
PYTHONPATH=../python $HOME/sire.app/bin/python amber.py
```
## Notebooks

We also provide several [Jupyter](http://jupyter.org) notebooks that show how
to work with BioSimSpace in an interactive environment. If you used the
regular [setup.py](../python/setup.py) installation, then you should be good
to go. If not, the following describes installation instructions for the
required dependencies.

To enable plotting functionality in the example notebooks you will also need
to install [Matplotlib](https://matplotlib.org):

```bash
$HOME/sire.app/bin/pip install matplotlib
```

Molecular visualisation requires [NGLView](https://github.com/arose/nglview),
which can be a bit of a pain to configure correctly.

First install Jupyter:
```bash
$HOME/sire.app/bin/pip install jupyter
```

Now install NGLView, overwriting any dependencies that were already installed above:
```bash
$HOME/sire.app/bin/pip install --no-cache-dir --ignore-installed nglview
```

We now need to enable the notebook extension:
```bash
$HOME/sire.app/bin/jupyter-nbextension enable nglview --py --sys-prefix
```

Versions of NGLView earlier than 1.1.2 have an issues that causes intermittent
problems when reading molecular structures from file. This can be fixed by
running the provided helper script:

```bash
./nglview_patch.sh $HOME/sire.app
```

To find out what version of NGLView you have installed, run:

```bash
$HOME/sire.app/bin/pip show nglview
```

In a similar manner we need to install and activate the Jupyter
[fileupload](https://pypi.python.org/pypi/fileupload) widget.

```bash
$HOME/sire.app/bin/pip install fileupload
$HOME/sire.app/bin/jupyter-nbextension install fileupload --py --sys-prefix
$HOME/sire.app/bin/jupyter-nbextension enable fileupload --py --sys-prefix
```

Finally, you may need to clean any existing Jupyter configuration folders, which
will likely be present if you already have Jupyter installed on your machine, i.e.
installed outside of your `sire.app`. For simplicity we'll delete these files,
although you may want to back them up:

```bash
rm -r ~/.jupyter ~/.local/share/jupyter
```

(In future, we probably want to set different `JUPYTER_CONFIG_DIR` and
`JUPYTER_PATH` environment variables for BioSimSpace, or provide a custom
`jupyter` executable that sets the variables before launching Jupyter.)

Phew!

You should now be able to launch Jupyter:

```bash
$HOME/sire.app/bin/jupyter notebook

or

PYTHONPATH=../python $HOME/sire.app/bin/jupyter notebook
```

This will open a tab in your browser, where you can choose a specific demo
notebook, such as `amber.ipynb`. If you find that you have issues visualising
systems with a large number of molecules, e.g. solvated systems, then you
might need to launch Jupyter as follows:

```bash
$HOME/sire.app/bin/jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000
```

Remember to shutdown the python kernel before closing any notebook tabs so that
any BioSimSpace processes that are running can be safely killed and cleaned up.

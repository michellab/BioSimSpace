# Demo

This directory contains example scripts that illustrate various ways of using
BioSimSpace. If you have already installed BioSimSpace into a Sire package then
you can run any of the scripts directly, i.e.

```bash
# Run an example AMBER workflow.
$HOME/sire.app/bin/python amber.py

# Run an example NAMD workflow.
$HOME/sire.app/bin/python namd.py
```

Alternatively, set the `PYTHONPATH` environment so that the python interpreter
can locate BioSimSpace:

```bash
PYTHONPATH=../python $HOME/sire.app/bin/python amber.py
```
## Notebooks

We also provide several [Jupyter](http://jupyter.org) notebooks that show how
to work with BioSimSpace in an interactive environment.

To enable plotting functionality in the example notebooks you will also need
to install [Matplotlib](https://matplotlib.org). This will already have been
installed if you have run any of the BioSimSpace demo scripts. If not, run
the following command:

```bash
$HOME/sire.app/bin/pip install matplotlib
```

Molecular visualisation requires [NGLview](https://github.com/arose/nglview),
which is a bit of a pain to configure. It is not possible to use the regular
`try_import` function from `Sire`, since `conda` will install the incorrect
version. At present, manual installation via `pip` is required.

First install Jupyter:
```bash
$HOME/sire.app/bin/pip install jupyter
```

Now install NGLview, overwriting any dependencies that were already installed above:
```bash
$HOME/sire.app/bin/pip --no-cache-dir --ignore-installed install nglview
```

We now need to enable the notebook extension:
```bash
$HOME/sire.app/bin/jupyter-nbextension enable nglview --py --sys-prefix
```

Unfortunately, the version of NGLview provided by `pip` (version 1.0, at this
time) contains a bug that causes intermittent issues when reading molecular
structures from file. This can be fixed by running the provided helper
script:

```bash
./nglview_patch.sh $HOME/sire.app
```

Finally, we need to clean any existing Jupyter configuration folders, which
may be present if you already have Jupyter installed on your machine, i.e.
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

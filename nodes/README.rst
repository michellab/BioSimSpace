Nodes
=====

A collection of example BioSimSpace nodes.

* *playground*: Experimental nodes: here be dragons!

For developers
--------------

Please install `nb-clean <https://pypi.org/project/nb-clean>`_ to ensure that
all Jupyter notebooks are cleaned of cell execution counts, metadata, and
outputs. This means that any commits will correspond to actual changes
in the notebook code and that users are presented with a fresh notebook
when it is launched.

Please also add ``***NO_CI***`` to commit messages to avoid unnecessarily
running the Azure pipelines build, e.g. rebuilding the BioSimSpace binary,
updating the website, etc.

BioSimSpace release instructions
*********************************

The following instructions describe how to create a new BioSimSpace release:

Step 1
======

Update the `CHANGELOG <https://github.com/openbiosim/biosimspace/blob/devel/doc/source/changelog.rst>`_
file with a summary of the changes for this relase. Feel free to link to
`GitHub issues <https://github.com/openbiosim/biosimspace/issues>`_ where relevant
and give credit for specific contributions.

Step 2
======

When you're happy, tag the commit that you want to be associated with the
release. The following will tag the latest commit:

.. code-block:: bash

    git tag -a 2023.1.0 -m "Tagging the 2023.1.0 release of BioSimSpace."

Step 3
======

Push the commit and tag to the ``devel`` branch on the remote:

.. code-block:: bash

    git push origin devel --follow-tags

This will trigger a new Azure Pipelines build which will create binaries
and Conda packages for the release. If you make a mistake and want to move
the tag to a later commit, simply delete the tag from the remote:

.. code-block:: bash

    git push origin :refs/tags/2023.1.0

Next, delete the Conda release package from the `Anaconda Cloud <https://anaconda.org/openbiosim/biosimspace/files>`_.
You can then move the tag to the latest commit:

.. code-block:: bash

    git tag -fa 2023.1.0

Finally, push the new commit and updated tag:

.. code-block:: bash

    git push origin devel --follow-tags

Step 4
======

Create a `GitHub release <https://github.com/openbiosim/biosimspace/releases>`_.
When drafting the release, simply choose the tag that you have created.

Step 5
======

Create a `pull request <https://github.com/openbiosim/biosimspace/pulls>`_ to
merge ``devel`` into the ``master`` branch.

That's it!

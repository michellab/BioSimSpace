BioSimSpace release instructions
*********************************

The following instructions describe how to create a new BioSimSpace release:

Step 1
======

Update the `CHANGELOG <https://github.com/michellab/BioSimSpace/blob/devel/doc/source/changelog.rst>`_
file with a summary of the changes for this relase. Feel free to link to
`GitHub issues <https://github.com/michellab/BioSimSpace/issues>`_ where relevant
and give credit for specific contributions.

Step 2
======

When you're happy, tag the commit that you want to be associated with the
release. The following will tag the latest commit:

.. code-block:: bash

    git tag -a 2019.1.0 -m "Tagging the 2019.1.0 release of BioSimSpace."

Step 3
======

Push the commit and tag to the ``devel`` branch on the remote:

.. code-block:: bash

    git push origin devel --follow-tags

This will trigger a new Azure Pipelines build which will create binaries
and Conda packages for the release. If you make a mistake and want to move
the tag to a later commit, simply delete the tag from the remote:

.. code-block:: bash

    git push origin :refs/tags/2019.1.0

Next, delete the Conda release package from the `Anaconda Cloud <https://anaconda.org/michellab/biosimspace/files>`_.
You can then move the tag to the latest commit:

.. code-block:: bash

    git tag -fa 2019.1.0

Finally, push the new commit and updated tag:

.. code-block:: bash

    git push origin devel --follow-tags

Step 4
======

Once the build has finished you can log into the `Oracle Cloud <https://cloud.oracle.com/home>`__
and create pre-authenticated URLs for the new release binaries. (These files
are currently located in the ``software_releases`` compartment of the ``Object Storage``
menu. When generating a download URL make sure to choose a sensible expiry
date. Copy the URL to your clipboard and add download links to the
`install <https://github.com/michellab/BioSimSpace/blob/devel/doc/source/install.rst>`_
page in the website documentation. (Note that the website won't be updated
until you next trigger a development build. Alternatively, you can move the
tag to this commit and re-push following the instructions above. The binaries
will be overwritten, but the URLs won't change.)

Step 5
======

Create a `GitHub release <https://github.com/michellab/BioSimSpace/releases>`_.
When drafting the release, simply choose the tag that you have created.

Step 6
======

Create a `pull request <https://github.com/michellab/BioSimSpace/pulls>`_ to
merge ``devel`` into the ``master`` branch.

That's it!

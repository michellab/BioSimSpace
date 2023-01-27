==================
Packaging releases
==================

:mod:`BioSimSpace` is now (almost) fully tested and deployed using GitHub actions.
The development process should be;

* New features are developed on forks or feature branches, called ``feature-{feature}``,
  either in the `main BioSimSpace repository <https://github.com/openbiosim/biosimspace>`__
  for authorised developers, or in personal forks for
  new developers.
* Bug fixes or issue fixes are developed on fix branches, called
  ``fix-issue-{number}`` (again in either the main repository or forks).
* Pull requests are issued from these branches to ``devel``. All merge conflicts
  must be fixed in the branch and all tests must pass before the pull
  request can be merged into ``devel``. NOTE THAT ONLY AUTHORISED
  DEVELOPERS CAN ACCEPT THE PULL REQUEST. Authorised developers will
  review the pull request as quickly as they can. They will be greatly
  helped if the feature is accompanied with tests, examples and/or tutorial
  instructions.
* Note that we encourage that new functionality in a pull request is
  documented, e.g. by adding to the tutorials or writing a detailed
  description for the website.

The result of this is that ``devel`` should contain the fully-working,
documented and tested code, and most up-to-date version of :mod:`BioSimSpace`.
However, this version should not be used for production runs.

.. note::

  The group of developers authorised to have access to the
  `main BioSimSpace repository <https://github.com/openbiosim/biosimspace>`__
  and to accept pull requests is not fixed,
  and will evolve over time. If you wish to join this group then
  please complete the tutorial and then demostrate your commitment
  by submitting good issues and pull requests from
  a personal fork of the repository. Please get in touch if you find
  this difficult, or follow
  `this workshop <https://chryswoods.com/introducing_git>`__
  and `this workshop <https://chryswoods.com/git_collaboration>`__ if you need
  to learn how to use Git, GitHub, feature branching, merging, pull
  requests etc.

Defining a release
------------------

We will release :mod:`BioSimSpace` regularly. Releases aim to be backwards
compatible and capable of being used for production runs, at least for
the functionality that is fully described in the tutorial.

.. note::

  It is the job of the release managers (currently
  `lohedges <https://github.com/lohedges>`__ and
  `chryswoods <https://github.com/chryswoods>`__) to decide when it is time
  to create a new release. If you are interested in helping join the release
  management group then please feel free to get in touch.

Creating a release
------------------

There are a number of stages to go through to create a release:

1. Make sure that all changes have been merged into the ``devel`` branch,
   and the GitHub Action has run fully, building :mod:`BioSimSpace` on all
   supported platforms, running all the unit tests correctly, and
   building and uploading the conda packages to the ``dev`` channel
   on `anaconda.org <https://anaconda.org/openbiosim/biosimspace/files>`__.

2. On your computer, checkout the ``main`` branch and then pull in the
   changes from ``devel``, e.g. using the command ``git pull origin main``.
   There shouldn't be any conflicts. If there are, resolve them so that
   ``main`` is identical to ``devel``.

3. Now update ``doc/source/changelog.rst`` with the changelog since
   the last release. Start the entry with a title that (will) link
   to the changes on GitHub, e.g.
   `https://github.com/openbiosim/biosimspace/compare/2022.2.1...2023.1.0 <https://github.com/openbiosim/biosimspace/compare/2022.2.1...2023.1.0>`__
   links to the changes between the ``2022.2.1`` and ``2023.1.0`` releases.
   You can change the version numbers in this URL to look at the
   difference for your versions. Follow a similar format for changes
   as already exist in this file. Try to link to tutorials that
   describe new functionality if available.

4. You can now commit your changes and push them to GitHub
   (e.g. ``git commit -a`` then ``git push``). GitHub Actions will
   only be triggered on pushes to ``devel``, so this won't do anything (yet).

5. Now create a git tag for this release. You can use the command with
   the format ``git tag -a {VERSION} -m "{VERSION} release"``, e.g.
   ``git tag -a 2023.1.1 -m "2023.1.1 release"`` would be the
   tag for the ``2023.1.1`` release.

6. Push this tag to GitHub, e.g. ``git push origin tag``, where
   ``tag`` is the new tag that you've created, e.g. ``2023.1.1``.

7. Now we are finally ready to build the packages. Do this by checking
   out the ``devel`` branch (``git checkout devel``) and then pulling
   the changes from ``main`` into ``devel`` (``git pull origin main``).
   This will create a merge commit. Accept the commit and then
   push this back up to GitHub (``git push``). This will trigger the
   GitHub actions that will automatically build and upload all(*)
   the conda packages.

8. (*) GitHub Actions don't currently build the ARM64/aarch64 packages.
   These have to be built and uploaded manually. On a MacOS/M1 or
   Linux/aarch64 computer you should create build environments for
   the Python versions that :mod:`BioSimSpace` should support. Activate
   this environment, and then checkout the ``main`` branch, run
   ``python actions/update_recipe.py`` and then run ``conda-build``
   via the command ``conda mambabuild -c conda-forge -c openbiosim/label/dev recipes/biosimspace``.
   This will result in a conda package in the ``conda-bld`` directory
   in the root directory of your conda environment. You then need
   to upload these packages to anaconda.org, e.g. via the command
   ``anaconda --token {PASSWORD} upload --user openbiosim --label main --label dev --force {/path/to/biosimspace-packages}``
   (modified as appropriate to include the anaconda password and the
   path to the built conda package).

9. On GitHub, you can now create a release by using the
   `Draft a New Release <https://github.com/OpenBioSim/biosimspace/releases/new>`__
   link. Choose the version number for your release from the
   tag you created earlier. The text should be simple,
   e.g. titled ``BioSimSpace {VERSION}``, with the body
   ``This is the {VERSION} release of BioSimSpace.``.

10. Next you should build the docker images for this release.
    Do this by following the instructions in the
    `containers repository <https://github.com/OpenBioSim/containers/blob/main/biosimspace/README.md>`__.
    You should make sure to run the extra command listed there
    to tag the container with the version number you used earlier.

11. Finally(!) you can now update the website. To do this, follow the
    instructions in the `website repository <https://github.com/OpenBioSim/biosimspace_website/blob/main/README.md>`__.

12. Bonus! Follow the instructions in the
    `containers repository <https://github.com/OpenBioSim/containers>`__
    to build the notebook container image and instruct
    `try.openbiosim.org <https://try.openbiosim.org>`__ to update
    and use that image.

13. Super-bonus! If you have time, please write a short news item piece
    that can be added to the `openbiosim website <https://openbiosim.org>`__
    to announce this new release.

(We are in the process of automating many of the above steps, so hope
that this process will become much easier in the future.)

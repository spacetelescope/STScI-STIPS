******************************
Developer Release Instructions
******************************

This page offers a complete walkthrough on how to prepare and release a new
version of STIPS on GitHub and PyPI, as well as how to sync the ReadTheDocs
documentation with the new release.

Pre-merge instructions
======================

Pull request checklist
----------------------

In addition to addressing any byproducts of the code review, you should be able
to answer "yes" to all of the following questions before merging a pull request
intended for a release:

* Have ``setup.cfg`` and ``stips/__init__.py`` been updated to reflect the release's new version number?

* If any dependencies have changed, have you updated them in all applicable files? These include:

    * ``setup.cfg``
    * ``environment.yml``
    * ``environment_dev.yml``
    * ``docs/installation.rst``
    * ``ret_data/retrieve_stips_data.py``

* Are you able to install the new version of STIPS locally with ``pip`` and ``conda``?

* Can you run the example notebooks in the ``notebooks/`` directory without encountering errors?

* Has ``CHANGES.rst`` been updated to include the current release and its associated updates at the top?

* Are all CI tests in the release's pull request passing?

* Has at least one developer given the pull request an official approval on GitHub?

Once you've completed the checklist, merge the pull request and begin the
release procedure in earnest.

Post-merge instructions
=======================

Tagging the merge commit
------------------------

Apply a tag to the merge commit, substituting in the proper version number.

.. code-block:: text

    git tag -s vX.Y.Z -m "Release vX.Y.Z"

.. note::
   The ``-s`` option allows you to sign the tag using the GNU Privacy Guard (GPG).
   While it's not strictly necessary, it's a good security measure.
   Follow GitHub's documentation for `creating a GPG key <https://docs.github.com/en/authentication/managing-commit-signature-verification/generating-a-new-gpg-key>`_
   and `adding it to your account <https://docs.github.com/en/authentication/managing-commit-signature-verification/adding-a-gpg-key-to-your-github-account>`_
   if you don't have one.

Next, push the tag online. (If your remote that points to ``spacetelescope``'s
online version of the repository is not named ``origin``, be sure to substitute
its name in below.)

.. code-block:: text

    git push origin vX.Y.Z

Release creation (on GitHub)
----------------------------

After the tests triggered by pushing the new tag succeed, follow these steps to
manually create a new release in the ``spacetelescope/STScI-STIPS`` GitHub repository.

#. Click "Releases" on main repository page or `follow this link <https://github.com/spacetelescope/STScI-STIPS/releases>`_.
#. Press the "Draft a new release" button above the list of releases.
#. Press "Choose a tag" and select the one you just pushed online. The target branch should already be selected as ``main``.
#. Press "Generate release notes" to list the commits included in this release.
#. In the "Release title" textbox, type ``Release X.Y.Z``, substituting in the proper version number.
#. Above the list of commits in the larger textbox for comments, write a single-sentence summary of the release's major updates.
#. Press "Publish" to complete this section.

Release creation (for PyPI)
---------------------------

Continue by building wheel and source distribution files to publish to PyPI so
users can ``pip`` install STIPS without cloning the repository. You will need a
PyPI account and maintainer status on `the STIPS project <https://pypi.org/project/stips/>`_
if you don't already have it there.

.. note::

  It's highly recommended that you clone a new copy of the ``STScI-STIPS``
  repository into a new directory on your system and work from there while
  completing this section.


Make sure you've installed the necessary packages for this procedure:

.. code-block:: text

    pip install --user --upgrade setuptools wheel twine

Then, identify the largest subdirectories in your STIPS directory:

.. code-block:: text

    du -h

Due to past commits that included large files, STIPS' ``.git/`` subdirectory
will likely be the largest at over 200 MB. Even with compression, working with
the repository as is generates distribution files that surpass PyPI's file size
limit for uploads.

Since the files in ``.git/`` aren't necessary for users to download as part of
a pip installation, delete them before building anything:

.. code-block:: text

    rm -r .git

Now, create the wheel and source distribution files:

.. code-block:: text

    python setup.py sdist bdist_wheel

The resulting wheel file and tarball are located in the ``dist`` directory.

Finally, upload them to PyPI:

.. code-block:: text

    python setup.py sdist bdist_wheel

Note that PyPI uploads now require an API token.
`Refer to their instructions <https://pypi.org/help/#apitoken>`_ if you haven't
yet set one up.

While the official release is now complete, keep reading for instructions on
updating the documentation on ReadTheDocs.

Documentation
=============

Navigate to `the active STIPS ReadTheDocs page <https://readthedocs.org/projects/stips/>`_.
(Note that the ReadTheDocs project name is ``stips``, matching the package name
but not the GitHub repository name.) Verify that new builds of ``latest`` and
``vX.Y.Z`` have been run successfully. If not, build them manually under the
"Build a version" header. If you followed earlier instructions, the release
commit will be the repository's current latest commit, so both versions should
be identical for the moment.

The versions of the documentation that should be visible to the public and
marked as "Active" on their "Edit" pages are ``latest`` and the new release,
``vX.Y.Z``. In the future, we may make past versions visible, too.

Troubleshooting: webhooks
-------------------------

ReadTheDocs is connected to the ``spacetelescope/STScI-STIPS`` GitHub repository
via a webhook, which can be found on `the repository's "Webhooks" page <https://github.com/spacetelescope/STScI-STIPS/settings/hooks>`_.
If for any reason the link happens to break, the webhook can be re-linked by
creating a new secret for the GitHub incoming webhook on `the ReadTheDocs "Integrations" page <https://readthedocs.org/dashboard/stips/integrations>`_
and then pasting it to the webhook on GitHub with the matching Payload URL.

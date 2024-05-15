************
Installation
************

STIPS is a simulation tool that depends on other modules such as PSF and exposure time calculators.  These underlying submodules need to be
installed for STIPS to function properly along with their supporting datasets.  There are multiple options for installation and they are listed
in this section along with instructions.

STIPS Requirements
##################

* ``pandeia>=3.1``: Exposure time calculator.

* ``webbpsf>=1.1.1``: Nancy Grace Roman PSF calculator. STIPS also requires that ``poppy``, a
  support package used by WebbPSF, have version ``>=1.0.3``.

* ``astropy``: STIPS uses Astropy in order to:

    * Read and write FITS files.

    * Read and write ASCII tables (specifically in the IPAC format).

    * Generate Sersic profile models (if any are in the generated scene).

* ``montage_wrapper``: STIPS uses ``montage`` to generate mosaics. It is
  only imported if STIPS is asked to generate a multi-detector image.

* ``numpy``: STIPS uses ``numpy`` extensively for almost everything that it does.

* ``photutils``: STIPS uses ``photutils`` to determine the flux inside the half-light radius
  in generated Sersic profiles.

* ``synphot>=1.1.1`` and ``stsynphot>=1.1.0``: STIPS uses ``synphot`` and
  ``stsynphot`` to generate bandpasses, count rates, and zero points. Note that
  the reference data must also be downloaded, as described below in :ref:`downloading-required-ref-data`.

* ``scipy``: STIPS uses SciPy to manipulate its internal images (zoom and rotate).

* ``esutil``: Used for retrieving data from sqlite databases in the form of ``numpy`` arrays.

.. warning::
   ``esutil`` is not installed by the STIPS ``setup.py`` file because its pip
   installation has errors. After installing STIPS, you must run ``pip install esutil --no-cache-dir``
   to have ``esutil`` available.

.. note::
   ``esutil`` is only needed if you are using ``stips.star_generator``.

Finally, STIPS requires a set of data files whose location is marked by setting the
environment variable ``stips_data``, which will be installed as part of these instructions.

Installing Using Conda and Source Code
######################################

STIPS can be installed using the source code and a Conda environment file.
If you do not have Anaconda or Miniconda installed, please visit the
`Anaconda docs <https://docs.anaconda.com/anaconda/install/>`_ for installation instructions.
We have included a Conda environment file for easily installing or updating Conda packages
to meet STIPS requirements.  Please follow the steps below to install STIPS:

.. _installing-as-a-user:

Installing as a User
********************

#. You will need to clone the STIPS source code from the
   `spacetelescope/STScI-STIPS <https://github.com/spacetelescope/STScI-STIPS.git>`_
   repository. ``cd`` into the directory where you would like to store the
   source code and run::

        git clone https://github.com/spacetelescope/STScI-STIPS.git

        cd STScI-STIPS

#. The environment file can be used in two ways:

   * To create a new Conda environment named ``stips``::

        conda env create -f environment.yml
        conda activate stips


   * Or, to install to or update an existing Conda environment::

        conda env update --name EXISTING-ENV --file environment.yml

Installing as a Developer
*************************

#. This step is identical to the first step of :ref:`installing-as-a-user`.

#. Follow the second step of :ref:`installing-as-a-user` but using the
   ``environment_dev.yml`` file instead of ``environment.yml``.

.. _downloading-required-ref-data:

Downloading Required Reference Data
************************************

STIPS, Pandeia, and WebbPSF need the reference datasets.
You will need to download the data and add them to your environmental path.

1. Add the following paths to your bash environmental path. It is recommended that you add the path to your ``.bash_profile`` file:

.. code-block:: text

	export stips_data="<absolute_path_to_this_folder>/ref_data/stips_data"
	export WEBBPSF_PATH="<absolute_path_to_this_folder>/ref_data/webbpsf-data"
	export PYSYN_CDBS="<absolute_path_to_this_folder>/ref_data/grp/redcat/trds"
	export pandeia_refdata="<absolute_path_to_this_folder>/ref_data/pandeia_data-x.x.x_roman"

.. note::

  Make sure that you have the correct version number for ``pandeia_refdata`` (replace the "x.x.x").

2. ``cd`` into the ``ref_data`` directory in your ``STScI-STIPS`` clone.

3. Run the following code (after ensuring your ``stips`` Conda environment is active)::

		python retrieve_stips_data.py


Testing Installation
*********************

To test if all the required files have been installed, please import STIPS in Python::

    bash-3.2$ python
    Python 3.11.9 | packaged by conda-forge | (main, Apr 19 2024, 18:45:13)
    [Clang 16.0.6 ] on darwin
    Type "help", "copyright", "credits" or "license" for more information.

    >>> import stips

    >>> print(stips.__env__report__)

You should receive an output of the following form:

.. code-block:: text

  STIPS Version x.y.z with Data Version x.y.z at /Some/Path/To/stips_data

  STIPS Grid Generated with x.y.z

  Pandeia version a.b.c with Data Version a.b.c. at /Some/Path/To/pandeia_refdata

  Webbpsf Version d.e.f with Data Version d.e.f at /Some/Path/To/webbpsf_data_path


Ignore the following warning message if it appears:

.. code-block:: text

  WARNING: stips_data environment variable not found. Falling back on local STIPS data.

************
Installation
************

STIPS is a simulation tool that depends on other modules such as PSF and exposure time calculators.  These underlying submodules need to be installed for STIPS to function properly along with their supporting datasets.  There are multiple options for installation and they are listed in this section along with instructions.

STIPS Requirements
##################

* `Pandeia`: Exposure time calculator.
* `WebbPSF`: James Webb and Nancy Grace Roman PSF calculator.
* `astropy`: STIPS uses astropy in order to:

	- Read and write FITS files.
	- Read and write ASCII tables (specifically in the IPAC format).
	- Generate Sersic profile models (if any are in the generated scene).

* `esutil`: Used for retrieving data from sqlite databases in the form of numpy arrays.
* `montage_wrapper`: STIPS uses montage to generate mosaics. It is only imported if
  STIPS is asked to generate a multi-detector image.
* `numpy`: STIPS uses numpy extensively for almost everything that it does.
* `photutils`: STIPS uses photutils to determine the flux inside the half-light radius
  in generated Sersic profiles.
* `synphot` and `stsynphot`: STIPS uses synphot and stsynphot to generate
  bandpasses, count rates, and zero points. Note that the reference data must
  also be downloaded, as described below in "Doanloading Required Data".
* `scipy`: STIPS uses scipy to manipulate its internal images (zoom and rotate).

Finally, STIPS requires a set of data files whose location is marked by setting the environment variable `stips_data`, which will be installed as part of these instructions.

Installing Using Conda and Source
##################################

STIPS can be installed using the source code and a Conda environment file.
If you do not have anaconda or miniconda installed, please visit the `anaconda docs <https://docs.anaconda.com/anaconda/install/>`_ for installation instructions.  We have included a Conda environment file for easily installing or updating Conda packages to meet STIPS requirements.  Please follow the steps below to install STIPS:

Installing
**********

1. You will need to clone the STIPS source code from the `spacetelescope/STScI-STIPS <https://github.com/spacetelescope/STScI-STIPS.git>`_ repository.  `cd` into the directory you would like to store the source code and run::

    git clone https://github.com/spacetelescope/STScI-STIPS.git

    cd STScI-STIPS

2. The environment file can be used in two ways:

    a. To create a new Conda environment named `stips` run::

        conda env create -f environment.yml

        conda activate stips


    b. To install to or update an existing (currently active) Conda environment::

        conda env update --file environment.yml


3. You can now install STIPS using the cloned source code as follows::

    python setup.py install


Downloading Required Reference Data
************************************

STIPS, Pandeia, and WebbPSF need the reference datasets.
You will need to download the data and add them to your environmental path

1. Add the following paths to your bash environmental path. It is recommended that you add the path to your `.bash_profile` file::

		export stips_data="<absolute_path_to_this_folder>/ref_data/stips_data"
		export WEBBPSF_PATH="<absolute_path_to_this_folder>/ref_data/webbpsf-data"
		export PYSYN_CDBS="<absolute_path_to_this_folder>/ref_data/grp/redcat/trds"
		export pandeia_refdata="<absolute_path_to_this_folder>/ref_data/pandeia_data-x.x.x_roman"

Make sure that you have the correct version number for `pandeia_refdata` (replace the "x.x.x").

2. `cd` into the "ref_data" directory in your STScI-STIPS clone

3. Run the following code (ensuring your stips environment is active)::

		python retrieve_stips_data.py


Testing Installation
*********************

To test if all the required files have been installed, please import STIPS in python::

    bash-3.2$ python
    Python 3.7.3 | packaged by conda-forge | (default, Dec  6 2019, 08:36:57)
    [Clang 9.0.0 (tags/RELEASE_900/final)] :: Anaconda, Inc. on darwin
    Type "help", "copyright", "credits" or "license" for more information.

    >>> import stips

		    print(stips.__env__report__)

You should receive an output of the following form::

		STIPS Version x.y.z with Data Version x.y.z at /Some/Path/To/stips_data

		STIPS Grid Generated with x.y.z

		Pandeia version a.b.c with Data Version a.b.c. at /Some/Path/To/pandeia_refdata

		Webbpsf Version d.e.f with Data Version d.e.f at /Some/Path/To/webbpsf_data_path

The following warning message can be ignored if it appears::

    WARNING: stips_data environment variable not found. Falling back on local STIPS data.

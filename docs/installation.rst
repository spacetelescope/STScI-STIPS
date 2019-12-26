************
Installation
************

STIPS is a simulation tool that depends on other modules such as PSF and exposure time calculators.
These underlying submodules need to be installed for STIPS to function properly along with their supporting datasets.
There are multiple options for installation and they are listed in this section along with instructions.

STIPS Requirements
##################

* `Pandeia`: Exposure time calculator.
* `WebbPSF`: James Webb and WFIRST PSF calculator.
* `astropy`: STIPS uses astropy in order to

	- Read and write FITS files
	- Read and write ASCII tables (specifically in the IPAC format)
	- Generate Sersic profile models (if any are in the generated scene)

* `esutil`: Used for retrieving data from sqlite databases in the form of numpy arrays
* `montage_wrapper`: STIPS uses montage to generate mosaics. It is only imported if
  STIPS is asked to generate a multi-detector image.
* `numpy`: STIPS uses numpy extensively for almost everything that it does
* `photutils`: STIPS uses photutils to determine the flux inside the half-light radius
  in generated Sersic profiles
* `pysynphot`: STIPS uses pysynphot to generate bandpasses, count rates, and
  zero points. Note that pysynphot's data files (also known as the CDBS data tree) must also be
  installed and available as indicated in pysynphot's documentation.
* `scipy`: STIPS uses scipy to manipulate its internal images (zoom and rotate)

Finally, STIPS requires a set of data files whose location is marked by setting the environment
variable `stips_data`. Currently these files are available as part of the STSCI-STIPS-UI GitHub
project, but they should eventually be made available as a (versioned) direct download.

Installing Using Conda and Source
##################################

STIPS can be installed using the source code and a Conda environment file.
If you do not have anaconda or miniconda installed, please visit the ` astroconda docs <https://astroconda.readthedocs.io/en/latest/getting_started.html#getting-started-jump>`_ for installation instructions.
We have included a Conda environment file for easily installing or updating Conda packages to meet STIPS requirements.
Please follow the steps below to install STIPS:

Installing
**********

1. You will need to clone the STIPS source code from the `spacetelescope/STScI-STIPS <https://github.com/spacetelescope/STScI-STIPS.git>`_ repository.`cd` into the directory you would like to store the source code and run::

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


Downloading Required Data
*************************

Pandeia and WebbPSF need the reference datasets.
You will need to download the data and add them to your environmental path

1. Check if wget is installed by running::

    wget

    # If it is not installed (usually on Mac OS),
    # it can be downloaded using Conda

    conda install wget

2. `cd` into a directory you would like to store the data in and run the following commands::

    # PySynphot reference data
    wget -qO- http://ssb.stsci.edu/cdbs/tarfiles/synphot1.tar.gz | tar xvz
    wget -qO- http://ssb.stsci.edu/cdbs/tarfiles/synphot2.tar.gz | tar xvz
    wget -qO- http://ssb.stsci.edu/cdbs/tarfiles/synphot5.tar.gz | tar xvz


    # Pandeia reference data
    wget -qO- https://stsci.box.com/shared/static/5j506xzg9tem2l7ymaqzwqtxne7ts3sr.gz | tar -xvz

    # WebbPSF reference data
    wget -qO- https://stsci.box.com/shared/static/qcptcokkbx7fgi3c00w2732yezkxzb99.gz | tar xvz

3. Add the data paths to your bash environmental path. It is recommended that you add the path to your `.bashrc` file::

    export WEBBPSF_PATH=/<path_to_data_dir>/webbpsf-data
    export PYSYN_CDBS=/<path_to_data_dir>/grp/hst/cdbs
    export pandeia_refdata=/<path_to_data_dir>/pandeia_data-x.x_wfirst

Make sure that you have the correct version number for `pandeia_refdata` (replace the "x.x").


Testing Installation
*********************

To test if all the required files have been installed, please import STIPS in python::

    bash-3.2$ python
    Python 3.7.3 | packaged by conda-forge | (default, Dec  6 2019, 08:36:57)
    [Clang 9.0.0 (tags/RELEASE_900/final)] :: Anaconda, Inc. on darwin
    Type "help", "copyright", "credits" or "license" for more information.

    >>> import stips

The following warning message can be ignored if it appears::

    WARNING: stips_data environment variable not found. Falling back on local STIPS data.

************
Installation
************

STIPS is a simulation tool that depends on other modules such as PSF and exposure time calculators.
These underlying submodules need to be installed for STIPS to function properly along with their supporting datasets.
There are multiple options for installation and they are listed in this section along with instructions.

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
  STIPS is asked to generate a multi-detector montage image.
* `numpy`: STIPS uses numpy extensively for almost everything that it does.
* `photutils`: STIPS uses photutils to determine the flux inside the half-light radius
  in generated Sersic profiles.
* `synphot` and `stsynphot`: STIPS uses synphot and stsynphot to generate 
  bandpasses, count rates, and zero points. Note that the reference data must
  also be downloaded, as described below in "Doanloading Required Data".
* `scipy`: STIPS uses scipy to manipulate its internal images (zoom and rotate).

Finally, STIPS requires a set of data files whose location is marked by setting the environment
variable `stips_data`. The current version of the STIPS data is located on box and can be downloaded via the link below.

Downloading STIPS Data
#######################

STIPS needs data for reference and calibration. The latest version of the STIPS data can be downloaded as follows, from within
an interactive python interpreter::

    # import STIPS
    import stips
    
    # Download all required reference data:
    stips.DownloadReferenceData()


Installing Using Conda and Source
##################################

STIPS can be installed using the source code and a Conda environment file.
If you do not have anaconda or miniconda installed, please visit the `astroconda docs <https://astroconda.readthedocs.io/en/latest/getting_started.html>`_ for installation instructions.
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

STIPS will automatically be installed by the environment file.


Downloading Required Data
*************************

Pandeia and WebbPSF need the reference datasets. The appropriate versions of this 
reference data for the version of STIPS you are using may be downloaded, along with the
STIPS reference data, by the following convenience function, run from python::

    import stips
    stips.DownloadReferenceData()

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


Installing Using Docker
#######################

Installing
**********

1. Start by installing the free `Docker Community Edition <https://www.docker.com/community-edition>`_ locally.
This will make the `docker` command available in your terminal. Note that after installing docker,
you must open the application once for docker to be available from the command line.

2. You will need to clone the STIPS source code from the `spacetelescope/STScI-STIPS <https://github.com/spacetelescope/STScI-STIPS.git>`_ repository.
`cd` into the directory you would like to store the source code and run::

    git clone https://github.com/spacetelescope/STScI-STIPS.git

    cd STScI-STIPS

3. Run the docker build command::

    docker build -t stips .



Testing Installation
*********************

To test if the Docker image was built correctly you can `exec` into the image and try to import STIPS::

    # cd into STScI-STIPS
    $ docker build -t stips .

    # Create Docker Image
    $ docker create -t -i stips bash

        8293abe302b0c4f07a04282e811824d74681b77d0174148cc8af68078c098fa6

    # Start Docker Image
    $ docker start -a -i 8293abe302b0

    (stips) root@8293abe302b0:~# python
    Python 3.7.3 | packaged by conda-forge | (default, Jul  1 2019, 21:52:21)
    [GCC 7.3.0] :: Anaconda, Inc. on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import stips

STIPS Python Requirements
=========================
STIPS currently runs under python 2.7 (and has been developed with versions 2.7.3 through 2.7.12).
Eventually, STIPS will be made compatible with python 3.

STIPS requires the following STScI-supported python packages be installed and running correctly
(along with any data files that the packages themselves require, as can be found in their
documentation):

* Pandeia (tested with versions 1.0 and 1.1.1)
* Webbpsf (tested with versions between 0.4.0 and 0.6.0)

STIPS also uses the following (more general) python packages:

* astropy (version 1.3.2): STIPS uses astropy in order to
	* read and write FITS files
	* read and write ASCII tables (specifically in the IPAC format)
	* generate Sersic profile models (if any are in the generated scene)
* esutil (version 0.6.0): Used for retrieving data from sqlite databases in the form of numpy arrays
* montage_wrapper (version 0.9.8): STIPS uses montage to generate mosaics. It is only imported if
  STIPS is asked to generate a multi-detector image.
* numpy (version 1.12.1): STIPS uses numpy extensively for almost everything that it does
* photutils (version 0.3.2): STIPS uses photutils to determine the flux inside the half-light radius
  in generated Sersic profiles
* pysynphot (version 0.9.8.5): STIPS uses pysynphot to generate bandpasses, count rates, and
  zero points. Note that pysynphot's data files (also known as the CDBS data tree) must also be
  installed and available as indicated in pysynphot's documentation.
* scipy (version 0.18.1): STIPS uses scipy to manipulate its internal images (zoom and rotate)

Finally, STIPS requires a set of data files whose location is marked by setting the environment
variable `stips_data`. Currently these files are available as part of the STSCI-STIPS-UI github
project, but they should eventually be made available as a (versioned) direct download.


STIPS Python 3 Status
---------------------
STIPS was originally designed with and for python 2, but work has recently begun on converting it to
python 3. Currently, STIPS is known to be importable from python 3, and the STIPS observation
module has been used with an internal-format catalogue. Any issues encountered when using STIPS with
python 3 should be reported on
[STIPS Github Issue 16](https://github.com/spacetelescope/STScI-STIPS/issues/16).

To install, in bash:

Once you have a version of "conda" on your machine, Obtain a copy of cdbs and pandeia _data-1.0 (need link to these), and note their paths:

export PYSYN_CDBS="[path to cdbs]"

export pandeia_refdata="[path to pandeia_data-1.0]"

conda create -n forSTIPS python=2.7.12 astropy=1.3.2 numpy=1.12.1 scipy=0.18.1 photutils=0.3.2 pysynphot=0.9.8.5 webbpsf=0.6

source activate forSTIPS

echo $WEBBPSF_PATH #if this is defined, then the environment should be correct

pip install "esutil==0.6.0"

pip install "montage-wrapper==0.9.9"

pip install "jwst_backgrounds==1.1.1"

pip install "pandeia.engine==1.0"

#to install in /local/tmp/Work  Substitute directory as needed

cd /local/tmp/Work

git clone https://github.com/spacetelescope/STScI-STIPS-UI.git

export stips_data="/local/tmp/Work/STScI-STIPS-UI/sim_input/stips_data"

git clone https://github.com/spacetelescope/STScI-STIPS.git

cd STScI-STIPS

python setup.py install
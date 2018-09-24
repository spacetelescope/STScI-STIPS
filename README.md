# STScI-STIPS

## Table of Contents

* [Overview](#overview)
* [Why use STIPS?](#why-use-stips)
* [STIPS Python Requirements](#stips-python-requirements)
* [STIPS Examples](#stips-examples)

## Overview

STIPS is the Space Telescope Imaging Product Simulator. It is designed to create simulations of 
full-detector post-pipeline astronomical scenes for any telescope. Currently STIPS has modules for
WFC3 IR (F110W and F160W only), JWST (NIRCam Short, NIRCam Long, and MIRI), and WFIRST (WFI). STIPS
has the ability to add instrumental distortion (if available) as well as calibration residuals
(currently flatfield residuals, dark current residuals, and cosmic ray residuals). It automatically
includes Poisson noise and readout noise. It does not include instrument saturation effects. In
addition, STIPS has the ability to generate its own scenes, consisting of stellar populations and
background galaxies (implemented as Sersic profiles).

## Why use STIPS?

STIPS is intended for cases where an ETC (e.g. Pandeia) does not provide enough detector area (e.g.
testing photometry code, quick looks at dither patterns or multi-detector observations of a scene).
For JWST and WFIRST, it obtains its background count levels and instrumental throughput levels from
Pandeia internally, so it should produce output within 10% of output produced by Pandeia.

If extremely good instrumental accuracy is needed, STIPS is not the ideal choice. Instead, the
various instrument design teams have produced much more detailed simulators. STIPS is intended to
run reasonably quickly, and to make scene generation and observation as easy as possible.

## STIPS Python Requirements

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

## STIPS Examples

**NOTE:** If you do not have environment variables pointing to the location of your Pandeia data,
STIPS data, and Webbpsf data, you must set these environment variables using `os.environ` or other
equivalent method prior to running any of these example scripts.

* Creating a scene with a single stellar population and a single galaxy population, then observing
  it with NIRCam Short F115W:
  
  		from stips.scene_module import SceneModule
  		from stips.observation_module import ObservationModule
  	
 		scm = SceneModule()

  		stellar = {'n_stars': 50000, 
  				   'age_low': 1.0e12, 'age_high': 1.0e12, 
  				   'z_low': -2.0, 'z_high': -2.0,
  				   'imf': 'salpeter', 'alpha': -2.35,
  				   'binary_fraction': 0.1,
  				   'distribution': 'invpow', 'clustered': True,
  				   'radius': 100.0, 'radius_units': 'pc',
  				   'distance_low': 20.0, 'distance_high': 20.0,
  				   'offset_ra': 0.0, 'offset_dec': 0.0}
  		stellar_cat_file = scm.CreatePopulation(stellar)
  		
  		galaxy = {'n_gals': 1000,
  				  'z_low': 0.0, 'z_high': 1.0,
  				  'rad_low': 0.01, 'rad_high': 2.0,
  				  'sb_v_low': 30.0, 'sb_v_high': 25.0,
  				  'distribution': 'uniform', 'clustered': False,
  				  'radius': 200.0, 'radius_units': 'arcsec',
  				  'offset_ra': 0.0, 'offset_dec': 0.0}
  		galaxy_cat_file = scm.CreateGalaxies(galaxy)
  		
  		obs = {'instrument': 'NIRCamShort', 
  		       'filters': ['F115W'], 
  		       'detectors': 1,
  			   'distortion': False,
  			   'oversample': 5,
  			   'pupil_mask': '',
  			   'background': 'avg',
  			   'observations_id': 1,
  			   'exptime': 1000,
  			   'offsets': [{'offset_id': 1, 'offset_centre': False, 'offset_ra': 0.0, 'offset_dec': 0.0, 'offset_pa': 0.0}]}
  		obm = ObservationModule(obs)
  		obm.nextObservation()
  		output_stellar_catalogues = obm.addCatalogue(stellar_cat_file)
  		output_galaxy_catalogues = obm.addCatalogue(galaxy_cat_file)
  		psf_file = obm.addError()
  		fits_file, mosaic_file, params = obm.finalize(mosaic=False)
  
  In this case, the output FITS file will be in the variable `fits_file`, and the output catalogues
  (showing the actual count rate and position of the sources observed) will be in the variables
  `output_stellar_catalogues` and `output_galaxy_catalogues`.
* Creating a scene from an existing source catalogue `input_sources.txt`, and observing it with the
  WFIRST WFI "J129" filter, offset by 0.5 degrees in RA, and rotated by 27 degrees:
  
  		from stips.observation_module import ObservationModule
  	
  		obs = {'instrument': 'WFI', 
  		       'filters': ['J129'], 
  		       'detectors': 1,
  			   'distortion': False,
  			   'oversample': 5,
  			   'pupil_mask': '',
  			   'background': 'avg',
  			   'observations_id': 1,
  			   'exptime': 1000,
  			   'offsets': [{'offset_id': 1, 'offset_centre': False, 'offset_ra': 0.5, 'offset_dec': 0.0, 'offset_pa': 27.0}]}
		scene_general = {'ra': 256.274799731, 'dec': 22.6899695529, 'pa': 0.0, 'seed': 1}
  		obm = ObservationModule(obs, scene_general=scene_general)
  		obm.nextObservation()
  		source_count_catalogues = obm.addCatalogue('input_sources.txt')
  		psf_file = obm.addError()
  		fits_file, mosaic_file, params = obm.finalize(mosaic=False)
  
  In this case, the output catalogue(s) will show the actual applied count rates. Whether there is
  only one output catalogue or two depends on the input catalogue format.

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

STIPS Examples
===============
.. note::

    If you do not have environment variables pointing to the location of your Pandeia data,
    STIPS data, and Webbpsf data, you must set these environment variables using `os.environ` or other
    equivalent method prior to running any of these example scripts.


* Creating a scene from an existing source catalogue `input_sources.txt`, and observing it with the
  WFIRST WFI "F129" filter, offset by 0.5 degrees in RA, and rotated by 27 degrees:

.. code-block:: python

    from stips.scene_module import SceneModule
    from stips.observation_module import ObservationModule


    scm = SceneModule()

    stellar = {
                'n_stars': 50000,
                'age_low': 1.0e12, 'age_high': 1.0e12,
                'z_low': -2.0, 'z_high': -2.0,
                'imf': 'salpeter', 'alpha': -2.35,
                'binary_fraction': 0.1,
                'distribution': 'invpow', 'clustered': True,
                'radius': 100.0, 'radius_units': 'pc',
                'distance_low': 20.0, 'distance_high': 20.0,
                'offset_ra': 0.0, 'offset_dec': 0.0

              }

    stellar_cat_file = scm.CreatePopulation(stellar)

    galaxy = {
                'n_gals': 1000,
                'z_low': 0.0, 'z_high': 0.2,
                'rad_low': 0.01, 'rad_high': 2.0,
                'sb_v_low': 30.0, 'sb_v_high': 25.0,
                'distribution': 'uniform', 'clustered': False,
                'radius': 200.0, 'radius_units': 'arcsec',
                'offset_ra': 0.0, 'offset_dec': 0.0
            }

    galaxy_cat_file = scm.CreateGalaxies(galaxy)


    obs = {
            'instrument': 'WFI',
            'filters': ['F129'],
            'detectors': 1,
            'distortion': False,
            'oversample': 5,
            'pupil_mask': '',
            'background': 'avg',
            'observations_id': 1,
            'exptime': 1000,
            'offsets': [{'offset_id': 1, 'offset_centre': False, 'offset_ra': 0.5, 'offset_dec': 0.0, 'offset_pa': 27.0}]
          }

    obm = ObservationModule(obs)
    obm.nextObservation()

    output_stellar_catalogues = obm.addCatalogue(stellar_cat_file)
    output_galaxy_catalogues = obm.addCatalogue(galaxy_cat_file)

    psf_file = obm.addError()

    fits_file, mosaic_file, params = obm.finalize(mosaic=False)

In this case, the output catalogue(s) will show the actual applied count rates. Whether there is
only one output catalogue or two depends on the input catalogue format.


* Creating a scene with a single stellar population and a single galaxy population, then observing
  it with NIRCam Short F115W:

.. code-block:: python

    from stips.scene_module import SceneModule
    from stips.observation_module import ObservationModule

    scm = SceneModule()

    stellar = {
               'n_stars': 50000,
               'age_low': 1.0e12, 'age_high': 1.0e12,
               'z_low': -2.0, 'z_high': -2.0,
               'imf': 'salpeter', 'alpha': -2.35,
               'binary_fraction': 0.1,
               'distribution': 'invpow', 'clustered': True,
               'radius': 100.0, 'radius_units': 'pc',
               'distance_low': 20.0, 'distance_high': 20.0,
               'offset_ra': 0.0, 'offset_dec': 0.0
               }

    stellar_cat_file = scm.CreatePopulation(stellar)

    galaxy = {
              'n_gals': 1000,
              'z_low': 0.0, 'z_high': 1.0,
              'rad_low': 0.01, 'rad_high': 2.0,
              'sb_v_low': 30.0, 'sb_v_high': 25.0,
              'distribution': 'uniform', 'clustered': False,
              'radius': 200.0, 'radius_units': 'arcsec',
              'offset_ra': 0.0, 'offset_dec': 0.0
              }

    galaxy_cat_file = scm.CreateGalaxies(galaxy)

    obs = {
            'instrument': 'NIRCamShort',
            'filters': ['F115W'],
            'detectors': 1,
            'distortion': False,
            'oversample': 5,
            'pupil_mask': '',
            'background': 'avg',
            'observations_id': 1,
            'exptime': 1000,
            'offsets': [{'offset_id': 1, 'offset_centre': False, 'offset_ra': 0.0, 'offset_dec': 0.0, 'offset_pa': 0.0}]
           }

    obm = ObservationModule(obs)
    obm.nextObservation()
    output_stellar_catalogues = obm.addCatalogue(stellar_cat_file)
    output_galaxy_catalogues = obm.addCatalogue(galaxy_cat_file)
    psf_file = obm.addError()
    fits_file, mosaic_file, params = obm.finalize(mosaic=False)


In this case, the output FITS file will be in the variable `fits_file`, and the output catalogues
(showing the actual count rate and position of the sources observed) will be in the variables
`output_stellar_catalogues` and `output_galaxy_catalogues`.


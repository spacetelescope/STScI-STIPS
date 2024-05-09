Additional STIPS Examples
=========================

Below is a condensed example of STIPS usage, similar to the usage in the :doc:`Basic Tutorial <basic_tutorial>`.

It creates a scene from an existing source catalog ``input_sources.txt``, then observes the
scene with the Roman WFI F129 filter offset by 0.5 degrees in RA and rotated by 27 degrees:

.. code-block:: python

    from stips.scene_module import SceneModule
    from stips.observation_module import ObservationModule


    scm = SceneModule()

    stellar = {
                'n_stars': 100,
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
                'n_gals': 10,
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

    obm.addError()

    fits_file, mosaic_file, params = obm.finalize(mosaic=False)

In this case, the output catalog(s) will show the actual applied count rates.
Whether there is only one output catalog or two depends on the input catalog format.


Fast Extended Sources
----------------------

As of version 2.2, STIPS includes an option to inject extended sources in scenes
using a Sersic-profile approximation. This approximation is ~8 times faster than
the current implementation, but it is also less accurate.

To activate this feature, users must turn on the ``fast_galaxy`` flag. This is how the syntax
looks, starting from the examples listed in the :doc:`STIPS Basic Tutorial <basic_tutorial>`.

.. code-block:: python

    observation_parameters = {
                              'instrument': 'WFI',
                              'filters': ['F129'],
                              'detectors': 1,
                              'distortion': False,
                              'background': 0.15,
                              'fast_galaxy': True,
                              'observations_id': 1,
                              'exptime': 1000,
                              'offsets': [offset]
                              }

.. note::

    We caution however that while this method is a useful approximation, the resulting
    integrated flux measurements can be off by a factor of ~2. Furthermore, the central
    pixel at the core of the galaxy should not be trusted, since this can be off by
    multiple orders of magnitude.


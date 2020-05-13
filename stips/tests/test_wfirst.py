from stips.scene_module import SceneModule
from stips.observation_module import ObservationModule

def test_wfirst_observation():
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

    psf_file = obm.addError()

    fits_file, mosaic_file, params = obm.finalize(mosaic=False)
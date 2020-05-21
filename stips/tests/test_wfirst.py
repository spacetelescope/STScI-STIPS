from stips.scene_module import SceneModule
from stips.observation_module import ObservationModule

import pytest


def create_catalogues():
    star_data = {
                    'n_stars': 100,
                    'age_low': 1.0e12, 'age_high': 1.0e12,
                    'z_low': -2.0, 'z_high': -2.0,
                    'imf': 'salpeter', 'alpha': -2.35,
                    'binary_fraction': 0.1,
                    'distribution': 'invpow', 'clustered': True,
                    'radius': 100.0, 'radius_units': 'pc',
                    'distance_low': 10.0, 'distance_high': 10.0,
                    'offset_ra': 0.0, 'offset_dec': 0.0
                }

    galaxy_data = {
                    'n_gals': 25,
                    'z_low': 0.0, 'z_high': 0.2,
                    'rad_low': 0.01, 'rad_high': 2.0,
                    'sb_v_low': 30.0, 'sb_v_high': 25.0,
                    'distribution': 'uniform', 'clustered': False,
                    'radius': 200.0, 'radius_units': 'arcsec',
                    'offset_ra': 0.0, 'offset_dec': 0.0
                  }

    scm = SceneModule()
    stellar_cat_file = scm.CreatePopulation(star_data)
    galaxy_cat_file = scm.CreateGalaxies(galaxy_data)
    
    return stellar_cat_file, galaxy_cat_file


def get_default_obs():
    obs = {
            'instrument': 'WFI',
            'filters': ['F129'],
            'detectors': 1,
            'distortion': False,
            'oversample': 1,
            'psf_grid_size': 1,
            'pupil_mask': '',
            'background': 'avg',
            'observations_id': 1,
            'exptime': 1000,
            'memmap': True,
            'parallel': True,
            'convolve_size': 4096,
            'offsets': [
                        {
                            'offset_id': 1, 
                            'offset_centre': False, 
                            'offset_ra': 0.5, 
                            'offset_dec': 0.0, 
                            'offset_pa': 27.0
                        }
                       ]
          }
    return obs


def test_wfirst_observation():

    stellar_cat_file, galaxy_cat_file = create_catalogues()
    
    obs = get_default_obs()

    obm = ObservationModule(obs)
    obm.nextObservation()
    output_stellar_catalogues = obm.addCatalogue(stellar_cat_file)
    output_galaxy_catalogues = obm.addCatalogue(galaxy_cat_file)
    psf_file = obm.addError()
    fits_file, mosaic_file, params = obm.finalize(mosaic=False)


@pytest.mark.veryslow
def test_wfirst_observation_deluxe():

    stellar_cat_file, galaxy_cat_file = create_catalogues()

    obs = get_default_obs()
    obs['psf_grid_size'] = 3
    obs['oversample'] = 5

    obm = ObservationModule(obs)
    obm.nextObservation()
    output_stellar_catalogues = obm.addCatalogue(stellar_cat_file)
    output_galaxy_catalogues = obm.addCatalogue(galaxy_cat_file)
    psf_file = obm.addError()
    fits_file, mosaic_file, params = obm.finalize(mosaic=False)


obs_data = [
    (
        {
            'psf_grid_size': 3,
        },
    ),
    (
        {
            'oversample': 3,
        },
    ),
    (
        {
            'memmap': False,
        },
    ),
    (
        {
            'filters': ['F106'],
        },
    ),
    (
        {
            'parallel': False,
        },
    )
]

@pytest.mark.veryslow
@pytest.mark.parametrize(("obs_changes"), obs_data)
def test_obs_parameters(obs_changes):
    
    stellar_cat_file, galaxy_cat_file = create_catalogues()

    obs = get_default_obs()
    for key in obs_changes[0]:
        obs[key] = obs_changes[0][key]

    obm = ObservationModule(obs)
    obm.nextObservation()
    output_stellar_catalogues = obm.addCatalogue(stellar_cat_file)
    output_galaxy_catalogues = obm.addCatalogue(galaxy_cat_file)
    psf_file = obm.addError()
    fits_file, mosaic_file, params = obm.finalize(mosaic=False)

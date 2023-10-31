from stips.scene_module import SceneModule
from stips.observation_module import ObservationModule

import pytest

from tempfile import TemporaryDirectory


def create_catalogues(out_path):
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

    scm = SceneModule(out_path=out_path)
    stellar_cat_file = scm.CreatePopulation(star_data)
    galaxy_cat_file = scm.CreateGalaxies(galaxy_data)

    return stellar_cat_file, galaxy_cat_file


def get_default_obs():
    obs = {
            'instrument': 'WFI',
            'filters': ['F158'],
            'detectors': 1,
            'distortion': False,
            'background': 'avg',
            'observations_id': 1,
            'exptime': 1000,
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


def test_roman_observation():

    dir_name = TemporaryDirectory()

    stellar_cat_file, galaxy_cat_file = create_catalogues(dir_name.name)

    obs = get_default_obs()

    obm = ObservationModule(obs, out_path=dir_name.name, cat_path=dir_name.name)
    obm.nextObservation()
    obm.addCatalogue(stellar_cat_file)
    obm.addCatalogue(galaxy_cat_file)
    obm.addError()
    obm.finalize(mosaic=False)


@pytest.mark.veryslow
def test_roman_observation_deluxe():

    dir_name = TemporaryDirectory()

    stellar_cat_file, galaxy_cat_file = create_catalogues(dir_name.name)

    obs = get_default_obs()

    obm = ObservationModule(obs, out_path=dir_name.name, cat_path=dir_name.name)
    obm.nextObservation()
    obm.addCatalogue(stellar_cat_file)
    obm.addCatalogue(galaxy_cat_file)
    obm.addError()
    obm.finalize(mosaic=False)


obs_data = [
    (
        {
            'filters': ['F106'],
        },
    )
]


@pytest.mark.veryslow
@pytest.mark.parametrize(("obs_changes"), obs_data)
def test_obs_parameters(obs_changes):

    dir_name = TemporaryDirectory()

    stellar_cat_file, galaxy_cat_file = create_catalogues(dir_name.name)

    obs = get_default_obs()
    for key in obs_changes[0]:
        obs[key] = obs_changes[0][key]

    obm = ObservationModule(obs, out_path=dir_name.name, cat_path=dir_name.name)
    obm.nextObservation()
    obm.addCatalogue(stellar_cat_file)
    obm.addCatalogue(galaxy_cat_file)
    obm.addError()
    obm.finalize(mosaic=False)

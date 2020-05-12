from stips.observation_module import ObservationModule
import numpy as np


if __name__ == '__main__':
    filename = 'example_H158.cat'
    seed = np.random.randint(9999)+1000
    scene_general = {'ra': 25.65,'dec': -37.21,'pa': 0.0, 'seed': seed}
    obs = {'instrument': 'WFI', 'filters': ['H158'], 'detectors': 1,'distortion': False, 'oversample': 1,'pupil_mask': '', 'background': 'custom', 'custom_background': 1.5, 'observations_id': 69, 'exptime': 10000.0,'offsets': [{'offset_id': 3, 'offset_centre': False,'offset_ra': 0.0, 'offset_dec': 0.0, 'offset_pa': 0.0}]}
    obm = ObservationModule(obs, scene_general=scene_general)
    obm.nextObservation()
    source_count_catalogues = obm.addCatalogue(str(filename))
    psf_file = obm.addError()
    fits_file, mosaic_file, params = obm.finalize(mosaic=False)

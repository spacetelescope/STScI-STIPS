import os, pytest

from stips import stips_data_base
from stips.utilities import SelectParameter
from stips.utilities.utilities import GetParameter

@pytest.fixture(autouse=True)
def pre_post_test():

    # Setup config file environment variable
    config_param = None
    if "stips_config" in os.environ:
        config_param = os.environ["stips_config"]
        del os.environ["stips_config"]
    
    # Setup stips_data_base by renaming any possible file
    if os.path.exists(os.path.join(stips_data_base, "stips_config.yaml")):
        os.rename(os.path.join(stips_data_base, "stips_config.yaml"),
                  os.path.join(stips_data_base, "stips_config_notused.yaml"))
    
    # this is where the test function runs
    yield
    
    # Teardown config file environment variable
    if config_param is not None:
        os.environ["stips_config"] = config_param
    
    # Teardown stips_data_base config file
    if os.path.exists(os.path.join(stips_data_base, "stips_config_notused.yaml")):
        os.rename(os.path.join(stips_data_base, "stips_config_notused.yaml"),
                  os.path.join(stips_data_base, "stips_config.yaml"))


def test_local_file(data_base):
    config_file = os.path.join(data_base, "override_config.yaml")
    with open(config_file, "w") as conf:
        conf.write("psf_grid_default_size : 17\n")
        conf.write("observation_distortion_enable : true")
    
    # A parameter in the local file
    assert SelectParameter('psf_grid_size', config_file=config_file) == 17
    assert GetParameter('psf_grid_size', config_file=config_file) is None

    assert SelectParameter('psf_grid_default_size', config_file=config_file) == 17
    assert GetParameter('psf_grid_default_size', config_file=config_file) == 17

    # A parameter not in the local file
    assert SelectParameter('observation_detector_oversample', config_file=config_file) == 1
    assert GetParameter('observation_detector_oversample', config_file=config_file) == 1

    assert SelectParameter('oversample', config_file=config_file) == 1
    assert GetParameter('oversample', config_file=config_file) is None
    
    if os.path.exists(config_file):
        os.remove(config_file)


def test_override_dict(data_base):
    config_dict = {
                    "psf_grid_default_size": 17,
                    "observation_distortion_enable": True
                  }
    
    # A parameter in the local file
    assert SelectParameter('psf_grid_size', override_dict=config_dict) == 17

    assert SelectParameter('psf_grid_default_size', override_dict=config_dict) == 17

    # A parameter not in the local file
    assert SelectParameter('observation_detector_oversample', override_dict=config_dict) == 1

    assert SelectParameter('oversample', override_dict=config_dict) == 1


def test_environment_variable(data_base):
    config_file = os.path.join(data_base, "override_config.yaml")
    with open(config_file, "w") as conf:
        conf.write("psf_grid_default_size : 17\n")
        conf.write("observation_distortion_enable : true")
    os.environ['stips_config'] = config_file
    
    # A parameter in the local file
    assert SelectParameter('psf_grid_size') == 17
    assert GetParameter('psf_grid_size') is None

    assert SelectParameter('psf_grid_default_size') == 17
    assert GetParameter('psf_grid_default_size') == 17

    # A parameter not in the local file
    assert SelectParameter('observation_detector_oversample') == 1
    assert GetParameter('observation_detector_oversample') == 1

    assert SelectParameter('oversample') == 1
    assert GetParameter('oversample') is None
    
    if os.path.exists(config_file):
        os.remove(config_file)
    if 'stips_config' in os.environ:
        del os.environ['stips_config']


def test_data_variable(data_base):
    config_file = os.path.join(stips_data_base, "stips_config.yaml")
    with open(config_file, "w") as conf:
        conf.write("psf_grid_default_size : 17\n")
        conf.write("observation_distortion_enable : true")
    
    # A parameter in the local file
    assert SelectParameter('psf_grid_size') == 17
    assert GetParameter('psf_grid_size') is None

    assert SelectParameter('psf_grid_default_size') == 17
    assert GetParameter('psf_grid_default_size') == 17

    # A parameter not in the local file
    assert SelectParameter('observation_detector_oversample') == 1
    assert GetParameter('observation_detector_oversample') == 1

    assert SelectParameter('oversample') == 1
    assert GetParameter('oversample') is None
    
    if os.path.exists(config_file):
        os.remove(config_file)

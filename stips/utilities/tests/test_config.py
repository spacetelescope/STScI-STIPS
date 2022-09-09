import os
import pytest

from stips import stips_data_base
# from stips.utilities import SelectParameter
# from stips.utilities.utilities import GetParameter


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
        conf.write("observation_distortion_enable : true")

    if os.path.exists(config_file):
        os.remove(config_file)


def test_environment_variable(data_base):
    config_file = os.path.join(data_base, "override_config.yaml")
    with open(config_file, "w") as conf:
        conf.write("observation_distortion_enable : true")
    os.environ['stips_config'] = config_file

    if os.path.exists(config_file):
        os.remove(config_file)
    if 'stips_config' in os.environ:
        del os.environ['stips_config']


def test_data_variable(data_base):
    config_file = os.path.join(stips_data_base, "stips_config.yaml")
    with open(config_file, "w") as conf:
        conf.write("observation_distortion_enable : true")

    if os.path.exists(config_file):
        os.remove(config_file)

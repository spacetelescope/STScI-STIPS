"""
General CGI form functions.

:Author: Pey Lian Lim

:Organization: Space Telescope Science Institute

"""
# External modules
import importlib
import inspect
import numpy as np
import os
import requests
import sys
import tarfile
import yaml
from astropy.io import ascii
from astropy.table import Table
from scipy.special import gamma, gammaincinv

from .. import __version__ as __stips__version__


def rind(x):
    """
    Convenience function to take a float, round it to the nearest integer, and convert it
    to the ``int`` type.
    """
    return np.round(x).astype(int)

def get_pandeia_background(wfi_filter, webapp = False):
    """
    Import Pandeia functions to calculate the image background in a given
    filter, in units of electrons per second.

    Parameters
    ----------
    wfi_filter : str
        Name of WFI filter
    webapp: bool
        Toggle strict API checking

    Returns
    -------
    total_back: float
        Total background in units of electrons per second
    """
    from pandeia.engine.calc_utils import build_default_calc
    from pandeia.engine.etc3D import setup, Scene
    from pandeia.engine.observation import Observation
    from pandeia.engine.signal import DetectorSignal
    from pandeia.engine.background import Background
    import scipy.signal as sg

    # Create default configuration file from Pandeia
    calc_input = build_default_calc('roman','wfi','imaging')
    calc_input['configuration']['instrument']['filter'] = wfi_filter.lower()

    # Setup Observation
    calc_config, instrument, strategy, scene_configuration, background, background_level, warnings = setup(calc_input, webapp=webapp)

    # Create Scence
    scene = Scene('roman', input=scene_configuration, webapp=webapp)

    # Create Observation
    observation = Observation(
        scene=scene,
        instrument=instrument,
        strategy=strategy,
        background=background,
        background_level=background_level,
        webapp=webapp
        )

    # Calculate signal from background, without IPC
    observation.instrument.the_detector.ipc = False
    my_detector_signal = DetectorSignal(observation, calc_config=calc_config, webapp=webapp, empty_scene=True)

    # Extract background rate
    rate_plus_bg = my_detector_signal.rate_plus_bg_list[0]
    rate_per_pix = rate_plus_bg['fp_pix']

    # Include IPC effects
    kernel = observation.instrument.get_ipc_kernel()
    detector = sg.fftconvolve(rate_per_pix, kernel, mode='same')

    # Total background
    total_back = np.mean(detector[2:-2,2:-2])

    return total_back

def remove_subfolder(tar_file, subfolder):
    """
    Utility function to take a tar file and remove a set of leading folders from each
    member in the tar file.
    """
    subfolder_len = len(subfolder)
    for member in tar_file.getmembers():
        if member.path.startswith(subfolder):
            member.path = member.path[subfolder_len:]
            yield member


def get_compressed_file(url, file_name, path="", dirs_to_remove=""):
    """
    Utility function to retrieve a .tgz compressed file from a given URL, and extract it
    to a provided path.
    """
    r = requests.get(url, allow_redirects=True)
    with open(file_name, 'wb') as output_file:
        output_file.write(r.content)
    with tarfile.open(file_name) as input_file:
        if len(dirs_to_remove) > 0:
            members = remove_subfolder(input_file, dirs_to_remove)
        else:
            members = input_file.getmembers()
        input_file.extractall(path=path, members=members)
    os.remove(file_name)


class classproperty(object):
    def __init__(self, f):
        self.f = classmethod(f)

    def __get__(self, *a):
        return self.f.__get__(*a)()


class __grid__(object):
    @classproperty
    def __pandeia__version__(self):
        with open(GetStipsData(os.path.join("grid", "VERSION.txt")), "r") as inf:
            line = inf.readline()
            items = line.strip().split()
            return items[1]

    @classproperty
    def __stips__version__(self):
        with open(GetStipsData(os.path.join("grid", "VERSION.txt")), "r") as inf:
            line = inf.readline()
            line = inf.readline()
            items = line.strip().split()
            return items[1]


class StipsEnvironment(object):
    @classproperty
    def __stips__version__(self):
        return __stips__version__

    @classproperty
    def __stips__data__location__(self):
        if 'stips_data' in os.environ:
            return os.environ['stips_data']
        return 'UNSET'

    @classproperty
    def __stips__data__version__(self):
        if 'stips_data' in os.environ:
            fname = os.path.join(os.environ['stips_data'], 'VERSION.txt')
            if os.path.isfile(fname):
                with open(fname, 'r') as inf:
                    line = inf.readline()
                    while len(line.strip()) == 0:
                        line = inf.readline()
                    return line.strip()
            return 'stips_data HAS NO VERSION FILE'
        return 'NO REFERENCE DATA SET'

    @classproperty
    def __stips__grid__version__(self):
        if 'stips_data' in os.environ:
            fname = os.path.join(os.environ['stips_data'], 'grid', 'VERSION.txt')
            if os.path.isfile(fname):
                with open(fname, 'r') as inf:
                    line = inf.readline()
                    line = inf.readline()
                    items = line.strip().split()
                    return items[1]
            return 'stips_data GRID HAS NO VERSION FILE'
        return 'NO REFERENCE DATA SET'

    @classproperty
    def __pandeia__version__(self):
        import pandeia.engine
        if hasattr(pandeia.engine, '__version__'):
            return pandeia.engine.__version__
        return 'UNKNOWN'

    @classproperty
    def __pandeia__data__location__(self):
        if 'pandeia_refdata' in os.environ:
            return os.environ["pandeia_refdata"]
        return 'UNSET'

    @classproperty
    def __pandeia__data__version__(self):
        if 'pandeia_refdata' in os.environ:
            fname = os.path.join(os.environ['pandeia_refdata'], 'VERSION_PSF')
            if os.path.isfile(fname):
                with open(fname, 'r') as inf:
                    line = inf.readline()
                    while len(line.strip()) == 0:
                        line = inf.readline()
                    return line.strip()
            return 'NO VERSION_PSF FILE FOUND'
        return 'NO REFERENCE DATA SET'

    @classproperty
    def __webbpsf__version__(self):
        import webbpsf
        if hasattr(webbpsf, '__version__'):
            return webbpsf.__version__
        return 'UNKNOWN'

    @classproperty
    def __webbpsf__data__location__(self):
        if 'WEBBPSF_PATH' in os.environ:
            return os.environ['WEBBPSF_PATH']
        return 'UNSET'

    @classproperty
    def __webbpsf__data__version__(self):
        if 'WEBBPSF_PATH' in os.environ:
            fname = os.path.join(os.environ['WEBBPSF_PATH'], 'version.txt')
            if os.path.isfile(fname):
                with open(fname, 'r') as inf:
                    line = inf.readline()
                    while len(line.strip()) == 0:
                        line = inf.readline()
                    return line.strip()
            return 'WEBBPSF_PATH HAS NO version.txt FILE'
        return 'NO REFERENCE DATA SET'

    @classproperty
    def __stips__environment__dict__(self):
        env = StipsEnvironment
        import astropy
        import photutils
        env_dict = {
                    'stips_version': env.__stips__version__,
                    'stips_data_location': env.__stips__data__location__,
                    'stips_data_version': env.__stips__data__version__,
                    'stips_grid_version': env.__stips__grid__version__,
                    'pandeia_version': env.__pandeia__version__,
                    'pandeia_data_location': env.__pandeia__data__location__,
                    'pandeia_data_version': env.__pandeia__data__version__,
                    'webbpsf_version': env.__webbpsf__version__,
                    'webbpsf_data_location': env.__webbpsf__data__location__,
                    'webbpsf_data_version': env.__webbpsf__data__version__,
                    'astropy_version': astropy.__version__,
                    'photutils_version': photutils.__version__
                   }
        return env_dict

    @classproperty
    def __stips__environment__report__(self):
        env = StipsEnvironment.__stips__environment__report__pretty__
        return env.replace("\n", " ").replace("\t", " ")

    @classproperty
    def __stips__environment__report__pretty__(self):
        env = StipsEnvironment.__stips__environment__dict__
        report = ""
        report += "STIPS Version {} with Data Version {} at {}.\n".format(env['stips_version'], env['stips_data_version'], env['stips_data_location'])
        report += "\tSTIPS Grid Generated with {}\n".format(env['stips_grid_version'])
        report += "Pandeia Version {} with Data Version {} at {}.\n".format(env['pandeia_version'], env['pandeia_data_version'], env['pandeia_data_location'])
        report += "Webbpsf Version {} with Data Version {} at {}.\n".format(env['webbpsf_version'], env['webbpsf_data_version'], env['webbpsf_data_location'])
        return report


def SetupDataPaths():
    """
    Set up the STIPS, synphot, webbpsf, and pandeia reference data environment variables.
    """
    for item in ["stips", "synphot", "webbpsf", "pandeia"]:
        var_name = GetParameter(item+"_data_name", use_data=False)
        if var_name not in os.environ:
            var_path = GetParameter(item+"_data", use_data=False)
#             print("Setting up {} to {}".format(var_name, var_path))
            if var_path == "$local":
                if item == "stips":
                    file_dir = os.path.dirname(os.path.abspath(__file__))
                    data_dir = os.path.join(file_dir, "..", "..", "ref_data", "stips_data")
                else:
                    data_dir = os.path.join(os.environ["stips_data"], "ref", var_name)
                var_path = os.path.normpath(data_dir)
            os.environ[var_name] = var_path
#             print("Set {} to {}".format(var_name, var_path))


def DownloadReferenceData():
    """
    Set up the STIPS, synphot, webbpsf, and pandeia reference data environment variables.
    """
    SetupDataPaths()

    # STIPS
    print("Checking STIPS data")
    stips_data_file = "stips_data-1.0.9.tgz"
    stips_url = "https://stsci.box.com/shared/static/4nebx2ndxr7c77lgocfbvxo7c2hyd3in.tgz"
    stips_data_path = os.environ[GetParameter("stips_data_name", use_data=False)]
    if not os.path.isdir(stips_data_path):
        print("Downloading STIPS data to {}".format(stips_data_path))
        os.makedirs(stips_data_path)
        get_compressed_file(stips_url, stips_data_file, stips_data_path, "stips_data/")
    else:
        print("Found at {}".format(stips_data_path))

    # synphot/stsynphot
    print("Checking synphot data")
    synphot_url = "https://ssb.stsci.edu/trds/tarfiles"
    synphot_data_path = os.environ[GetParameter("synphot_data_name", use_data=False)]
    if not os.path.isdir(synphot_data_path):
        print("Downloading synphot data to {}".format(synphot_data_path))
        os.makedirs(synphot_data_path)
        for i in range(1, 8):
            file_name = "synphot{}.tar.gz".format(i)
            print("\tDownloading {}".format(file_name))
            url = synphot_url+"/"+file_name
            get_compressed_file(url, file_name, synphot_data_path, "grp/redcat/trds/")
    else:
        print("Found at {}".format(synphot_data_path))

    # webbpsf
    print("Checking webbpsf data")
    webbpsf_url = "https://stsci.box.com/shared/static/t90gqazqs82d8nh25249oq1obbjfstq8.gz"
    webbpsf_data_path = os.environ[GetParameter("webbpsf_data_name", use_data=False)]
    webbpsf_data_file = "webbpsf_data.tar.gz"
    if not os.path.isdir(webbpsf_data_path):
        print("Downloading webbpsf data to {}".format(webbpsf_data_path))
        os.makedirs(webbpsf_data_path)
        get_compressed_file(webbpsf_url, webbpsf_data_file, webbpsf_data_path,
                            "webbpsf-data/")
    else:
        print("Found at {}".format(webbpsf_data_path))

    # pandeia
    print("Checking pandeia data")
    pandeia_data_file = "pandeia_data-3.1_roman.tar.gz"
    pandeia_url = "https://stsci.box.com/shared/static/cmljh0lsffz4345064eso7lix70f9477.gz"
    pandeia_data_path = os.environ[GetParameter("pandeia_data_name", use_data=False)]
    if not os.path.isdir(pandeia_data_path):
        print("Downloading pandeia data to {}".format(pandeia_data_path))
        os.makedirs(pandeia_data_path)
        get_compressed_file(pandeia_url, pandeia_data_file, pandeia_data_path,
                            "pandeia_data-3.1_roman/")
    else:
        print("Found at {}".format(pandeia_data_path))

    # Done.


def GetStipsDataDir():
    """
    Get the STIPS data directory path.
    """
    stips_version = __stips__version__.strip().replace(".", "")
    if "dev" in stips_version:
        stips_version = stips_version[:stips_version.find("dev")]
    if "stips_data" not in os.environ:
        msg = "ERROR: stips_data environment variable not found. STIPS "
        msg += "requires the data directory to function. Please download the "
        msg += "STIPS data from <https://stsci.box.com/v/stips-data-{}> "
        msg += "and set the stips_data environment variable to point to it.\n"
        sys.stderr.write(msg.format(stips_version))
        raise EnvironmentError("stips_data environment variable not found.")
    stips_data_base = os.environ["stips_data"]
    if not os.path.exists(stips_data_base):
        msg = "ERROR: stips_data directory at {} not ".format(stips_data_base)
        msg += "found. STIPS requires the data directory to function. "
        msg += "Please make sure that the STIPS data directory exists.\n"
        sys.stderr.write(msg)
        raise FileNotFoundError("${stips_data} does not exist.")
    return stips_data_base


def GetStipsData(to_retrieve):
    """
    Retrieve a file from the stips_data directory. Will also print out a warning if the directory
    can't be found.
    """
    stips_data_base = GetStipsDataDir()
    stips_version = __stips__version__.strip().replace(".", "")
    if "dev" in stips_version:
        stips_version = stips_version[:stips_version.find("dev")]
    retrieval_file = os.path.join(stips_data_base, to_retrieve)
    if not os.path.exists(retrieval_file):
        msg = "ERROR: STIPS data file {} not found. ".format(retrieval_file)
        msg += "STIPS requires the data directory to function. Please download "
        msg += "a new copy of the STIPS data from "
        msg += "<https://stsci.box.com/v/stips-data-{}>.\n"
        sys.stderr.write(msg.format(stips_version))
        raise FileNotFoundError("File {} does not exist.".format(retrieval_file))
    return retrieval_file


def SelectParameter(name, override_dict=None, config_file=None):
    """
    If override_dict contains the key name, return override_dict[name].
    Otherwise, if the parameter name is present in the configuration file,
    return the value found in the configuration file. Otherwise, if an alternate
    name (as defined in a local dictionary) is found in the configuration file,
    return the value for that name. Otherwise, return None.

    Parameters
    ----------
    name : str
        Name of parameter

    override_dict : dict, default None
        Dictionary that may override a configuration value

    config_file : str, default None
        Supplied configuration file

    Returns
    -------
    value : obj
        The value found (None if no value is found)
    """
    name_mappings = {
                        'background': 'observation_default_background',
                        'background_location': 'observation_jbt_location',
                        'cat_path': 'input_location',
                        'cat_type': 'catalogue_type',
                        'convolve': 'residual_convolve_psf',
                        'convolve_size': 'psf_convolution_max_size',
                        'cores': 'parallel_ncores',
                        'distortion': 'observation_distortion_enable',
                        'jbt_location': 'observation_jbt_location',
                        'out_path': 'output_location',
                        'psf_grid_size': 'psf_grid_default_size',
                        'seed': 'random_seed'
                    }

    if override_dict is not None:
        if name in override_dict:
            return override_dict[name]
        elif name in name_mappings and name_mappings[name] in override_dict:
            return override_dict[name_mappings[name]]

    value = GetParameter(name, config_file)
    if value is not None:
        return value
    elif name in name_mappings:
        return GetParameter(name_mappings[name], config_file)

    return None


def GetParameter(param, config_file=None, use_provided=True, use_cwd=True,
                 use_environ=True, use_data=True):
    """
    Retrieve a parameter from the STIPS configuration file. This function looks
    for the STIPS configuration file as follows (returning the first file found)

    - If a file is provided to the function, check that file
    - If there is a stips_config.yaml in the current directory, check that file
    - If there is a stips_config environment variable, check that file
    - If there is a "stips_config.yaml" file in stips_data, check that file
    - Check the internal data/stips_config.yaml file.

    Parameters
    ----------
    param : str
        Name of parameter

    config_file : str, default None
        Supplied configuration file

    Returns
    -------
    value : obj
        The value found (None if no value is found)
    """
    file_used = "local"
    settings = None
    conf_file = None
    if use_data:
        stips_data_dir = GetStipsDataDir()
    local_dir = os.path.dirname(os.path.abspath(__file__))
    local_data_dir = os.path.join(local_dir, "..", "data")
    local_config_file = os.path.join(local_data_dir, "stips_config.yaml")
    local_config = os.path.normpath(local_config_file)
    cwd_config = os.path.join(os.getcwd(), "stips_config.yaml")
    if use_data:
        data_config = os.path.join(stips_data_dir, "stips_config.yaml")
    env_config = os.environ.get('stips_config', None)
    if use_environ and env_config is not None:
        if not os.path.isfile(env_config):
            env_config = os.path.join(env_config, "stips_config.yaml")

    if use_provided and config_file is not None and os.path.isfile(config_file):
        conf_file = config_file
        file_used = "provided"
    elif use_cwd and os.path.isfile(cwd_config):
        conf_file = cwd_config
        file_used = "cwd"
    elif use_environ and 'stips_config' in os.environ and os.path.isfile(env_config):
        conf_file = env_config
        file_used = "environ"
    elif use_data and os.path.isfile(data_config):
        file_used = "data"
        conf_file = data_config
    elif os.path.isfile(local_config):
        conf_file = local_config

    if conf_file is not None:
        with open(conf_file, 'r') as config:
            settings = yaml.safe_load(config)

    if settings is not None and param in settings:
        return TranslateParameter(param, settings[param], data_dir=use_data)
    elif param not in settings and file_used == "provided":
        # Try without the supplied config file in case it doesn't include
        #   the full set of parameters
        return GetParameter(param, use_provided=False)
    elif param not in settings and file_used == "cwd":
        # Try turning off the cwd flag as well
        return GetParameter(param, use_provided=False, use_cwd=False)
    elif param not in settings and file_used == "environ":
        # Try without the environment variable config in case it doesn't include
        #   the full set of parameters
        return GetParameter(param, use_provided=False, use_cwd=False,
                            use_environ=False)
    elif param not in settings and file_used == "data":
        # Try without the stips_data config in case it doesn't include
        #   the full set of parameters
        return GetParameter(param, use_provided=False, use_cwd=False,
                            use_environ=False, use_data=False)

    return None


def TranslateParameter(param, value, data_dir=True):
    """
    Check if a parameter is in a dictionary of special values and, if so,
    substitute in the proper value.

    Parameters
    ----------
    param : str
        Name of parameter

    value : obj
        Supplied value

    Returns
    -------
    value : obj
        The value as translated
    """
    translations = {
                    'input_location': {'$CWD': os.getcwd()},
                    'output_location': {'$CWD': os.getcwd()}
                   }
    # It is possible that this function will be called before the STIPS data directory
    # has been set up.
    if data_dir:
        translations['psf_cache_location'] = {'$DATA': GetStipsDataDir()}

    if param in translations:
        if value in translations[param]:
            if callable(translations[param][value]):
                return translations[param][value]()
            return translations[param][value]
        return value

    return value


def OffsetPosition(in_ra, in_dec, delta_ra, delta_dec):
    """
    Offset a position given in decimal degrees.

    Parameters
    ----------
    in_ra: float
        Initial RA (decimal degrees).

    in_dec: float
        Initial DEC (demical degrees).

    delta_ra: float
        Offset in RA (decimal degrees).

    delta_dec: float
        Offset in DEC (decimal degrees).

    Returns
    -------
    ra: float
        Offset RA.

    dec: float
        Offset DEC.
    """
    ra = in_ra
    dec = in_dec + delta_dec
    if dec > 90.:
        dec = 180 - dec
        ra = 180 + ra
    if dec < -90.:
        dec = -180 - dec
        ra = 180 + ra
    ra = ra + delta_ra
    if ra > 360.:
        ra = ra - 360.
    if ra < 0.:
        ra = ra + 360.
    return ra, dec


def InstrumentList(excludes=[]):
    """
    Looks through the available instrument files, and imports all instruments that do not match
    one of the names in 'excludes', or whose telescope does not match one of the names in 'excludes'

    Parameters
    ----------

    excludes: list of strings
        Class and Telescope names to exclude when building the array.

    Returns
    -------
    instruments: dictionary
        dictionary of instrument classes (by instrument name)
    """
    instruments = {}
    from ..instruments import __all__ as all_instruments
    for f in all_instruments:
        if f == "__init__":
            continue
        module = importlib.import_module(".instruments."+f, "stips")
        for name, obj in inspect.getmembers(module):
            if inspect.isclass(obj) and hasattr(obj, "TELESCOPE") and hasattr(obj, "INSTRUMENT") and hasattr(obj, "DETECTOR"):
                telescope = getattr(obj, "TELESCOPE").lower()
                instrument = getattr(obj, "INSTRUMENT").lower()
                detector = getattr(obj, "DETECTOR").lower()
                add = True
                for exclude in excludes:
                    excl = exclude.lower()
                    if excl == telescope or excl == instrument or excl == detector or excl == name.lower():
                        add = False
                if add:
                    instruments[name] = obj
    return instruments


def read_metadata(filename, n_lines=100000, format='ascii.ipac'):
    lines = []
    for i, line in enumerate(open(filename, 'r')):
        lines.append(line.rstrip(os.linesep))
        if i == n_lines:
            break
    t = Table.read(lines, format=format)
    return t


def read_table(filename, n_chunk=100000, format="ipac"):
    """
    Chunk reader to (hopefully) not take up so much memory when reading very large tables.
    """
    names = None
    lines = []
    for i, line in enumerate(open(filename, 'r')):
        lines.append(line.rstrip(os.linesep))
        if i % n_chunk == n_chunk - 1:
            if i < n_chunk:  # first
                chunk = ascii.read(lines, format=format, guess=False)
                names = chunk.colnames
                yield chunk
            else:
                chunk = ascii.read(lines, format='no_header', names=names, guess=False)
                yield chunk
            lines = []
    if lines:
        if names is not None:
            yield ascii.read(lines, format='no_header', names=names, guess=False)
        else:
            yield ascii.read(lines, format=format, guess=False)


def b(n):
    # Normalisation constant
    return gammaincinv(2*n, 0.5)


def sersic_lum(Ie, re, n):
    # total luminosity (integrated to infinity)
    bn = b(n)
    g2n = gamma(2*n)
    return Ie * re**2 * 2*np.pi*n * np.exp(bn)/(bn**(2*n)) * g2n


if __name__ == '__main__':
    pass

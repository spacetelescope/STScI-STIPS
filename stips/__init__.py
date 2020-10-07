from __future__ import absolute_import

import os

# make use of astropy affiliate framework to set __version__, __githash__, and
# add the test() helper function

from ._astropy_init import *

__all__ = ['astro_image', 'errors', 'galaxy_module', 'instruments',
           'observation_module', 'scene_module', 'stellar_module',
           'utilities', 'version', 'test']

try:
    from .version import version as __version__
except ImportError:
    __version__ = 'dev'
version = __version__

# Local Definitions

from .astro_image import AstroImage
from .instruments import Instrument
from .observation_module import ObservationModule
from .scene_module import SceneModule
from .utilities import __grid__
from .utilities import CachedJbtBackground
from .utilities import GetStipsData
from .utilities import internet
from .utilities import StipsEnvironment

__grid__pandeia__version__ = StipsEnvironment.__pandeia__version__
__grid__stips__version__ = StipsEnvironment.__stips__grid__version__
__env__report__ = StipsEnvironment.__stips__environment__report__pretty__
__env__dict__ = StipsEnvironment.__stips__environment__dict__

try:
     stips_data_base = os.environ["stips_data"]
except KeyError:
     stips_data_base = None

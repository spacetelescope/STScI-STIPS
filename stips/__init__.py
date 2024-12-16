from __future__ import absolute_import

import os

__all__ = ['observation_module', 'scene_module']

__version__ = "2.2.2"
version = __version__

from .utilities import SetupDataPaths

SetupDataPaths()

# Local Definitions
from .observation_module import ObservationModule
from .scene_module import SceneModule
from .utilities import GetStipsData
from .utilities import StipsEnvironment
from .utilities import DownloadReferenceData

__grid__pandeia__version__ = StipsEnvironment.__pandeia__version__
__grid__stips__version__ = StipsEnvironment.__stips__grid__version__
__env__report__ = StipsEnvironment.__stips__environment__report__pretty__
__env__dict__ = StipsEnvironment.__stips__environment__dict__

try:
    stips_data_base = os.environ["stips_data"]
except KeyError:
    stips_data_base = None

#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import builtins
import glob

from setuptools import setup
from setuptools.config import read_configuration

# Store the package name in a built-in variable so it's easy
# to get from other parts of the setup infrastructure
PACKAGENAME = read_configuration('setup.cfg')['metadata']['name']

# Freeze build information in version.py. Note that this gets information
# about the package (name and version) from the setup.cfg file.
version = read_configuration('setup.cfg')['metadata']['version']

setup(
        version=version, 
        package_data={"stips": ['data/*']}#,
#        scripts=glob.glob("stips/commands/*")
)

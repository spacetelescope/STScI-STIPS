from __future__ import absolute_import

import os
from setuptools import setup, find_packages

with open(os.path.join(os.getcwd(), "stips", "version", "__init__.py")) as f:
    exec(f.readline().strip())

setup(name='stips',
      version=__version__,
      description='Space Telescope Imaging Product Simulator',
      classifiers=[
                   'Development Status :: 4 - Beta',
                   'License :: OSI Approved :: GPLv3 License',
                   'Programming Language :: Python :: 2.7.3',
                   'Topic :: Software Development :: Programmatic Representation',
      ],
      url='http://www.stsci.edu',
      author='Brian York',
      author_email='york@stsci.edu',
      license='GPLv3',
      packages=find_packages(),
      zip_safe=False,
      install_requires=['astropy', 'automodinit', 'esutil', 'numpy', 'photutils', 'pysynphot', 
                        'scipy', 'jwst_backgrounds', 'webbpsf', 'pandeia.engine'],
      include_package_data=True)

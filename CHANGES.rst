#############
Release Notes
#############

Version History and Change Log
------------------------------

Version 2.2.2
=============
- Updated pandeia.engine version to be >= 3.1 instead of == 3.1

Version 2.2.1
=============
- Fixed a bug on the version of STIPS in __init__.py
- Updated WebbPSF version to be >= 1.1.1 instead of == 1.1.1
- Updated synphot>=1.1.1 and stsynphot>=1.1.0 from forced to specific version.

Version 2.2.0
=============
- Added a functionality to allow for faster simulations of extended sources.
- Fixed a bug related to the simulation of extended sources.
- Fixed a bug related to the WCS.
- Updated pandeia to version 3.1 and its corresponding reference files.
- Changed the minimum required Python version to 3.10

Version 2.1.0
=============
- Updated the Flux PHOTFNU values to match Pandeia
- Updated the background and noise estimates to use Pandeia
- Updated webbpsf to version 1.1
- Updated pandeia to version 2.0

Version 2.0.0
=============

- STIPS now uses an ePSF, and calculates source appearance on a per-source basis
- Dropped support for python 2, and python <= 3.7
- Dropped support for HST and JWST instruments
- Dropped support for jbt background tool, as its results are JWST-specific
- Updated pandeia to 1.7
- Updated webbpsf to 1.0.0
- Dropped astropy helpers code


Version 1.0.8
=============
*2020 May 22*

**STIPS Improvements**

- PSFs generated with `webbpsf` are now PSF grids. [:pr:``, :user:`york-stsci`]
- STIPS has the option to keep all data in memory. [:pr:``, :user:`york-stsci`]
- Data files have been removed from the repo to allow STIPS to be uploaded to PyPi (pip). The data has been migrated to a STScI box folder. Users can now download the data and set their `stips_data` to allow access to the data that once lived in the repository. [:pr:`59`, :user:`york-stsci`]
- Travis was setup to run and pass tests. Outdated tests have been removed to allow Travis to pass. [:pr:`65`, :user:`robelgeda`]
- Test data has been moved to a STScI box folder. [:pr:`64`, :user:`robelgeda`]
- `dev` folder added for any developer related tools. [:pr:`69`, :user:`robelgeda`]
- Frozen environments saved at `dev/conda_envs`. [:pr:`69`, :user:`robelgeda`]

Version 1.0.7
=============
*2020 January 8*

**STIPS Improvements**

- Cookie cutter template used to create better package infrastructure. [:pr:`40`, :user:`robelgeda`]
- Docker file added for ease of install. [:pr:`48`, :user:`robelgeda`]
- environment.yml added for easy conda env build. [:pr:`42`, :user:`robelgeda`]
- Read the Docs documentation established. [:pr:`55`, :user:`robelgeda`]
- F062 filter added to Roman WFI. [:pr:`51`, :user:`york-stsci`]
- Update STIPS to use WbbPSF 0.9.0. [:pr:`51`, :user:`york-stsci`]
- Travis CI initiated for unit and regression testing. [:pr:`40`, :user:`robelgeda`]
- Python version set to 3.7 [:pr:`40`, :user:`robelgeda`]
- Licenses updated [:pr:`40`, :user:`robelgeda`]

**General bug fixes and small changes**

- Updated astro_image.py to use a PC matrix rather than a CD matrix for the image WCS, which hopefully will result in astropy actually giving you a correctly formatted FITS WCS. [:pr:`46`, :user:`york-stsci`]
- Adding WCS information to PSF files. PSF files will now have the following:
    - RA equal to the observation RA at which they were produced
    - DEC equal to the observation DEC at which they were produced
    - PA equal to the observation PA at which they were produced
    - CDELT keywords equal to the PIXELSCL keyword, but adjusted to degrees rather than arcsec.
    - [:pr:`47`, :user:`york-stsci`]

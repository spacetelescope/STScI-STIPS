Documentation
=============

.. note::

  As of version 2.0, STIPS no longer offers support for HST or JWST.

Overview
--------
STIPS is the Space Telescope Imaging Product Simulator. It is designed to create
simulations of full-detector post-pipeline astronomical scenes for the Nancy Grace Roman
Space Telescope's Wide-Field Instrument (WFI). STIPS has the ability to add
instrumental distortion (if available) as well as calibration residuals from flatfields,
dark currents, and cosmic rays. It automatically includes Poisson noise and readout noise.
It does not include instrument saturation effects.

Why use STIPS?
--------------
STIPS is intended to produce quick simulations of Level 2 (L2) images, and is provided for
cases where `Pandeia <https://pypi.org/project/pandeia.engine/>`_ does not
provide a large enough simulation area (e.g., full-detector or multiple-detector
observations). STIPS obtains its Roman instrument and filter values from
Pandeia, so it should produce output within 10% of output produced by Pandeia.

STIPS does not start with Level 1 (L1) images and propagate instrumental calibrations
through the simulations. While it does have the ability to add error residuals (representing
the remaining uncertainty after pipeline calibration), these residuals are not validated
against actual pipeline calibrations of L1 images. STIPS is not the ideal choice if
extremely good instrumental accuracy is needed. Pandeia is the preferred tool for
high-accuracy observations.

Developed by Brian York (`@york-stsci <https://github.com/york-stsci>`_),
Robel Geda (`@robelgeda <https://github.com/robelgeda>`_),
and O. Justin Otor (`@ojustino <https://github.com/ojustino>`_).
Python ePSF code developed by
Sebastian Gomez (`@gmzsebastian <https://github.com/gmzsebastian>`_) based on Fortran code
developed by Andrea Bellini (`@AndreaBellini <https://github.com/AndreaBellini>`_).

.. figure:: roman_figures/stips_demo.png

  Fig. 1: Simulated WFI image of a star cluster and background galaxies.

Using STIPS
-----------

.. toctree::
  :maxdepth: 2

  installation
  using_stips
  basic_tutorial
  examples
  bugs
  release
  help

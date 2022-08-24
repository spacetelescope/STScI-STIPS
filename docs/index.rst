Documentation
=============

.. note::

  STIPS 2.0 no longer supports HST or JWST. STIPS 1.0.8 is the most recent version to 
  offer support for those telescopes and instruments.

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
STIPS is intended for cases where an exposure time calculator (e.g. Pandeia) does not
provide enough detector area (e.g. testing photometry code, quick looks at dither patterns
or multi-detector observations of a scene). For Roman, it obtains its background count
levels and instrumental throughput levels from Pandeia internally, so it should produce
output within 10% of output produced by Pandeia.

STIPS is not the ideal choice if extremely good instrumental accuracy is needed. Pandeia
is the preferred tool for high-accuracy observations.

Developed by Brian York (`@york-stsci <https://github.com/york-stsci>`_),
Robel Geda (`@robelgeda <https://github.com/robelgeda>`_), 
and O. Justin Otor (`@ojustino <https://github.com/ojustino>`_). 
Python ePSF code developed by 
Sebastian Gomez (`@gmzsebastian <https://github.com/gmzsebastian>`_).

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
  help

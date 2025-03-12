# STScI-STIPS

[![Build Status](https://github.com/spacetelescope/STScI-STIPS/actions/workflows/ci_workflows.yml/badge.svg?branch=main)](https://github.com/spacetelescope/STScI-STIPS/actions/workflows/ci_workflows.yml?query=branch:main)
[![STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![DOI](https://zenodo.org/badge/94350619.svg)](https://zenodo.org/badge/latestdoi/94350619)

For documentation and installation instructions please visit https://stips.readthedocs.io

## Table of Contents

* [Overview](#overview)
* [Why use STIPS?](#why-use-stips)

## Overview

STIPS is the Space Telescope Imaging Product Simulator. It is designed to create
simulations of full-detector post-pipeline astronomical scenes for the Nancy Grace Roman
Space Telescope's Wide-Field Instrument (WFI). STIPS has the ability to add
instrumental distortion (if available) as well as calibration residuals from flatfields,
dark currents, and cosmic rays. It automatically includes Poisson noise and readout noise.
It does not include instrument saturation effects.

## Why use STIPS?

STIPS is intended to produce quick simulations of Level 2 (L2) images, and is provided for
cases where [Pandeia](https://pypi.org/project/pandeia.engine) does not
provide a large enough simulation area (e.g., full-detector or multiple-detector
observations). STIPS obtains its Roman instrument and filter values from
Pandeia, so it should produce output within 10% of output produced by Pandeia.

STIPS does not start with Level 1 (L1) images and propagate instrumental calibrations
through the simulations. While it does have the ability to add error residuals (representing
the remaining uncertainty after pipeline calibration), these residuals are not validated
against actual pipeline calibrations of L1 images. STIPS is not the ideal choice if
extremely good instrumental accuracy is needed. Pandeia is the preferred tool for
high-accuracy observations.

Developed by Brian York ([@york-stsci](https://github.com/york-stsci)),
Robel Geda ([@robelgeda](https://github.com/robelgeda)), and
O. Justin Otor ([@ojustino](https://github.com/ojustino)).
Python ePSF code developed by
Sebastian Gomez ([@gmzsebastian](https://github.com/gmzsebastian)) based on Fortran code
developed by Andrea Bellini ([@AndreaBellini](https://github.com/AndreaBellini)).

![Alt text](docs/roman_figures/stips_demo.png?raw=true "Roman WFI Image of a Star Cluster and Background Galaxies")

## Citation

If you use STIPS, please cite the [STIPS PASP Paper](https://ui.adsabs.harvard.edu/abs/2024PASP..136l4502S/abstract) that describes the code.

```
@ARTICLE{2024PASP..136l4502S,
       author = {{Stips Development Team} and {Gomez}, Sebastian and {Bellini}, Andrea and {Al-Kowsi}, Hanna and {Desjardins}, Tyler and {Geda}, Robel and {Han}, Eunkyu and {Otor}, O. Justin and {Riedel}, Adric and {Ryan}, Russell and {Spitzer}, Isaac and {York}, Brian},
        title = "{STIPS: The Nancy Grace Roman Space Telescope Imaging Product Simulator}",
      journal = {\pasp},
     keywords = {Astronomical techniques, Telescopes, Astronomical methods, 1684, 1689, 1043, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2024,
        month = dec,
       volume = {136},
       number = {12},
          eid = {124502},
        pages = {124502},
          doi = {10.1088/1538-3873/ad9524},
archivePrefix = {arXiv},
       eprint = {2411.11978},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024PASP..136l4502S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
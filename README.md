# STScI-STIPS

[![Build Status](https://travis-ci.com/spacetelescope/STScI-STIPS.svg?branch=master)](https://travis-ci.com/spacetelescope/STScI-STIPS)
[![STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)

For documentation and installation instructions please visit https://stsci-stips.readthedocs.io

## Table of Contents

* [Overview](#overview)
* [Why use STIPS?](#why-use-stips)

## Overview

STIPS is the Space Telescope Imaging Product Simulator. It is designed to create simulations of 
full-detector post-pipeline astronomical scenes for any telescope. Currently STIPS has modules for
WFC3 IR (F110W and F160W only), JWST (NIRCam Short, NIRCam Long, and MIRI), and Roman (WFI). STIPS
has the ability to add instrumental distortion (if available) as well as calibration residuals
(currently flatfield residuals, dark current residuals, and cosmic ray residuals). It automatically
includes Poisson noise and readout noise. It does not include instrument saturation effects. In
addition, STIPS has the ability to generate its own scenes, consisting of stellar populations and
background galaxies (implemented as Sersic profiles).

## Why use STIPS?

STIPS is intended for cases where an ETC (e.g. Pandeia) does not provide enough detector area (e.g.
testing photometry code, quick looks at dither patterns or multi-detector observations of a scene).
For JWST and Roman, it obtains its background count levels and instrumental throughput levels from
Pandeia internally, so it should produce output within 10% of output produced by Pandeia.

If extremely good instrumental accuracy is needed, STIPS is not the ideal choice. Instead, the
various instrument design teams have produced much more detailed simulators. STIPS is intended to
run reasonably quickly, and to make scene generation and observation as easy as possible.

Developed by Brian York ([@york-stsci](https://github.com/york-stsci)) and
Robel Geda ([@robelgeda](https://github.com/robelgeda)).


![Alt text](docs/roman_figures/stips_demo.png?raw=true "Roman WFI Image of a Star Cluster and Background Galaxies")


The STIPS PSF Grid
==================

Overview
--------

As of version 2.0.0, STIPS uses WebbPSF PSF grids to implement an ePSF. The ePSF
uses a 3x3 grid of PSFs to handle PSF variation across the detector, and an oversample of
4x4 to handle sub-pixel positioning. Each point source has its PSF calculated individually
as it is added to the detector.

Internal Image Size
-------------------

The size of the STIPS pre-convolution image is set by the size of the detector. In
addition, the image is padded by half of the PSF size on each side in order to include
sources that fall off the detector, but whose PSF may land on the detector.

PSF Image Size
--------------

The size of the STIPS generated PSF image depends on the source brightness. For most
sources, a 45x45-pixel PSF is created via an ePSF oversampled to a factor of 4. For bright
sources (magnitude < the bright limit, default 14), a 91x91-pixel PSF is created to
provide flux from the extended wings. For extra-bright sources
(magnitude < the extra-bright limit, default 3), a 181x181-pixel PSF is created to provide
even more wing flux to the model.

The STIPS PSF Grid
==================
.. note::

    The STIPS PSF Grid creation and configuration process is still in flux. At the moment, there is a somewhat confusing conflation of PSF convolution size (which should be mostly to do with available memory and computation speed) and how many individual PSF interpolations are used from a given grid. In the future, these considerations will be split apart, so the documentation here is subject to change.

Overview
--------

Starting with STIPS 1.0.8, STIPS supports (and requires) the use of WebbPSF PSF grids when creating PSFs for convolution. When the STIPS image convolution step is performed, it will divide up the image into squares (with size set by the PSF convolution size configuration item), and a PSF will be generated for each square, with the PSF detector location set to the centre of that square. As a result, different regions of the detector will have different PSFs applied to them.

Internal Image Size
-------------------

The size of the STIPS pre-convolution image is set by the size of the detector. In addition, the image is padded by half of the PSF size on each side in order to include sources that fall off the detector, but whose PSF may land on the detector.

PSF Image Size
--------------

The size of the STIPS generated PSF image depends on the source brightness. For most sources, a 45x45-pixel PSF is created via an ePSF oversampled to a factor of 4. For bright sources (magnitude < the bright limit, default 14), a 91x91-pixel PSF is created to provide flux from the extended wings. For extra-bright sources (magnitude < the extra-bright limit, default 3), a 181x181-pixel PSF is created to provide even more wing flux to the model.

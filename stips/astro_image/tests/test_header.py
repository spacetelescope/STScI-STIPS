import os

import numpy as np

from astropy.wcs import WCS

from stips.astro_image import AstroImage

def new_empty_astro_image(data_base):
    fits_path = os.path.join(data_base, "astro_image", "zero_image.fits")
    print(fits_path)
    return AstroImage.initFromFits(fits_path, psf=False)


def test_add_rdnoise(data_base):
    # Test add one point
    ai = new_empty_astro_image(data_base)  # all data values = 0

    ai.introduceReadnoise(17.5)
    assert ai.hdu.header['rdnoise'] == 17.5

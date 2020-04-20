import os

import numpy as np

from astropy.wcs import WCS

from stips.astro_image import AstroImage
from stips import stips_data_base

def new_empty_astro_image():
    fits_path = os.path.join(stips_data_base, "test", "astro_image", "zero_image.fits")
    return AstroImage.initFromFits(fits_path)


def test_add_points():
    # Test add one point
    ai = new_empty_astro_image()  # all data values = 0

    ai.addPoints([49], [49], [1000])
    assert(ai.hdu.data[49, 49] == 1000)
    assert(ai.sum == 1000)

    # do this one more time to check if the arithmetic is good
    ai.addPoints([49], [49], [1000])
    assert(ai.hdu.data[49, 49] == 2000)
    assert(ai.sum == 2000)

    # Test add a list of points
    ai = new_empty_astro_image()  # all data values = 0
    x = np.arange(0, 10)
    y = np.arange(0, 10)
    value = x*100 + y
    ai.addPoints(x, y, value)

    for i in range(0, 10):
        xi = x[i]
        yi = y[i]
        vi = value[i]
        assert(ai.hdu.data[yi, xi] == vi)
    assert(ai.sum == value.sum())

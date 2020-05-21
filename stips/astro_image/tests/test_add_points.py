import os

import numpy as np

from astropy.wcs import WCS

from stips.astro_image import AstroImage

def new_empty_astro_image(data_base):
    fits_path = os.path.join(data_base, "astro_image", "zero_image.fits")
    print(fits_path)
    return AstroImage.initFromFits(fits_path)


def test_add_points(data_base):
    # Test add one point
    ai = new_empty_astro_image(data_base)  # all data values = 0

    ai.addPoints([49], [49], [1000])
    assert(ai.hdu.data[49, 49] == 1000)
    assert(ai.sum == 1000)

    # do this one more time to check if the arithmetic is good
    ai.addPoints([49], [49], [1000])
    assert(ai.hdu.data[49, 49] == 2000)
    assert(ai.sum == 2000)

    # Test add a list of points
    ai = new_empty_astro_image(data_base)  # all data values = 0
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

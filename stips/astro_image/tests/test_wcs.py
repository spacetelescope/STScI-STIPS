import os, pytest

import numpy as np

from astropy.wcs import WCS

from stips.astro_image import AstroImage

@pytest.mark.veryslow
def test_astro_image_rotation(data_base):
    fits_path = os.path.join(data_base, "test", "wcs_test.fits")

    # Test 90 degree rotation
    ai = AstroImage.initFromFits(fits_path)
    default_wcs = ai._getWcs()  # WCS(ai.hdu)
    default_data = ai.hdu.data[:]

    ai.rotate(90)
    rotated_wcs = ai._getWcs()  # WCS(ai.hdu)
    rotated_data = ai.hdu.data[:]

    arr_size_x = ai.hdu.shape[0]
    arr_size_y = ai.hdu.shape[1]

    for x in range(arr_size_x):
        for y in range(arr_size_y):
            wcs_location = default_wcs.all_pix2world([[y, x]], 0)
            default_value = default_data[y][x]

            y, x = np.round(rotated_wcs.all_world2pix(wcs_location, 0)).astype(int)[0]
            assert (0 <= x < arr_size_x)
            assert (0 <= y < arr_size_y)
            rotated_value = rotated_data[y][x]

            if not np.isclose(default_value, rotated_value):
                print("Failed Values: ", default_value, rotated_value)
            assert np.isclose(default_value, rotated_value)

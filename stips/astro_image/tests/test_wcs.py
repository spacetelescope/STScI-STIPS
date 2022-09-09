from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import os
import pytest
import random

from stips.astro_image import AstroImage


@pytest.mark.veryslow
def test_astro_image_rotation_full(data_base):
    fits_path = os.path.join(data_base, "test", "wcs_test.fits")

    # Test 90 degree rotation
    ai = AstroImage.initFromFits(fits_path, psf=False)
    default_wcs = ai._getWcs()  # WCS(ai.hdu)
    default_data = ai.hdu.data[:]

    ai.rotate(90)
    rotated_wcs = ai._getWcs()  # WCS(ai.hdu)
    rotated_data = ai.hdu.data[:]

    arr_size_x = ai.hdu.shape[0]
    arr_size_y = ai.hdu.shape[1]

    for x in range(arr_size_x):
        print("Starting Iteration {} of {}".format(x, arr_size_x))
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


def test_astro_image_rotation_sample(data_base):
    fits_path = os.path.join(data_base, "test", "wcs_test.fits")

    # Test 90 degree rotation
    ai = AstroImage.initFromFits(fits_path, psf=False)
    default_wcs = ai._getWcs()  # WCS(ai.hdu)
    default_data = ai.hdu.data[:]

    ai.rotate(90)
    rotated_wcs = ai._getWcs()  # WCS(ai.hdu)
    rotated_data = ai.hdu.data[:]

    arr_size_x = ai.hdu.shape[0]
    arr_size_y = ai.hdu.shape[1]

    random.seed()

    coords = []

    for x in range(1000):
        coords.append((random.randrange(arr_size_x),
                       random.randrange(arr_size_y)))

    for x, y in coords:
        wcs_location = default_wcs.all_pix2world([[y, x]], 0)
        default_value = default_data[y][x]

        y, x = np.round(rotated_wcs.all_world2pix(wcs_location, 0)).astype(int)[0]
        assert (0 <= x < arr_size_x)
        assert (0 <= y < arr_size_y)
        rotated_value = rotated_data[y][x]

        if not np.isclose(default_value, rotated_value):
            print("Failed Values: ", default_value, rotated_value)
        assert np.isclose(default_value, rotated_value)


# RA, DEC, PA
wcs_data = [
     (18.5,   22.3, 27.9),
     (139.1,  -5.0,  0.0),
     (228.0, -55.8, 45.0)
]


@pytest.mark.parametrize(("ra", "dec", "pa"), wcs_data)
def test_fits_wcs(ra, dec, pa):
    input_image = AstroImage(ra=ra, dec=dec, pa=pa, psf=False)
    input_image.toFits('test_wcs.fits')

    pc_11 = np.cos(np.radians(pa))
    pc_12 = -np.sin(np.radians(pa))
    pc_21 = np.sin(np.radians(pa))
    pc_22 = np.cos(np.radians(pa))
    pc_matrix = [[pc_11, pc_12], [pc_21, pc_22]]

    with fits.open("test_wcs.fits") as input_file:
        img_wcs = WCS(input_file[0].header)
        assert np.isclose(img_wcs.wcs.crval[0], ra)
        assert np.isclose(img_wcs.wcs.crval[1], dec)
        assert np.isclose(img_wcs.wcs.pc[0][0], pc_matrix[0][0])
        assert np.isclose(img_wcs.wcs.pc[0][1], pc_matrix[0][1])
        assert np.isclose(img_wcs.wcs.pc[1][0], pc_matrix[1][0])
        assert np.isclose(img_wcs.wcs.pc[1][1], pc_matrix[1][1])
        assert np.isclose(input_file[0].header["CRVAL1"], ra)
        assert np.isclose(input_file[0].header["CRVAL2"], dec)
        assert np.isclose(input_file[0].header["PC1_1"], pc_matrix[0][0])
        assert np.isclose(input_file[0].header["PC1_2"], pc_matrix[0][1])
        assert np.isclose(input_file[0].header["PC2_1"], pc_matrix[1][0])
        assert np.isclose(input_file[0].header["PC2_2"], pc_matrix[1][1])

    os.remove("test_wcs.fits")

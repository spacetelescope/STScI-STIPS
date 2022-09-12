import numpy as np
import os
import pytest

from stips.astro_image import AstroImage
from stips import stips_data_base
from stips.utilities.testing import makeWCS, verifyData


currentdir = os.path.join(stips_data_base, "test", "astro_image")


testRA_data = [
    ([2., 4.], ["RA---TAN", "DEC--TAN"], 2., 4.,  -8., 352.),
    ([4., 2.], ["DEC--TAN", "RA---TAN"], 2., 4., 229., 229.)
]


@pytest.mark.parametrize(("coords", "coord_types", "ra", "dec", "new_ra", "out_ra"), testRA_data)
def test_RA(coords, coord_types, ra, dec, new_ra, out_ra):
    my_wcs = makeWCS(coords=coords, coord_types=coord_types)
    image = AstroImage(wcs=my_wcs, psf=False)
    assert image.ra == ra
    assert image.dec == dec
    image.ra = new_ra
    assert image.ra == out_ra


testDEC_data = [
    ([2., 4.], ["RA---TAN", "DEC--TAN"], 2., 4., -8., -8.),
    ([4., 2.], ["DEC--TAN", "RA---TAN"], 2., 4., 29., 29.)
]


@pytest.mark.parametrize(("coords", "coord_types", "ra", "dec", "new_dec", "out_dec"), testDEC_data)
def test_DEC(coords, coord_types, ra, dec, new_dec, out_dec):
    my_wcs = makeWCS(coords=coords, coord_types=coord_types)
    image = AstroImage(wcs=my_wcs, psf=False)
    assert image.ra == ra
    assert image.dec == dec
    image.dec = new_dec
    assert image.dec == out_dec


testPA_data = [
    (37.,   37.,  [.01, .17], -19.7, 340.3),
    (-4.5, 355.5, [.19, .22],  38.9,  38.9)
]


@pytest.mark.parametrize(("pa_in", "pa", "scale", "new_pa", "pa_out"), testPA_data)
def test_PA(pa_in, pa, scale, new_pa, pa_out):
    image = AstroImage(pa=pa_in, scale=scale, psf=False)
    assert abs(image.pa - pa) < 1.e-4
    assert abs(image.xscale - scale[0]*3600) < 1.e-4
    assert abs(image.yscale - scale[1]*3600) < 1.e-4
    image.pa = new_pa
    assert abs(image.pa - pa_out) < 1.e-4
    assert abs(image.xscale - scale[0]*3600) < 1.e-4
    assert abs(image.yscale - scale[1]*3600) < 1.e-4


updateHeader_data = [
    ([("TIMEINFO", 2.5), ("PLOTINFO", "TEST")],
     [("TIMEINFO", 2.5), ("PLOTINFO", "TEST")]),
    ([("TIMEINFO", 2.5), ("TIMEINFO", 4.5)],
     [("TIMEINFO", 4.5)])
]


@pytest.mark.parametrize(("header_in", "header_out"), updateHeader_data)
def test_updateHeader(header_in, header_out):
    image = AstroImage(psf=False)
    for k, v in header_in:
        image.updateHeader(k, v)
    for k, v in header_out:
        assert image.header[k] == v


addHistory_data = [
    ("ITEM1", "ITEM2", "ITEM3", "ITEM4"),
    ("TEST1", "TEST2", "TEST3", "TEST4")
]


@pytest.mark.parametrize(("history"), addHistory_data)
def test_addHistory(history):
    image = AstroImage(psf=False)
    for item in history:
        image.addHistory(item)
    for item in history:
        assert item in image.history
    for item in image.history:
        assert item in history


rotate_data = [
    (
        np.fromfile(os.path.join(currentdir, "rotate_01"), sep=' ').reshape(1000, 1000),
        39.8,
        np.fromfile(os.path.join(currentdir, "rotate_02"), sep=' ').reshape(1000, 1000)
    ),
    (
        np.fromfile(os.path.join(currentdir, "rotate_01"), sep=' ').reshape(1000, 1000),
        90.,
        np.fromfile(os.path.join(currentdir, "rotate_03"), sep=' ').reshape(1000, 1000)
    ),
    (
        np.fromfile(os.path.join(currentdir, "rotate_01"), sep=' ').reshape(1000, 1000),
        -28.9,
        np.fromfile(os.path.join(currentdir, "rotate_04"), sep=' ').reshape(1000, 1000)
    ),
    (
        np.fromfile(os.path.join(currentdir, "rotate_01"), sep=' ').reshape(1000, 1000),
        180.,
        np.fromfile(os.path.join(currentdir, "rotate_05"), sep=' ').reshape(1000, 1000)
    )
]


@pytest.mark.parametrize(("input", "angle", "output"), rotate_data)
def test_rotate(input, angle, output):
    # If psf=True, then the image will create a default PSF, and pad its
    #   borders (and thus increase its size) sufficiently to include off-image
    #   regions of half the PSF width on each side. That will, in turn, result
    #   in the image sizes not matching the pre-made output arrays. Because this
    #   test is not related to PSFs in any way, no PSF is created.
    image = AstroImage(data=input, psf=False)
    image.rotate(angle)
    verifyData(image.hdu.data, output)

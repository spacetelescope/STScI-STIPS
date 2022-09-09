import os

from stips.astro_image import AstroImage


def new_empty_astro_image(data_base):
    fits_path = os.path.join(data_base, "astro_image", "zero_image.fits")
    return AstroImage.initFromFits(fits_path, psf=False)


def test_add_rdnoise(data_base):
    # Test add one point
    ai = new_empty_astro_image(data_base)  # all data values = 0

    assert 'rdnoise' not in ai.hdu.header

    ai.introduceReadnoise(17.5)
    assert ai.hdu.header['rdnoise'] == 17.5

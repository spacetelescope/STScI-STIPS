import numpy as np
from astropy import wcs


def makeGaussian(size, fwhm=3, center=None):
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    return np.exp(-4*np.log(2)*((x-x0)**2 + (y-y0)**2)/fwhm**2)


def makeWCS(coords=[0.,0.],coord_types=["RA---TAN","DEC--TAN"],xsize=1,ysize=1,pa=0.,scale=[1.,1.],sip=None):
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = coord_types
    w.wcs.crpix = [int(np.floor(xsize/2.)),int(np.floor(ysize/2.))]
    w.wcs.crval = coords
    w.wcs.cdelt = np.array(scale) / 3600.
    cpa = np.cos(np.radians(pa%360.))
    spa = np.sin(np.radians(pa%360.))
    w.wcs.pc = np.array([[cpa,-spa],[spa,cpa]])
    if sip is not None:
        w.sip = sip
    return w


def verifyData(dat1,dat2):
    assert dat1.shape[0] == dat2.shape[0]
    assert dat1.shape[1] == dat2.shape[1]
#     np.testing.assert_allclose(dat1,dat2,atol=1e-3)
    for y in range(dat1.shape[0]):
        for x in range(dat1.shape[1]):
            verifyPoint(dat1,dat2,x,y)


def verifyPoint(dat1,dat2,x,y):
    assert abs(dat1[y,x] - dat2[y,x]) < 1e-3


def verifyImage(im1,im2):
    assert im1.out_path == im2.out_path
    assert im1.name == im2.name
    assert im1.xsize == im2.xsize
    assert im1.ysize == im2.ysize
    np.testing.assert_allclose((im1.xscale,im1.yscale),(im2.xscale,im2.yscale),atol=1e-3)
    assert im1.distorted == im2.distorted
    np.testing.assert_allclose((im1.ra,im1.dec,im1.pa),(im2.ra,im2.dec,im2.pa),atol=1e-3)
    assert im1.history == im2.history
    assert im1.zeropoint == im2.zeropoint
    assert im1.header['EQUINOX'] == im2.header['EQUINOX']
    np.testing.assert_allclose((im1.header['PA_APER']),(im2.header['PA_APER']),atol=1e-3)
    assert im1.header['VAFACTOR'] == im2.header['VAFACTOR']
    np.testing.assert_allclose((im1.header['ORIENTAT']),(im2.header['ORIENTAT']),atol=1e-3)
    np.testing.assert_allclose((im1.header['RA_APER']),(im2.header['RA_APER']),atol=1e-3)
    np.testing.assert_allclose((im1.header['DEC_APER']),(im2.header['DEC_APER']),atol=1e-3)
    assert im1.header['NAXIS1'] == im2.header['NAXIS1']
    assert im1.header['NAXIS2'] == im2.header['NAXIS2']


def verifyParameters(image,results):
    assert image.out_path == results['out_path']
    assert image.name == results['name']
    assert image.xsize == results['xsize']
    assert image.ysize == results['ysize']
    np.testing.assert_allclose((image.xscale,image.yscale),(results['xscale'],results['yscale']),atol=1e-3)
    assert image.distorted == results['distorted']
    np.testing.assert_allclose((image.ra,image.dec,image.pa),(results['ra'],results['dec'],results['pa']),atol=1e-3)
    assert image.history == results['history']
    assert image.zeropoint == results['zeropoint']
    assert image.header['EQUINOX'] == results['equinox']
    np.testing.assert_allclose((image.header['PA_APER']),(results['pa_aper']),atol=1e-3)
    assert image.header['VAFACTOR'] == results['vafactor']
    np.testing.assert_allclose((image.header['ORIENTAT']),(results['orientat']),atol=1e-3)
    np.testing.assert_allclose((image.header['RA_APER']),(results['ra_aper']),atol=1e-3)
    np.testing.assert_allclose((image.header['DEC_APER']),(results['dec_aper']),atol=1e-3)
    assert image.header['NAXIS1'] == results['naxis1']
    assert image.header['NAXIS2'] == results['naxis2']

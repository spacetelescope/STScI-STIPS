from __future__ import division

import os

import pytest

import numpy as np

import astropy.io.fits as pyfits
import astropy.wcs as wcs

from stips import AstroImage
from stips import stips_data_base
from stips.utilities.testing import (make_gaussian, make_wcs, verify_data,
                                     verify_image, verify_parameters,
                                     verify_point)

currentdir = os.path.join(stips_data_base, "test", "astro_image")

####################
# Setup Input Data #
####################

creation_data = [
    (
        {},
        {'out_path':os.getcwd(), 'name':"", 'xsize':1, 'ysize':1,
         'xscale':1.,'yscale':1., 'distorted':False, 'ra':0.,
         'dec':0., 'pa':0.,'history':[], 'zeropoint':0.,
         'equinox':2000.0, 'pa_aper':0.,'vafactor':0., 'orientat':0.,
         'ra_aper':0., 'dec_aper':0.,'naxis1':1, 'naxis2':1}
    ),
    (
        {'detname':'Thomas',
         'data':np.zeros((37,23)), 'scale':[0.37,7.52], 'ra':327.5,
         'dec':-38.4, 'pa':28.17},
        {'out_path':os.getcwd(), 'name':'Thomas', 'xsize':23, 'ysize':37,
         'xscale':0.37,'yscale':7.52, 'distorted':False, 'ra':327.5,
         'dec':-38.4, 'pa':28.17,'history':[], 'zeropoint':0.,
         'equinox':2000.0, 'pa_aper':28.17,'vafactor':0., 'orientat':28.17,
         'ra_aper':327.5, 'dec_aper':-38.4,'naxis1':23, 'naxis2':37}
    )
]


initFromFits_data = [
    (
        os.path.join(currentdir,"iabf01c5q_flt.fits"),
        {"ext":1},
        {'out_path':os.getcwd(), 'name':"", 'xsize':1014, 'ysize':1014,
         'xscale':0.1354,'yscale':0.1210, 'distorted':True, 'ra':5.6604,
         'dec':-72.0678, 'pa':149.6388,
         'history':["Created from FITS file iabf01c5q_flt.fits"],
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':149.6388,'vafactor':0.,
         'orientat':149.6388, 'ra_aper':5.6604, 'dec_aper':-72.0678,
         'naxis1':1014, 'naxis2':1014}
    ),
    (
        os.path.join(currentdir,"iabf01c5q_flt.fits"),
        {"ext":1,"zeropoint":-21.1},
                {'out_path':os.getcwd(), 'name':"", 'xsize':1014, 'ysize':1014,
                 'xscale':0.1354,'yscale':0.1210, 'distorted':True, 'ra':5.6604,
                 'dec':-72.0678, 'pa':149.6388,
                 'history':["Created from FITS file iabf01c5q_flt.fits"],
                 'zeropoint':-21.1, 'equinox':2000.0, 'pa_aper':149.6388,'vafactor':0.,
                 'orientat':149.6388, 'ra_aper':5.6604, 'dec_aper':-72.0678,
                 'naxis1':1014, 'naxis2':1014}
    )
]


initDataFromFits_data = [ 
    (
        os.path.join(currentdir,"iabf01c5q_flt.fits"),
        {"ext":1},
        {'out_path':os.getcwd(), 'name':"", 'xsize':1014, 'ysize':1014, 'xscale':1.,
         'yscale':1., 'distorted':False, 'ra':0., 'dec':0., 'pa':0.,
         'history':["Data imported from FITS file iabf01c5q_flt.fits"], 
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':0.,'vafactor':0., 
         'orientat':0., 'ra_aper':0., 'dec_aper':0., 'naxis1':1014, 'naxis2':1014}
    ),
    (
        os.path.join(currentdir,"iabf01c5q_flt.fits"),
        {"ext":1, "scale":[0.1354, 0.1210], "ra":5.6604, "dec":-72.0678,
         "pa":146.742, "zeropoint":-21.1},
        {'out_path':os.getcwd(), 'name':"", 'xsize':1014, 'ysize':1014, 
         'xscale':0.1354,'yscale':0.1210, 'distorted':False, 'ra':5.6604,
         'dec':-72.0678, 'pa':146.742,
         'history':["Data imported from FITS file iabf01c5q_flt.fits"], 
         'zeropoint':-21.1, 'equinox':2000.0, 'pa_aper':146.742,'vafactor':0., 
         'orientat':146.742, 'ra_aper':5.6604, 'dec_aper':-72.0678,
         'naxis1':1014, 'naxis2':1014}
    )
]


initFromPoints_data = [
    (
        np.array((3.,7.,9.,17.,22)),
        np.array((0.,2.,7.5,12.,14.)),
        np.array((3.,1.3e-7,2.9e5,0.2,0.7)),
        [23,15],
        {},
        {'out_path':os.getcwd(), 'name':"", 'xsize':23, 'ysize':15, 
         'xscale':1.,'yscale':1., 'distorted':False, 'ra':0., 
         'dec':0., 'pa':0.,'history':['Adding 5 point sources'], 
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':0.,'vafactor':0., 
         'orientat':0., 'ra_aper':0., 'dec_aper':0.,'naxis1':23, 'naxis2':15}
    ),
    (
        np.array((3.,7.,9.,17.,22)),
        np.array((0.,2.,7.5,12.,14.)),
        np.array((3.,1.3e-7,2.9e5,0.2,0.7)),
        [30,30],
        {'xsize':30,'ysize':30,'dec':6.2,'pa':12.22,'scale':[0.11,0.15]},
        {'out_path':os.getcwd(), 'name':"", 'xsize':30, 'ysize':30, 
         'xscale':0.11,'yscale':0.15, 'distorted':False, 'ra':0., 
         'dec':6.2, 'pa':12.22,'history':['Adding 5 point sources'], 
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':12.22,'vafactor':0., 
         'orientat':12.22, 'ra_aper':0., 'dec_aper':6.2,'naxis1':30, 'naxis2':30}
    ),
    (
        np.array((3.,7.,9.,17.,22)),
        np.array((0.,2.,7.5,12.,14.)),
        np.array((3.,1.3e-7,2.9e5,0.2,0.7)),
        [23,30],
        {'xsize':10,'ysize':30,'dec':6.2,'pa':12.22,'scale':[0.11,0.15]},
        {'out_path':os.getcwd(), 'name':"", 'xsize':23, 'ysize':30, 
         'xscale':0.11,'yscale':0.15, 'distorted':False, 'ra':0., 
         'dec':6.2, 'pa':12.22,'history':['Adding 5 point sources'], 
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':12.22,'vafactor':0., 
         'orientat':12.22, 'ra_aper':0., 'dec_aper':6.2,'naxis1':23, 'naxis2':30}
    ),
    (
        np.array((3.,7.,9.,17.,22)),
        np.array((0.,2.,7.5,12.,14.)),
        np.array((3.,1.3e-7,2.9e5,0.2,0.7)),
        [30,15],
        {'xsize':30,'ysize':10,'dec':6.2,'pa':12.22,'scale':[0.11,0.15]},
        {'out_path':os.getcwd(), 'name':"", 'xsize':30, 'ysize':15, 
         'xscale':0.11,'yscale':0.15, 'distorted':False, 'ra':0., 
         'dec':6.2, 'pa':12.22,'history':['Adding 5 point sources'], 
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':12.22,'vafactor':0., 
         'orientat':12.22, 'ra_aper':0., 'dec_aper':6.2,'naxis1':30, 'naxis2':15}
    )
]



initFromProfile_data = [
    (
        17.8,-3.1,7.7,0.25,3.5,28.7,0.3,
        np.array((0.0060756)).reshape(1,1),
        {},
        {'out_path':os.getcwd(), 'name':"", 'xsize':1, 'ysize':1,
         'xscale':1.,'yscale':1., 'distorted':False, 'ra':0.,
         'dec':0., 'pa':0.,'zeropoint':0., 'equinox':2000.0, 'pa_aper':0.,
         'history':['Adding Sersic profile at (17.800000,-3.100000) with flux 7.700000, index 0.250000, Re 3.500000, Phi 28.700000, and axial ratio 0.300000'],
         'vafactor':0., 'orientat':0., 'ra_aper':0., 'dec_aper':0.,
         'naxis1':1, 'naxis2':1}
    ),
    (
        17.8,5.1,7.7,0.25,3.5,38.7,0.3,
        np.fromfile(os.path.join(currentdir,"initFromProfile_01"),sep=' ').reshape(25,30),
        {'xsize':30.,'ysize':25.},
        {'out_path':os.getcwd(), 'name':"", 'xsize':30, 'ysize':25,
         'xscale':1.,'yscale':1., 'distorted':False, 'ra':0.,
         'dec':0., 'pa':0.,'history':['Adding Sersic profile at (17.800000,5.100000) with flux 7.700000, index 0.250000, Re 3.500000, Phi 38.700000, and axial ratio 0.300000'],
         'zeropoint':0., 'equinox':2000.0, 'pa_aper':0.,'vafactor':0.,
         'orientat':0., 'ra_aper':0., 'dec_aper':0.,'naxis1':30, 'naxis2':25}
    )
]


testRA_data = [
    ([2.,4.],["RA---TAN","DEC--TAN"],2.,4., -8., 352.),
    ([4.,2.],["DEC--TAN","RA---TAN"],2.,4.,229.,229.)
]


testDEC_data = [
    ([2.,4.],["RA---TAN","DEC--TAN"],2.,4., -8., -8.),
    ([4.,2.],["DEC--TAN","RA---TAN"],2.,4.,29.,29.)
]


testPA_data = [
    (37.,37.,[.01,.17],-19.7,340.3),
    (-4.5,355.5,[.19,.22],38.9,38.9)
]



updateHeader_data = [
    ([("TIMEINFO",2.5),("PLOTINFO","TEST")],
     [("TIMEINFO",2.5),("PLOTINFO","TEST")]),
    ([("TIMEINFO",2.5),("TIMEINFO",4.5)],
     [("TIMEINFO",4.5)])
]

addHistory_data = [
    ("ITEM1","ITEM2","ITEM3","ITEM4"),
    ("TEST1","TEST2","TEST3","TEST4")
]

addCatalogue_data =   [
    (
        {},
        os.path.join(currentdir,"cat_test.txt"),
        np.array([0.]).reshape(1,1)
    ),
    (
        {'xsize':1000,'ysize':1000},
        os.path.join(currentdir,"cat_test.txt"),
        np.fromfile(os.path.join(currentdir,"addCatalogue_01"),sep=' ').reshape(1000,1000)
    )
]

convolve_data = [
(
    np.fromfile(os.path.join(currentdir,"convolve_01"),sep=' ').reshape(1000,1000),
    make_gaussian(9,3),
    np.fromfile(os.path.join(currentdir,"convolve_02"),sep=' ').reshape(1000,1000)
    ),
    (
    np.fromfile(os.path.join(currentdir,"convolve_01"),sep=' ').reshape(1000,1000),
    np.fromfile(os.path.join(currentdir,"convolve_04"),sep=' ').reshape(1000,1000),
    np.fromfile(os.path.join(currentdir,"convolve_03"),sep=' ').reshape(1000,1000)
)
]

addWithAlignment_data =   [
    (
        np.zeros((10,10),dtype='float32'),
        np.arange(100,dtype='float32').reshape(10,10),
        {},
        np.arange(100,dtype='float32').reshape(10,10)
    ),
    (
        np.zeros((10,10),dtype='float32'),
        np.arange(100,dtype='float32').reshape(10,10),
        {'dec':(5./3600.)},
        np.fromfile(os.path.join(currentdir,"addWithAlignment_01"),sep=' ').reshape(10,10)
    ),
    (
        np.zeros((10,10),dtype='float32'),
        np.arange(100,dtype='float32').reshape(10,10),
        {'pa':90.},
        np.fromfile(os.path.join(currentdir,"addWithAlignment_02"),sep=' ').reshape(10,10)
    ),
    (
        np.zeros((10,10),dtype='float32'),
        np.arange(100,dtype='float32').reshape(10,10),
        {'pa':38.23},
        np.fromfile(os.path.join(currentdir,"addWithAlignment_03"),sep=' ').reshape(10,10)
    ),
    (
        np.zeros((10,10),dtype='float32'),
        np.arange(100,dtype='float32').reshape(10,10),
        {'pa':27.5,'ra':(-5./3600.)},
        np.fromfile(os.path.join(currentdir,"addWithAlignment_04"),sep=' ').reshape(10,10)
    )
]


rescale_data =   [
    (
        np.full((10,10),3.),
        np.array([0.5,0.5]),
        np.full((20,20),.75)
    ),
    (
        np.full((10,10),3.),
        np.array([2.,2.]),
        np.full((5,5),12.)
    ),
    (
        np.full((10,10),3.),
        np.array([1.,2.]),
        np.full((10,5),6.)
    ),
    (
        np.full((10,10),3.),
        np.array([1.25,.75]),
        np.full((8,13),2.8846)
    )
]

crop_data =   [
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [5.,5.], [0,0],
        np.arange(100,dtype='float32').reshape(10,10)[3:8,3:8]
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [3.,3.], [4,4],
        np.arange(100,dtype='float32').reshape(10,10)[:3,:3]
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [12.,12.], [0,0],
        np.fromfile(os.path.join(currentdir,"crop_01"),sep=' ').reshape(12,12)
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [12.,12.], [4,4],
        np.fromfile(os.path.join(currentdir,"crop_02"),sep=' ').reshape(12,12)
    )
]
bin_data =   [
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [2,2],
        np.fromfile(os.path.join(currentdir,"bin_01"),sep=' ').reshape(5,5)
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [3,3],
        np.fromfile(os.path.join(currentdir,"bin_02"),sep=' ').reshape(3,3)
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        [2,3],
        np.fromfile(os.path.join(currentdir,"bin_03"),sep=' ').reshape(3,5)
    ),
    (
        np.arange(90,dtype='float32').reshape(10,9),
        [3,2],
        np.fromfile(os.path.join(currentdir,"bin_04"),sep=' ').reshape(5,3)
    )
]

rotate_data =   [
    (
        np.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
        39.8,
        np.fromfile(os.path.join(currentdir,"rotate_02"),sep=' ').reshape(1000,1000)
    ),
    (
        np.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
        90.,
        np.fromfile(os.path.join(currentdir,"rotate_03"),sep=' ').reshape(1000,1000)
    ),
    (
        np.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
        -28.9,
        np.fromfile(os.path.join(currentdir,"rotate_04"),sep=' ').reshape(1000,1000)
    ),
    (
        np.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
        180.,
        np.fromfile(os.path.join(currentdir,"rotate_05"),sep=' ').reshape(1000,1000)
    )
]
addWithOffset_data =   [
    (
        np.arange(100,dtype='float32').reshape(10,10),
        np.full((10,10),3.),
        (4,-2),
        {},
        np.fromfile(os.path.join(currentdir,"addWithOffset_01"),sep=' ').reshape(10,10)
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        np.full((10,10),3.),
        (0,0),
        {'scale':[2.,2.]},
        np.fromfile(os.path.join(currentdir,"addWithOffset_02"),sep=' ').reshape(10,10)
    ),
    (
        np.arange(100,dtype='float32').reshape(10,10),
        np.full((10,10),3.),
        (0,0),
        {'scale':[0.5,0.5]},
        np.fromfile(os.path.join(currentdir,"addWithOffset_03"),sep=' ').reshape(10,10)
    )
]


#########
# Tests #
#########

@pytest.mark.parametrize(("inputs","results"), creation_data)
def test_creation(inputs,results):
    image = AstroImage(**inputs)
    data = np.zeros((image.ysize, image.xsize))
    verify_parameters(image, results)
    verify_data(image.data, data)



@pytest.mark.parametrize(("file","inputs","results"), initFromFits_data)
def test_init_from_fits(file,inputs,results):
    image = AstroImage.initFromFits(file,**inputs)
    fits = pyfits.open(file)
    data = fits[inputs["ext"]].data
    fits.close()
    verify_parameters(image, results)
    verify_data(image.data, data)


@pytest.mark.parametrize(("file","inputs","results"), initDataFromFits_data)
def test_init_data_from_fits(file,inputs,results):
    image = AstroImage.initDataFromFits(file,**inputs)
    fits = pyfits.open(file)
    data = fits[inputs["ext"]].data
    fits.close()
    verify_parameters(image, results)
    verify_data(image.data, data)


@pytest.mark.parametrize(("xs","ys","rates","datasize","inputs","results"), initFromPoints_data)
def test_init_from_points(xs,ys,rates,datasize,inputs,results):
    image = AstroImage.initFromPoints(xs,ys,rates,**inputs)
    verify_parameters(image,results)
    dat = np.zeros((datasize[1],datasize[0]),dtype='float32')
    dat[0,3] = 3.
    dat[2,7] = 1.3e-7
    dat[7,9] = 2.9e5
    dat[12,17] = 0.2
    dat[14,22] = 0.7
    verify_data(image.data, dat)


@pytest.mark.parametrize(("x","y","f","n","re","phi","ar","dat","inputs","results"), initFromProfile_data)
def test_init_from_profile(x,y,f,n,re,phi,ar,dat,inputs,results):
    image = AstroImage.initFromProfile(x,y,f,n,re,phi,ar,**inputs)
    verify_parameters(image, results)
    verify_data(image.data, dat)


@pytest.mark.parametrize(("coords","coord_types","ra","dec","new_ra","out_ra"), testRA_data)
def test_ra(coords,coord_types,ra,dec,new_ra,out_ra):
    my_wcs = make_wcs(coords=coords, coord_types=coord_types)
    image = AstroImage(wcs=my_wcs)
    assert image.ra == ra
    assert image.dec == dec
    image.ra = new_ra
    assert image.ra == out_ra


@pytest.mark.parametrize(("coords","coord_types","ra","dec","new_dec","out_dec"), testDEC_data)
def test_dec(coords,coord_types,ra,dec,new_dec,out_dec):
    my_wcs = make_wcs(coords=coords,coord_types=coord_types)
    image = AstroImage(wcs=my_wcs)
    assert image.ra == ra
    assert image.dec == dec
    image.dec = new_dec
    assert image.dec == out_dec


@pytest.mark.parametrize(("pa_in","pa","scale","new_pa","pa_out"), testPA_data)
def test_pa(pa_in,pa,scale,new_pa,pa_out):
    image = AstroImage(pa=pa_in,scale=scale)
    assert abs(image.pa - pa) < 1.e-4
    assert abs(image.xscale - scale[0]) < 1.e-4
    assert abs(image.yscale - scale[1]) < 1.e-4
    image.pa = new_pa
    assert abs(image.pa - pa_out) < 1.e-4
    assert abs(image.xscale - scale[0]) < 1.e-4
    assert abs(image.yscale - scale[1]) < 1.e-4


@pytest.mark.parametrize(("header_in","header_out"), updateHeader_data)
def test_update_header(header_in,header_out):
    image = AstroImage()
    for k,v in header_in:
        image.updateHeader(k,v)
    for k,v in header_out:
        assert image.header[k] == v


@pytest.mark.parametrize("history", addHistory_data)
def test_add_history(history):
    image = AstroImage()
    for item in history:
        image.addHistory(item)
    for item in history:
        assert item in image.history
    for item in image.history:
        assert item in history


@pytest.mark.parametrize(("inputs","catalogue","final"), addCatalogue_data)
def test_add_catalogue(inputs, catalogue, final):
    image = AstroImage(**inputs)
    out_cat = image.addCatalogue(catalogue)
    if os.path.exists(out_cat):
        os.remove(out_cat)
    verify_data(image.data,final)


@pytest.mark.parametrize(("input","kernel","result"), convolve_data)
def test_convolve(input, kernel, result):
    input_image = AstroImage(data=input)
    kernel_image = AstroImage(data=kernel)
    input_image.convolve(kernel_image)
    verify_data(input_image.data,result)


@pytest.mark.parametrize(("input","angle","output"), rotate_data)
def test_rotate(input, angle, output):
    image = AstroImage(data=input)
    image.rotate(angle)
    verify_data(image.data,output)


@pytest.mark.parametrize(("in1","in2","offset","inputs","result"), addWithOffset_data)
def test_add_with_offset(in1, in2, offset, inputs, result):
    im1 = AstroImage(data=in1,**inputs)
    im2 = AstroImage(data=in2)
    im1.addWithOffset(im2,offset[0],offset[1])
    verify_data(im1.data,result)


@pytest.mark.parametrize(("in1","in2","inputs","result"), addWithAlignment_data)
def test_add_with_alignment(in1, in2, inputs, result):
    im1 = AstroImage(data=in1)
    im2 = AstroImage(data=in2,**inputs)
    im1.addWithAlignment(im2)
    verify_data(im1.data,result)


@pytest.mark.parametrize(("input","scale","result"), rescale_data)
def test_rescale(input, scale, result):
    im1 = AstroImage(data=input)
    im1.rescale(scale)
    verify_data(im1.data,result)
    im1.rescale([1.,1.])
    verify_data(im1.data,input)


@pytest.mark.parametrize(("input","size","offset","result"), crop_data)
def test_crop(input, size, offset, result):
    im1 = AstroImage(data=input)
    im1.crop(size[0],size[1],offset[0],offset[1])
    verify_data(im1.data,result)


@pytest.mark.parametrize(("input","bin","result"), bin_data)
def test_bin(input, bin, result):
    im1 = AstroImage(data=input)
    im1.bin(bin[0],bin[1])
    verify_data(im1.data,result)


#######################
# Saved Test Snippets #
#######################
#     @property
#     def hdu(self):
#         """Output AstroImage as a FITS Primary HDU"""
#         self.logger.info("Creating Primary HDU from AstroImage %s",self.name)
#         hdu = pyfits.PrimaryHDU(self.data,header=self.wcs.to_header())
#         for k,v in self.header.iteritems():
#             hdu.header[k] = v
#         for item in self.history:
#             hdu.header.add_history(item)
#         self.logger.info("Created Primary HDU from AstroImage %s",self.name)
#         return hdu
#
#     @property
#     def imageHdu(self):
#         """Output AstroImage as a FITS Extension HDU"""
#         self.logger.info("Creating Extension HDU from AstroImage %s",self.name)
#         hdu = pyfits.ImageHDU(self.data,header=self.wcs.to_header(),name=self.name)
#         for k,v in self.header.iteritems():
#             hdu.header[k] = v
#         for item in self.history:
#             hdu.header.add_history(item)
#         self.logger.info("Created Extension HDU from AstroImage %s",self.name)
#         return hdu
#
#     def toFits(self,outFile):
#         """Create a FITS file from the current state of the AstroImage data."""
#         self.logger.info("Writing AstroImage %s to FITS",self.name)
#         hdulist = pyfits.HDUList([self.hdu])
#         hdulist.writeto(outFile)

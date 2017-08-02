from __future__ import division

import astropy.io.fits as pyfits
import astropy.wcs as wcs

from AstroImage import AstroImage

def makeGaussian(size, fwhm=3, center=None):
    x = numpy.arange(0, size, 1, float)
    y = x[:,numpy.newaxis]
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    return numpy.exp(-4*numpy.log(2)*((x-x0)**2 + (y-y0)**2)/fwhm**2)

def makeWCS(coords=[0.,0.],coord_types=["RA---TAN","DEC--TAN"],xsize=1,ysize=1,pa=0.,scale=[1.,1.],sip=None):
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = coord_types
    w.wcs.crpix = [int(numpy.floor(xsize/2.)),int(numpy.floor(ysize/2.))]
    w.wcs.crval = coords
    w.wcs.cdelt = numpy.array(scale) / 3600.
    cpa = numpy.cos(numpy.radians(pa%360.))
    spa = numpy.sin(numpy.radians(pa%360.))
    w.wcs.pc = numpy.array([[cpa,-spa],[spa,cpa]])
    if sip is not None:
        w.sip = sip
    return w

def verifyData(dat1,dat2):
    assert dat1.shape == dat2.shape
#     numpy.testing.assert_allclose(dat1,dat2,atol=1e-3)
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
    numpy.testing.assert_allclose((im1.xscale,im1.yscale),(im2.xscale,im2.yscale),atol=1e-3)
    assert im1.distorted == im2.distorted
    numpy.testing.assert_allclose((im1.ra,im1.dec,im1.pa),(im2.ra,im2.dec,im2.pa),atol=1e-3)
    assert im1.history == im2.history
    assert im1.zeropoint == im2.zeropoint
    assert im1.header['EQUINOX'] == im2.header['EQUINOX']
    numpy.testing.assert_allclose((im1.header['PA_APER']),(im2.header['PA_APER']),atol=1e-3)
    assert im1.header['VAFACTOR'] == im2.header['VAFACTOR']
    numpy.testing.assert_allclose((im1.header['ORIENTAT']),(im2.header['ORIENTAT']),atol=1e-3)
    numpy.testing.assert_allclose((im1.header['RA_APER']),(im2.header['RA_APER']),atol=1e-3)
    numpy.testing.assert_allclose((im1.header['DEC_APER']),(im2.header['DEC_APER']),atol=1e-3)
    assert im1.header['NAXIS1'] == im2.header['NAXIS1']
    assert im1.header['NAXIS2'] == im2.header['NAXIS2']

def verifyParameters(image,results):
    assert image.out_path == results['out_path']
    assert image.name == results['name']
    assert image.xsize == results['xsize']
    assert image.ysize == results['ysize']
    numpy.testing.assert_allclose((image.xscale,image.yscale),(results['xscale'],results['yscale']),atol=1e-3)
    assert image.distorted == results['distorted']
    numpy.testing.assert_allclose((image.ra,image.dec,image.pa),(results['ra'],results['dec'],results['pa']),atol=1e-3)
    assert image.history == results['history']
    assert image.zeropoint == results['zeropoint']
    assert image.header['EQUINOX'] == results['equinox']
    numpy.testing.assert_allclose((image.header['PA_APER']),(results['pa_aper']),atol=1e-3)
    assert image.header['VAFACTOR'] == results['vafactor']
    numpy.testing.assert_allclose((image.header['ORIENTAT']),(results['orientat']),atol=1e-3)
    numpy.testing.assert_allclose((image.header['RA_APER']),(results['ra_aper']),atol=1e-3)
    numpy.testing.assert_allclose((image.header['DEC_APER']),(results['dec_aper']),atol=1e-3)
    assert image.header['NAXIS1'] == results['naxis1']
    assert image.header['NAXIS2'] == results['naxis2']

creation_data =         [
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
                                 'data':numpy.zeros((37,23)), 'scale':[0.37,7.52], 'ra':327.5, 
                                 'dec':-38.4, 'pa':28.17},
                                {'out_path':os.getcwd(), 'name':'Thomas', 'xsize':23, 'ysize':37, 
                                 'xscale':0.37,'yscale':7.52, 'distorted':False, 'ra':327.5, 
                                 'dec':-38.4, 'pa':28.17,'history':[], 'zeropoint':0., 
                                 'equinox':2000.0, 'pa_aper':28.17,'vafactor':0., 'orientat':28.17, 
                                 'ra_aper':327.5, 'dec_aper':-38.4,'naxis1':23, 'naxis2':37}
                            )
                        ]

@pytest.mark.parametrize(("inputs","results"), creation_data)
def test_creation(inputs,results):
    image = AstroImage(**inputs)
    data = numpy.zeros((image.ysize,image.xsize))
    verifyParameters(image,results)
    verifyData(image.data,data)

initFromFits_data = [ (
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

@pytest.mark.parametrize(("file","inputs","results"), initFromFits_data)
def test_initFromFits(file,inputs,results):
    image = AstroImage.initFromFits(file,**inputs)
    fits = pyfits.open(file)
    data = fits[inputs["ext"]].data
    fits.close()
    verifyParameters(image,results)
    verifyData(image.data,data)

initDataFromFits_data =     [ (
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

@pytest.mark.parametrize(("file","inputs","results"), initDataFromFits_data)
def test_initDataFromFits(file,inputs,results):
    image = AstroImage.initDataFromFits(file,**inputs)
    fits = pyfits.open(file)
    data = fits[inputs["ext"]].data
    fits.close()
    verifyParameters(image,results)
    verifyData(image.data,data)

initFromPoints_data =   [
                            (
                                numpy.array((3.,7.,9.,17.,22)),
                                numpy.array((0.,2.,7.5,12.,14.)),
                                numpy.array((3.,1.3e-7,2.9e5,0.2,0.7)),
                                [23,15],
                                {},
                                {'out_path':os.getcwd(), 'name':"", 'xsize':23, 'ysize':15, 
                                 'xscale':1.,'yscale':1., 'distorted':False, 'ra':0., 
                                 'dec':0., 'pa':0.,'history':['Adding 5 point sources'], 
                                 'zeropoint':0., 'equinox':2000.0, 'pa_aper':0.,'vafactor':0., 
                                 'orientat':0., 'ra_aper':0., 'dec_aper':0.,'naxis1':23, 'naxis2':15}
                            ),
                            (
                                numpy.array((3.,7.,9.,17.,22)),
                                numpy.array((0.,2.,7.5,12.,14.)),
                                numpy.array((3.,1.3e-7,2.9e5,0.2,0.7)),
                                [30,30],
                                {'xsize':30,'ysize':30,'dec':6.2,'pa':12.22,'scale':[0.11,0.15]},
                                {'out_path':os.getcwd(), 'name':"", 'xsize':30, 'ysize':30, 
                                 'xscale':0.11,'yscale':0.15, 'distorted':False, 'ra':0., 
                                 'dec':6.2, 'pa':12.22,'history':['Adding 5 point sources'], 
                                 'zeropoint':0., 'equinox':2000.0, 'pa_aper':12.22,'vafactor':0., 
                                 'orientat':12.22, 'ra_aper':0., 'dec_aper':6.2,'naxis1':30, 'naxis2':30}
                            ),
                            (
                                numpy.array((3.,7.,9.,17.,22)),
                                numpy.array((0.,2.,7.5,12.,14.)),
                                numpy.array((3.,1.3e-7,2.9e5,0.2,0.7)),
                                [23,30],
                                {'xsize':10,'ysize':30,'dec':6.2,'pa':12.22,'scale':[0.11,0.15]},
                                {'out_path':os.getcwd(), 'name':"", 'xsize':23, 'ysize':30, 
                                 'xscale':0.11,'yscale':0.15, 'distorted':False, 'ra':0., 
                                 'dec':6.2, 'pa':12.22,'history':['Adding 5 point sources'], 
                                 'zeropoint':0., 'equinox':2000.0, 'pa_aper':12.22,'vafactor':0., 
                                 'orientat':12.22, 'ra_aper':0., 'dec_aper':6.2,'naxis1':23, 'naxis2':30}
                            ),
                            (
                                numpy.array((3.,7.,9.,17.,22)),
                                numpy.array((0.,2.,7.5,12.,14.)),
                                numpy.array((3.,1.3e-7,2.9e5,0.2,0.7)),
                                [30,15],
                                {'xsize':30,'ysize':10,'dec':6.2,'pa':12.22,'scale':[0.11,0.15]},
                                {'out_path':os.getcwd(), 'name':"", 'xsize':30, 'ysize':15, 
                                 'xscale':0.11,'yscale':0.15, 'distorted':False, 'ra':0., 
                                 'dec':6.2, 'pa':12.22,'history':['Adding 5 point sources'], 
                                 'zeropoint':0., 'equinox':2000.0, 'pa_aper':12.22,'vafactor':0., 
                                 'orientat':12.22, 'ra_aper':0., 'dec_aper':6.2,'naxis1':30, 'naxis2':15}
                            )
                        ]

@pytest.mark.parametrize(("xs","ys","rates","datasize","inputs","results"), initFromPoints_data)
def test_initFromPoints(xs,ys,rates,datasize,inputs,results):
    image = AstroImage.initFromPoints(xs,ys,rates,**inputs)
    verifyParameters(image,results)
    dat = numpy.zeros((datasize[1],datasize[0]),dtype='float32')
    dat[0,3] = 3.
    dat[2,7] = 1.3e-7
    dat[7,9] = 2.9e5
    dat[12,17] = 0.2
    dat[14,22] = 0.7
    verifyData(image.data,dat)

initFromProfile_data =   [
                            (
                                17.8,-3.1,7.7,0.25,3.5,28.7,0.3,
                                numpy.array((0.0060756)).reshape(1,1),
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
                                numpy.fromfile(os.path.join(currentdir,"initFromProfile_01"),sep=' ').reshape(25,30),
                                {'xsize':30.,'ysize':25.},
                                {'out_path':os.getcwd(), 'name':"", 'xsize':30, 'ysize':25, 
                                 'xscale':1.,'yscale':1., 'distorted':False, 'ra':0., 
                                 'dec':0., 'pa':0.,'history':['Adding Sersic profile at (17.800000,5.100000) with flux 7.700000, index 0.250000, Re 3.500000, Phi 38.700000, and axial ratio 0.300000'], 
                                 'zeropoint':0., 'equinox':2000.0, 'pa_aper':0.,'vafactor':0., 
                                 'orientat':0., 'ra_aper':0., 'dec_aper':0.,'naxis1':30, 'naxis2':25}
                            )
                        ]

@pytest.mark.parametrize(("x","y","f","n","re","phi","ar","dat","inputs","results"), initFromProfile_data)
def test_initFromProfile(x,y,f,n,re,phi,ar,dat,inputs,results):
    image = AstroImage.initFromProfile(x,y,f,n,re,phi,ar,**inputs)
    verifyParameters(image,results)
    verifyData(image.data,dat)

testRA_data = [    ([2.,4.],["RA---TAN","DEC--TAN"],2.,4., -8., 352.),
                    ([4.,2.],["DEC--TAN","RA---TAN"],2.,4.,229.,229.)
                ]

@pytest.mark.parametrize(("coords","coord_types","ra","dec","new_ra","out_ra"), testRA_data)
def test_RA(coords,coord_types,ra,dec,new_ra,out_ra):
    my_wcs = makeWCS(coords=coords,coord_types=coord_types)
    image = AstroImage(wcs=my_wcs)
    assert image.ra == ra
    assert image.dec == dec
    image.ra = new_ra
    assert image.ra == out_ra

testDEC_data = [    ([2.,4.],["RA---TAN","DEC--TAN"],2.,4., -8., -8.),
                    ([4.,2.],["DEC--TAN","RA---TAN"],2.,4.,29.,29.)
                ]

@pytest.mark.parametrize(("coords","coord_types","ra","dec","new_dec","out_dec"), testDEC_data)
def test_DEC(coords,coord_types,ra,dec,new_dec,out_dec):
    my_wcs = makeWCS(coords=coords,coord_types=coord_types)
    image = AstroImage(wcs=my_wcs)
    assert image.ra == ra
    assert image.dec == dec
    image.dec = new_dec
    assert image.dec == out_dec

testPA_data =   [   
                    (37.,37.,[.01,.17],-19.7,340.3),
                    (-4.5,355.5,[.19,.22],38.9,38.9)
                ]

@pytest.mark.parametrize(("pa_in","pa","scale","new_pa","pa_out"), testPA_data)
def test_PA(pa_in,pa,scale,new_pa,pa_out):
    image = AstroImage(pa=pa_in,scale=scale)
    assert abs(image.pa - pa) < 1.e-4
    assert abs(image.xscale - scale[0]) < 1.e-4
    assert abs(image.yscale - scale[1]) < 1.e-4
    image.pa = new_pa
    assert abs(image.pa - pa_out) < 1.e-4
    assert abs(image.xscale - scale[0]) < 1.e-4
    assert abs(image.yscale - scale[1]) < 1.e-4

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

updateHeader_data = [
                        ([("TIMEINFO",2.5),("PLOTINFO","TEST")],
                         [("TIMEINFO",2.5),("PLOTINFO","TEST")]),
                        ([("TIMEINFO",2.5),("TIMEINFO",4.5)],
                         [("TIMEINFO",4.5)])
                    ]

@pytest.mark.parametrize(("header_in","header_out"), updateHeader_data)
def test_updateHeader(header_in,header_out):
    image = AstroImage()
    for k,v in header_in:
        image.updateHeader(k,v)
    for k,v in header_out:
        assert image.header[k] == v

addHistory_data = [
                    ("ITEM1","ITEM2","ITEM3","ITEM4"),
                    ("TEST1","TEST2","TEST3","TEST4")
                  ]

@pytest.mark.parametrize(("history"), addHistory_data)
def test_addHistory(history):
    image = AstroImage()
    for item in history:
        image.addHistory(item)
    for item in history:
        assert item in image.history
    for item in image.history:
        assert item in history

addCatalogue_data =   [
                            (
                                {},
                                os.path.join(currentdir,"cat_test.txt"),
                                numpy.array([0.]).reshape(1,1)
                            ),
                            (
                                {'xsize':1000,'ysize':1000},
                                os.path.join(currentdir,"cat_test.txt"),
                                numpy.fromfile(os.path.join(currentdir,"addCatalogue_01"),sep=' ').reshape(1000,1000)
                            )
                        ]

@pytest.mark.parametrize(("inputs","catalogue","final"), addCatalogue_data)
def test_addCatalogue(inputs,catalogue,final):
    image = AstroImage(**inputs)
    out_cat = image.addCatalogue(catalogue)
    if os.path.exists(out_cat):
        os.remove(out_cat)
    verifyData(image.data,final)

convolve_data = [
                    (
                        numpy.fromfile(os.path.join(currentdir,"convolve_01"),sep=' ').reshape(1000,1000),
                        makeGaussian(9,3),
                        numpy.fromfile(os.path.join(currentdir,"convolve_02"),sep=' ').reshape(1000,1000)
                    ),
                    (
                        numpy.fromfile(os.path.join(currentdir,"convolve_01"),sep=' ').reshape(1000,1000),
                        numpy.fromfile(os.path.join(currentdir,"convolve_04"),sep=' ').reshape(1000,1000),
                        numpy.fromfile(os.path.join(currentdir,"convolve_03"),sep=' ').reshape(1000,1000)
                    )
                ]

@pytest.mark.parametrize(("input","kernel","result"), convolve_data)
def test_convolve(input,kernel,result):
    input_image = AstroImage(data=input)
    kernel_image = AstroImage(data=kernel)
    input_image.convolve(kernel_image)
    verifyData(input_image.data,result)

rotate_data =   [
                            (
                                numpy.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
                                39.8,
                                numpy.fromfile(os.path.join(currentdir,"rotate_02"),sep=' ').reshape(1000,1000)
                            ),
                            (
                                numpy.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
                                90.,
                                numpy.fromfile(os.path.join(currentdir,"rotate_03"),sep=' ').reshape(1000,1000)
                            ),
                            (
                                numpy.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
                                -28.9,
                                numpy.fromfile(os.path.join(currentdir,"rotate_04"),sep=' ').reshape(1000,1000)
                            ),
                            (
                                numpy.fromfile(os.path.join(currentdir,"rotate_01"),sep=' ').reshape(1000,1000),
                                180.,
                                numpy.fromfile(os.path.join(currentdir,"rotate_05"),sep=' ').reshape(1000,1000)
                            )
                ]

@pytest.mark.parametrize(("input","angle","output"), rotate_data)
def test_rotate(input,angle,output):
    image = AstroImage(data=input)
    image.rotate(angle)
    verifyData(image.data,output)

addWithOffset_data =   [
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                numpy.full((10,10),3.),
                                (4,-2),
                                {},
                                numpy.fromfile(os.path.join(currentdir,"addWithOffset_01"),sep=' ').reshape(10,10)
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                numpy.full((10,10),3.),
                                (0,0),
                                {'scale':[2.,2.]},
                                numpy.fromfile(os.path.join(currentdir,"addWithOffset_02"),sep=' ').reshape(10,10)
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                numpy.full((10,10),3.),
                                (0,0),
                                {'scale':[0.5,0.5]},
                                numpy.fromfile(os.path.join(currentdir,"addWithOffset_03"),sep=' ').reshape(10,10)
                            )
                        ]

@pytest.mark.parametrize(("in1","in2","offset","inputs","result"), addWithOffset_data)
def test_addWithOffset(in1,in2,offset,inputs,result):
    im1 = AstroImage(data=in1,**inputs)
    im2 = AstroImage(data=in2)
    im1.addWithOffset(im2,offset[0],offset[1])
    verifyData(im1.data,result)

addWithAlignment_data =   [
                            (
                                numpy.zeros((10,10),dtype='float32'),
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                {},
                                numpy.arange(100,dtype='float32').reshape(10,10)
                            ),
                            (
                                numpy.zeros((10,10),dtype='float32'),
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                {'dec':(5./3600.)},
                                numpy.fromfile(os.path.join(currentdir,"addWithAlignment_01"),sep=' ').reshape(10,10)
                            ),
                            (
                                numpy.zeros((10,10),dtype='float32'),
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                {'pa':90.},
                                numpy.fromfile(os.path.join(currentdir,"addWithAlignment_02"),sep=' ').reshape(10,10)
                            ),
                            (
                                numpy.zeros((10,10),dtype='float32'),
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                {'pa':38.23},
                                numpy.fromfile(os.path.join(currentdir,"addWithAlignment_03"),sep=' ').reshape(10,10)
                            ),
                            (
                                numpy.zeros((10,10),dtype='float32'),
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                {'pa':27.5,'ra':(-5./3600.)},
                                numpy.fromfile(os.path.join(currentdir,"addWithAlignment_04"),sep=' ').reshape(10,10)
                            )
                        ]

@pytest.mark.parametrize(("in1","in2","inputs","result"), addWithAlignment_data)
def test_addWithAlignment(in1,in2,inputs,result):
    im1 = AstroImage(data=in1)
    im2 = AstroImage(data=in2,**inputs)
    im1.addWithAlignment(im2)
    verifyData(im1.data,result)

rescale_data =   [
                            (
                                numpy.full((10,10),3.),
                                numpy.array([0.5,0.5]),
                                numpy.full((20,20),.75)
                            ),
                            (
                                numpy.full((10,10),3.),
                                numpy.array([2.,2.]),
                                numpy.full((5,5),12.)
                            ),
                            (
                                numpy.full((10,10),3.),
                                numpy.array([1.,2.]),
                                numpy.full((10,5),6.)
                            ),
                            (
                                numpy.full((10,10),3.),
                                numpy.array([1.25,.75]),
                                numpy.full((8,13),2.8846)
                            )
                        ]

@pytest.mark.parametrize(("input","scale","result"), rescale_data)
def test_rescale(input,scale,result):
    im1 = AstroImage(data=input)
    im1.rescale(scale)
    verifyData(im1.data,result)
    im1.rescale([1.,1.])
    verifyData(im1.data,input)

crop_data =   [
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [5.,5.], [0,0],
                                numpy.arange(100,dtype='float32').reshape(10,10)[3:8,3:8]
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [3.,3.], [4,4],
                                numpy.arange(100,dtype='float32').reshape(10,10)[:3,:3]
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [12.,12.], [0,0],
                                numpy.fromfile(os.path.join(currentdir,"crop_01"),sep=' ').reshape(12,12)
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [12.,12.], [4,4],
                                numpy.fromfile(os.path.join(currentdir,"crop_02"),sep=' ').reshape(12,12)
                            )
                        ]

@pytest.mark.parametrize(("input","size","offset","result"), crop_data)
def test_crop(input,size,offset,result):
    im1 = AstroImage(data=input)
    im1.crop(size[0],size[1],offset[0],offset[1])
    verifyData(im1.data,result)

bin_data =   [
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [2,2],
                                numpy.fromfile(os.path.join(currentdir,"bin_01"),sep=' ').reshape(5,5)
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [3,3],
                                numpy.fromfile(os.path.join(currentdir,"bin_02"),sep=' ').reshape(3,3)
                            ),
                            (
                                numpy.arange(100,dtype='float32').reshape(10,10),
                                [2,3],
                                numpy.fromfile(os.path.join(currentdir,"bin_03"),sep=' ').reshape(3,5)
                            ),
                            (
                                numpy.arange(90,dtype='float32').reshape(10,9),
                                [3,2],
                                numpy.fromfile(os.path.join(currentdir,"bin_04"),sep=' ').reshape(5,3)
                            )
                        ]

@pytest.mark.parametrize(("input","bin","result"), bin_data)
def test_bin(input,bin,result):
    im1 = AstroImage(data=input)
    im1.bin(bin[0],bin[1])
    verifyData(im1.data,result)

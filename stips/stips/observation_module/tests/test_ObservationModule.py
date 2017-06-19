from __future__ import division

import cPickle,inspect,numpy,os,pytest,sys
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
topdir = os.path.dirname(parentdir)
sys.path.insert(0,parentdir)

os.environ['TERM'] = 'xterm'
os.environ['PYSYN_CDBS'] = os.path.join(os.getcwd(),'sim_input','cdbs')
os.environ['WEBBPSF_PATH'] = os.path.join(os.getcwd(),"sim_input","webbpsf-data")
os.environ['WEBBPSF_SKIP_CHECK'] = '1'

import astropy.io.fits as pyfits
from astropy.table import Table,Column
import test_AstroImage

from AstroImage import AstroImage
from ObservationModule import ObservationModule

def verifyParameters(obs,results):
    assert obs.instrument_name == results['instrument_name']
    assert obs.filters == results['filters']
    assert obs.offsets == results['offsets']
    assert obs.coadd == results['coadd']
    assert obs.exptime == results['exptime']
    assert obs.oversample == results['oversample']
    assert obs.psf_commands == results['psf_commands']
    assert obs.ra == results['ra']
    assert obs.dec == results['dec']
    assert obs.pa == results['pa']
    assert obs.seed == results['seed']
    assert obs.id == results['id']
    assert obs.prefix == results['prefix']
    assert obs.in_path == results['in_path']
    assert obs.cat_path == results['out_path']
    assert obs.out_path == results['out_path']
    assert obs.poisson == results['poisson']
    assert obs.flat == results['flat']
    assert obs.dark == results['dark']
    assert obs.cosmic == results['cosmic']
    assert obs.version == results['version']

creation_data =         [
                            (
                                {},
                                {},
                                {'instrument_name':'NIRCamShort','filters':[],'offsets':[],'coadd':1,
                                 'exptime':1.,'oversample':1,'psf_commands':"",'ra':0.,'dec':0.,
                                 'pa':0.,'seed':0,'id':0,'prefix':'sim','in_path':os.getcwd(),
                                 'cat_path':os.getcwd(),'out_path':os.getcwd(),'poisson':True,
                                 'flat':True,'dark':True,'cosmic':True,'version':"0.0"}
                            ),
                            (
                                {'instrument':'MIRI','filters':['F1130W'],'offset':[{'offset_ra':1.3,
                                 'offset_dec':-3.9,'offset_pa':38.3}],'coadd':32,'exptime':500.,
                                 'oversample':12,'pupil_mask':'Commanding','id':7},
                                {'out_prefix':'testing_alpha','in_path':'/path/to/somewhere',
                                 'out_path':'/path/to/somewhere/else','general':{'ra':7.3,'dec':-9.1,
                                 'pa':-12.5,'seed':3210},'residuals':{'do_err_flat':False,
                                 'do_err_dark':False,'do_err_cray':False},'version':'9.5'},
                                {'instrument_name':'MIRI','filters':['F1130W'],
                                 'offsets':[{'offset_ra':1.3,'offset_dec':-3.9,'offset_pa':38.3}],
                                 'coadd':32,'exptime':16000.,'oversample':12,
                                 'psf_commands':"Commanding",'ra':7.3,'dec':-9.1,'pa':-12.5,
                                 'seed':3210,'id':7,'prefix':'testing_alpha',
                                 'in_path':'/path/to/somewhere','cat_path':'/path/to/somewhere/else',
                                 'out_path':'/path/to/somewhere/else','poisson':True,'flat':False,
                                 'dark':False,'cosmic':False,'version':"9.5"}
                            )
                        ]

@pytest.mark.parametrize(("obs","inputs","results"), creation_data)
def test_creation(obs,inputs,results):
    obs = ObservationModule(obs,**inputs)
    verifyParameters(obs,results)

prep_data = [
                (
                    {'input_image':os.path.join(currentdir,"iabf01c5q_flt.fits"),'image_scale':0.3,
                     'image_wcs':True,'image_ext':1},
                    AstroImage.initFromFits(os.path.join(currentdir,"iabf01c5q_flt.fits"),ext=1)
                ),
                (
                    {'input_image':os.path.join(currentdir,"iabf01c5q_flt.fits"),'image_scale':1.2,
                     'image_wcs':False,'image_ext':1},
                    AstroImage.initDataFromFits(os.path.join(currentdir,"iabf01c5q_flt.fits"),ext=1,
                                                             scale=[1.2,1.2])
                )
            ]

@pytest.mark.parametrize(("in_img","astro_img"), prep_data)
def test_prepImage(in_img,astro_img):
    obs = ObservationModule({})
    obs.prepImage(in_img)
    test_AstroImage.verifyImage(obs.images[in_img['input_image']],astro_img)
    test_AstroImage.verifyData(obs.images[in_img['input_image']].data,astro_img.data)

obs_data =  [
                (
                    #observation
                    {'filters':["F070W","F115W","F150W"],
                     'offset':[{'offset_ra':0.,'offset_dec':0.,'offset_pa':0.,'offset_centre':False},
                               {'offset_ra':33.7,'offset_dec':-12.1,'offset_pa':39.1,'offset_centre':False}]},
                    #kwargs
                    {},
                    #expected results per-observation
                    [("F070W",0.,0.,0.),("F070W",0.0093611,-0.0033611,39.1),
                     ("F115W",0.,0.,0.),("F115W",0.0093611,-0.0033611,39.1),
                     ("F150W",0.,0.,0.),("F150W",0.0093611,-0.0033611,39.1),None]
                ),
                (
                    #observation
                    {'filters':['F140M','F162M'],
                     'offset':[{'offset_ra':0.,'offset_dec':0.,'offset_pa':0.,'offset_centre':False},
                               {'offset_ra':120.,'offset_dec':-180.,'offset_pa':39.1,'offset_centre':False}]},
                    #kwargs
                    {'general':{'ra':39.,'dec':-15.,'pa':45.}},
                    #expected results per-observation
                    [("F140M",39.,-15.,45.),("F140M",39.0333333,-15.05,84.1),
                     ("F162M",39.,-15.,45.),("F162M",39.0333333,-15.05,84.1),None]
                ),
                (
                    #observation
                    {'filters':["F070W","F115W"],
                     'offset':[{'offset_ra':0.,'offset_dec':0.,'offset_pa':0.,'offset_centre':True},
                               {'offset_ra':360.,'offset_dec':-360.,'offset_pa':45.,'offset_centre':True}]},
                    #kwargs
                    {},
                    #expected results per-observation
                    [("F070W",359.99984166666667,-0.135725,0.),("F070W",0.09984167,-0.235725,45.),
                     ("F115W",359.99984166666667,-0.135725,0.),("F115W",0.09984167,-0.235725,45.),None]
                )
            ]

@pytest.mark.parametrize(("obs","kwargs","results"), obs_data)
def test_nextObservation(obs,kwargs,results):
    om = ObservationModule(obs,**kwargs)
    for res in results:
        obs_num = om.nextObservation()
        if res is None:
            assert obs_num is None
        else:
            filter,ra,dec,pa = res
            assert filter == om.instrument.filter
            numpy.testing.assert_allclose((ra),(om.instrument.ra),atol=1e-3)
            numpy.testing.assert_allclose((dec),(om.instrument.dec),atol=1e-3)
            numpy.testing.assert_allclose((pa),(om.instrument.pa),atol=1e-3)

#     -----------
#     def addCatalogue(self, catalogue):
#         """
#         Add a catalogue to the internal image.
#         
#         Parameters
#         ----------
#         
#         self: obj
#             Class instance.
#         
#         catalogue: string
#             Name of catalogue file
#         """
#         self.logger.info("Running catalogue %s",catalogue)
#         cats = self.instrument.addCatalogue(catalogue)
#         self.logger.info('Finished catalogue %s',catalogue)
#         return cats
# 
#     -----------
#     def addImage(self, img, units):
#         """
#         Create an internal image from a user-provided file
#         
#         Parameters
#         ----------
#         
#         self: obj
#             Class instance.
#         
#         file: string
#             File to read
#         
#         units: string
#             Unit type of file
#         """
#         self.logger.info('Adding image %s to observation',file)
#         self.instrument.addImage(self.images[img],units)
#         self.logger.info('Image Added')
# 
#     -----------
#     def addError(self):
#         """
#         Add internal sources of error to the image
#         
#         Parameters
#         ----------
#         
#         self: obj
#             Class instance.
#         """
#         psf_name = "%s_%d_psf.fits" % (self.imgbase,self.obs_count)
#         self.instrument.psf.toFits(psf_name)
#         self.logger.info("Adding Error")
#         readnoise is always true
#         self.instrument.addError(self.poisson,True,self.flat,self.dark,self.cosmic,self.exptime,self.coadd)
#         self.logger.info("Finished Adding Error")
#         return psf_name
#         
#     -----------
#     def finalize(self):
#         """
#         Finalize FITS file
#         
#         Parameters
#         ----------
#         
#         self: obj
#             Class instance.
#         """
#         self.instrument.toFits("%s_%d.fits" % (self.imgbase,self.obs_count))
#         mosaics = self.instrument.toMosaic("%s_%d_mosaic.fits"%(self.imgbase,self.obs_count))
#         return "%s_%d.fits"%(self.imgbase,self.obs_count),mosaics,self.params
#     
#     -----------
#     def totalObservations(self):
#         """
#         Return the total number of observations
#         """
#         return len(self.observations)

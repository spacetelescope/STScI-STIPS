from __future__ import absolute_import, division
__filetype__ = "detector"

#External Modules
import os

#Local Modules
from ..astro_image import AstroImage
from .instrument import Instrument
from .hst_instrument import HstInstrument
from ..utilities import OffsetPosition, GetStipsData

class WFC3IR(HstInstrument):
    __classtype__ = "detector"
    """
    The NircamShort class contains the necessary constants and modifications to run NircamShort
    observations.
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        instrument: string, which instrument this is
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Still only a base class. Init is super init.
        """
        #Initialize superclass
        super(WFC3IR,self).__init__(**kwargs)
        
        #For WFC3 at the moment, there is no way to oversample the PSF. Thus no oversample.
        # self.oversample = kwargs.get('oversample', 1)


    def resetPSF(self):
        if self.filter not in self.FILTERS:
            raise ValueError("Filter %s is not a valid WFC3IR filter" % (self.filter))
        psf_path = GetStipsData(os.path.join('psf_data', 'PSF_WFC3IR_{}.fits'.format(self.filter)))
        self.psf = AstroImage.initDataFromFits(psf_path, detname="WFC3IRPSF", logger=self.logger)


    @property
    def bandpass(self):
        import pysynphot as ps
        ps.setref(**self.REFS)
        obsmode = "wfc3,ir," + self.filter.lower()
        return ps.ObsBandpass(obsmode)
            
    def generateReadnoise(self):
        """
        Readnoise formula that is similar to JWST ETC.
        """
        return 12. #max ramp, lowest noise
    
    @classmethod
    def handleDithers(cls,form):
        """
        WFC3-specific dither patterns
        """
        dither_pattern = form['dither_type']
        if dither_pattern == "BOX-UVIS":
            initial_dithers = cls.DITHER_OFFSETS[dither_pattern]
        elif dither_pattern == "BLOB":
            dither_points = int(form['dither_points'])
            initial_dithers = []
            while len(initial_dithers) <= dither_points:
                initial_dithers += cls.DITHER_OFFSETS[dither_pattern]
            initial_dithers = initial_dithers[:dither_points]
        else: #dither_pattern == "NO" and instrument potentially matters again
            initial_dithers = [(0.,0.)]
        dither_subpixel = form['dither_subpixel']
        subpixel_dithers = cls.DITHER_OFFSETS['SUBPIXEL'][dither_subpixel]
        return cls.doSubpixel(initial_dithers,subpixel_dithers)
            
    INSTRUMENT = "WFC3IR"
    DETECTOR = "WFC3IR"
    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    # Instrument Offset is from JWST V2,V3 centre
    INSTRUMENT_OFFSET = (0.,0.,0.) #*****at some point get the exact SIAF for HST
    # Detector Offsets are RELATIVE TO INSTRUMENT OFFSET
    DETECTOR_OFFSETS = ((0.,0.,0.),) #*****at some point get the exact SIAF for HST
    OFFSET_NAMES = (("WFC3IR"),)
    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
    N_DETECTORS = [1]
    DETECTOR_SIZE = (1014,1014) #pixels
    PIXEL_SIZE = 18.0 #um
    SCALE = [0.1355,0.1211]
    FILTERS = ('F110W','F160W')
    DEFAULT_FILTER = 'F160W'
    FLATFILE = 'err_flat_wfc3ir.fits'
    DARKFILE = 'err_rdrk_wfc3ir.fits' # IREF, IHB
    BACKGROUND = {  'none': {'F110W': 0., 'F160W': 0.},
                    'avg':  {'F110W': 1.401E+00, 'F160W': 7.521E-01}
                 }
    BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
    BACKGROUNDS = ['None', 'Average zodiacal background', 'Median zodiacal background', 'Maximum zodiacal background', 'Minimum zodiacal background']
    BGTEXT = {'none': 'None', 'avg': 'Average zodiacal background', 
              'med': 'Median zodiacal background', 'max': 'Maximum zodiacal background', 
              'min': 'Minimum zodiacal background', 'custom': 'Custom thermal background rate'}
    #PHOTFNU has units of Jy
    PHOTFNU = {'F110W':6.760E-08, 'F160W':1.505E-07}
    #PHOTPLAM has units of um
    PHOTPLAM = {'F110W':1.1534, 'F160W':1.5369}
    DITHERS = ("SUBPIXEL ONLY","BOX-UVIS","BLOB")
    DITHER_POINTS = {
                        "SUBPIXEL ONLY": ["0"],
                        "BOX-UVIS": ["4"],
                        "BLOB": ["1", "2"]
                     }
    DITHER_SIZE = {
                    "SUBPIXEL ONLY": ["STANDARD"],
                    "BOX-UVIS": ["Standard"],
                    "BLOB": ["Standard"]
                  }
    DITHER_SUBPIXEL = {
                        "SUBPIXEL ONLY": ["LINE","LINE-3PT","BOX-MIN"],
                        "BOX-UVIS": ["NONE","LINE","LINE-3PT","BOX-MIN"],
                        "BLOB": ["NONE","LINE","LINE-3PT","BOX-MIN"]
                      }
    DITHER_OFFSETS = {
                        "BOX-UVIS":[(-11.071,-17.744),(11.947,-17.457),(11.071,17.744),(-11.947,17.457)],
                        "BLOB":[(-1.930,-1.729),(1.930,1.729)],
                        "SUBPIXEL": {
                                        "NONE":     [(0.000,0.000)],
                                        "BOX-MIN":  [(0.000,0.000),(0.542,0.182),(0.339,0.485),(-0.203,0.303)],
                                        "LINE":     [(0.000,0.000),(0.474,0.424)],
                                        "LINE-3PT": [(0.000,0.000),(0.451,0.403),(0.902,0.806)]
                                    }
                     }

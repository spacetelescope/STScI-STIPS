from __future__ import absolute_import, division
__filetype__ = "detector"

#External Modules
import os

import  numpy as np

from astropy import wcs
from astropy.io import fits as pyfits

#Local Modules
from ..astro_image import AstroImage
from .instrument import Instrument
from .wfirst_instrument import WfirstInstrument
from ..utilities import OffsetPosition

from .. import __version__ as stips_version

class WFI(WfirstInstrument):
    __classtype__ = "detector"
    """
    The WFI class contains the necessary constants and modifications to run WFI
    observations.
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        instrument: string, which instrument this is
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Initializes detectors (single detector for now).
        """
        self.classname = self.__class__.__name__
        #Initialize superclass
        super(WFI, self).__init__(**kwargs)
        
        #Set oversampling
        self.oversample = kwargs.get('oversample', self.OVERSAMPLE_DEFAULT)
        
        #Set PSF grid points
        self.grid_size = kwargs.get('grid_size', self.PSF_GRID_SIZE_DEFAULT)


    def generateReadnoise(self):
        """
        Readnoise formula that is similar to HST ETC.
        """
        return 12. #max ramp, lowest noise
    
    @classmethod
    def handleDithers(cls,form):
        """
        For now, handle WFI as equivalent to WFC3IR in dither patterns.
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
        else: #dither_pattern == "SUBPIXEL ONLY" and instrument potentially matters again
            initial_dithers = [(0.,0.)]
        dither_subpixel = form['dither_subpixel']
        subpixel_dithers = cls.DITHER_OFFSETS['SUBPIXEL'][dither_subpixel]
        return cls.doSubpixel(initial_dithers,subpixel_dithers)
    
    INSTRUMENT = "WFI"
    DETECTOR = "WFI"

    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    DETECTOR_OFFSETS = ((0.,0.,0.), (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), 
                        (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), 
                        (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), 
                        (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), (0.,0.,0.), 
                        (0.,0.,0.), (0.,0.,0.))
    OFFSET_NAMES = ("SCA01", "SCA02", "SCA03", "SCA04", "SCA05", "SCA06", 
                    "SCA07", "SCA08", "SCA09", "SCA10", "SCA11", "SCA12", 
                    "SCA13", "SCA14", "SCA15", "SCA16", "SCA17", "SCA18")
    N_OFFSET = {1: (0., 0., 0.), 2: (0., 0., 0.), 8: (0., 0., 0.),
                4: (0., 0., 0.), 5: (0., 0., 0.), 6: (0., 0., 0.),
                7: (0., 0., 0.), 8: (0., 0., 0.), 9: (0., 0., 0.), 
                10: (0., 0., 0.), 11: (0., 0., 0.), 12: (0., 0., 0.), 
                13: (0., 0., 0.), 14: (0., 0., 0.), 15: (0., 0., 0.), 
                16: (0., 0., 0.), 17: (0., 0., 0.), 18: (0., 0., 0.)}

    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use    
    N_DETECTORS = [1]
    INSTRUMENT_OFFSET = (0.,0.,0.) #Presumably there is one, but not determined
    DETECTOR_SIZE = (4096,4096) #pixels
    PIXEL_SIZE = 18.0 #um (Assume for now)
    SCALE = [0.11,0.11] #Assume for now
    FILTERS = ('F062', 'F087', 'F106', 'F129', 'F158', 'F184', 'F146', 'F149') #W149 needs to go away at some point.
    DEFAULT_FILTER = 'F184' #Assume for now
    
    PSF_INSTRUMENT = "WFI"
    
    # Reference Files
    FLATFILE = 'err_flat_wfi.fits' #Use for the moment
    DARKFILE = 'err_rdrk_wfi.fits' # IREF, IHB (use for the moment)

    # Background Values
    BACKGROUND = {  'none': {'F062': 0., 'F087': 0.,'F106': 0.,'F129': 0.,'F158': 0.,'F184': 0., 'F146': 0., 'F149': 0.},
                    'avg':  {'F062': 1.401E+00, 'F087': 1.401E+00, 'F106': 1.401E+00, 'F129': 7.000E-01,
                             'F158': 7.521E-01, 'F184': 8.500E-01, 'F146': 7.000E-01, 'F149': 7.000E-01}
                 }
    BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
    BACKGROUNDS = ['None', 'Average zodiacal background', 'Median zodiacal background', 'Maximum zodiacal background', 'Minimum zodiacal background']
    BGTEXT = {'none': 'None', 'avg': 'Average zodiacal background', 
              'med': 'Median zodiacal background', 'max': 'Maximum zodiacal background', 
              'min': 'Minimum zodiacal background', 'custom': 'Custom thermal background rate'}
    #PHOTFNU has units of Jy
    #For now, just assuming similar PHOTFNU to WFC3IR.
    PHOTFNU = {     'F062': 1.0e-8, 'F087':1.0e-8, 'F106':1.0e-8, 'F129':1.0e-8, 'F158':1.0e-8, 'F184':1.0e-8, 'F146':1.0e-8, 'F149':1.0e-8}
    #PHOTPLAM has units of um
    #For now, just put them in the middle
    PHOTPLAM = {'F062': 0.6700, 'F087':0.8735, 'F106':1.0595, 'F129':1.2925, 'F158':1.577, 'F184':1.5815, 'F146':1.4635, 'F149':1.4635}
    #For now, just put in HST-style dithers.
    DITHERS = ("SUBPIXEL ONLY","BOX-UVIS","BLOB") #Assume for now
    DITHER_POINTS = {
                        "SUBPIXEL ONLY": ["0"],
                        "BOX-UVIS": ["4"],
                        "BLOB": ["1", "2"]
                     }
    DITHER_SIZE = {
                    "SUBPIXEL ONLY": ["STABDARD"],
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

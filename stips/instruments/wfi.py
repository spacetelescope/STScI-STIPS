from __future__ import absolute_import, division
__filetype__ = "detector"

#External Modules
import os

import  numpy as np

from astropy.io import fits as pyfits

#Local Modules
from ..astro_image import AstroImage
from .instrument import Instrument
from .wfirst_instrument import WfirstInstrument
from ..utilities import OffsetPosition

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
        super(WFI,self).__init__(**kwargs)
        
        #Set oversampling
        self.oversample = kwargs.get('oversample', 1)

    def resetPSF(self):
        import webbpsf
        if self.filter not in self.FILTERS:
            raise ValueError("Filter %s is not a valid WFI filter" % (self.filter))
        have_psf = False
        if os.path.exists(os.path.join(self.out_path, "psf_cache")):
            if os.path.exists(os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("WFI", self.filter, self.oversample))):
                with pyfits.open(os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("WFI", self.filter, self.oversample))) as psf:
                    if psf[0].header['VERSION'] >= webbpsf.__version__ and (self.psf_commands is None or self.psf_commands == ''):
                        self.psf = AstroImage(data=psf[0].data, detname="WFI {} PSF".format(self.filter), logger=self.logger)
                        have_psf = True
        if not have_psf:
            base_state = self.getState()
            self.updateState(base_state+"<br /><span class='indented'>Generating PSF</span>")
            from webbpsf import wfirst
            ins = wfirst.WFI()
            if self.psf_commands is not None and self.psf_commands != '':
                for attribute,value in self.psf_commands.iteritems():
                    setattr(ins,attribute,value)
            ins.filter = self.filter
            max_safe_size = int(np.floor(30. * self.PHOTPLAM[self.filter] / (2. * self.SCALE[0])))
            max_ins_size = max(self.DETECTOR_SIZE) * self.oversample
            max_conv_size = int(np.floor(2048 / self.oversample))
            self._log("info", "PSF choosing between {}, {}, and {}".format(max_safe_size, max_ins_size, max_conv_size))
            if hasattr(ins, 'calc_psf'):
                ins.calcPSF = ins.calc_psf
            psf = ins.calcPSF(oversample=self.oversample, fov_pixels=min(max_safe_size, max_ins_size, max_conv_size), normalize='last')
            self._log("info", "PSF Total Flux: {}".format(np.sum(psf[0].data)))
            psf[0].header['VERSION'] = webbpsf.__version__
            if os.path.exists(os.path.join(self.out_path, "psf_cache")):
                dest = os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("WFI", self.filter, self.oversample))
                pyfits.writeto(dest, psf[0].data, header=psf[0].header, overwrite=True)
            self.psf = AstroImage(data=psf[0].data, detname="WFI %s PSF" % (self.filter), logger=self.logger)
            self.updateState(base_state)
        
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
    DETECTOR_OFFSETS = ((0.,0.,0.),) #There will be 18, but simulate as single
    OFFSET_NAMES = (("WFIRST-WFI"),)
    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
    N_DETECTORS = [1]
    INSTRUMENT_OFFSET = (0.,0.,0.) #Presumably there is one, but not determined
    DETECTOR_SIZE = (4096,4096) #pixels
    PIXEL_SIZE = 18.0 #um (Assume for now)
    SCALE = [0.11,0.11] #Assume for now
    DIST_A =   [[   0.,             0.,             0.],
                [   0.,             0.,             0.],
                [   0.,             0.,             0.]]
    DIST_B =   [[   0.,             0.,             0.],
                [   0.,             0.,             0.],
                [   0.,             0.,             0.]]
    DIST_AP =  [[   0.,             0.,             0.],
                [   0.,             0.,             0.],
                [   0.,             0.,             0.]]
    DIST_BP =  [[   0.,             0.,             0.],
                [   0.,             0.,             0.],
                [   0.,             0.,             0.]]
    FILTERS = ('R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'W149') #W149 needs to go away at some point.
    DEFAULT_FILTER = 'F184' #Assume for now
    FLATFILE = 'err_flat_wfi.fits' #Use for the moment
    DARKFILE = 'err_rdrk_wfi.fits' # IREF, IHB (use for the moment)
    BACKGROUND = {  'none': {'R062': 0., 'Z087': 0.,'Y106': 0.,'J129': 0.,'H158': 0.,'F184': 0., 'W146': 0., 'W149': 0.},
                    'avg':  {'R062': 1.401E+00, 'Z087': 1.401E+00, 'Y106': 1.401E+00, 'J129': 7.000E-01,
                             'H158': 7.521E-01, 'F184': 8.500E-01, 'W146': 7.000E-01, 'W149': 7.000E-01}
                 }
    BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
    BACKGROUNDS = ['None', 'Average zodiacal background', 'Median zodiacal background', 'Maximum zodiacal background', 'Minimum zodiacal background']
    BGTEXT = {'none': 'None', 'avg': 'Average zodiacal background', 
              'med': 'Median zodiacal background', 'max': 'Maximum zodiacal background', 
              'min': 'Minimum zodiacal background', 'custom': 'Custom thermal background rate'}
    #PHOTFNU has units of Jy
    #For now, just assuming similar PHOTFNU to WFC3IR.
    PHOTFNU = {     'R062': 1.0e-8, 'Z087':1.0e-8, 'Y106':1.0e-8, 'J129':1.0e-8, 'H158':1.0e-8, 'F184':1.0e-8, 'W146':1.0e-8, 'W149':1.0e-8}
    #PHOTPLAM has units of um
    #For now, just put them in the middle
    PHOTPLAM = {'R062': 0.6700, 'Z087':0.8735, 'Y106':1.0595, 'J129':1.2925, 'H158':1.577, 'F184':1.5815, 'W146':1.4635, 'W149':1.4635}
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

from __future__ import absolute_import, division
__filetype__ = "detector"

#External Modules
import os

import numpy as np

from astropy.io import fits as pyfits

#Local Modules
from ..astro_image import AstroImage
from .instrument import Instrument
from .jwst_instrument import JwstInstrument
from ..utilities import OffsetPosition

class NIRCamBase(JwstInstrument):
    """
    Base class with data common to NIRCam in all instances. Uses NIRCamShort dither patterns.
    """

    def __init__(self, **kwargs):
        """
        Init does the following:
            - super().__init__()
            - looks for self.oversample (if present)
            - finds filter and verifies that it is a valid filter
            - finds target co-ordinates and PA (if present)
            - looks for input and output paths
            - sets flat residual reference file
            - sets dark residual reference file
            - determines instrument/filter zeropoint
            - creates PSF
            - creates detectors and specifies their relative offset
        """
        self.classname = self.__class__.__name__
        #Initialize superclass
        super(NIRCamBase,self).__init__(**kwargs)

        #Set oversampling
        self.oversample = kwargs.get('oversample', 1)

        #Adjust # of detectors based on keyword:
        n_detectors = int(kwargs.get('detectors', len(self.DETECTOR_OFFSETS)))
        self.DETECTOR_OFFSETS = self.DETECTOR_OFFSETS[:n_detectors]
        self.OFFSET_NAMES = self.OFFSET_NAMES[:n_detectors]
        self.CENTRAL_OFFSET = self.N_OFFSET[n_detectors]
        self._log('info', "{} with {} detectors. Central offset {}".format(self.DETECTOR, n_detectors, self.CENTRAL_OFFSET))

    def resetPSF(self):
        import webbpsf
        if self.filter not in self.FILTERS:
            raise ValueError("Filter %s is not a valid %s filter" % (self.filter,self.classname))
        have_psf = False
        if os.path.exists(os.path.join(self.out_path, "psf_cache")):
            if os.path.exists(os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("NIRCam", self.filter, self.oversample))):
                with pyfits.open(os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("NIRCam", self.filter, self.oversample))) as psf:
                    if psf[0].header['VERSION'] >= webbpsf.__version__ and (self.psf_commands is None or self.psf_commands == ''):
                        self.psf = AstroImage(data=psf[0].data, detname="NIRCam {} PSF".format(self.filter), logger=self.logger)
                        have_psf = True
        if not have_psf:
            base_state = self.getState()
            self.updateState(base_state+"<br /><span class='indented'>Generating PSF</span>")
            ins = webbpsf.NIRCam()
            if self.psf_commands is not None and self.psf_commands != '':
                for attribute,value in self.psf_commands.iteritems():
                    setattr(ins,attribute,value)
            ins.filter = self.filter
            max_safe_size = int(np.floor(30. * self.PHOTPLAM[self.filter] / (2. * self.SCALE[0])))
            max_ins_size = max(self.DETECTOR_SIZE) * self.oversample
            max_conv_size = int(np.floor(2048 / self.oversample))
            self._log("info", "PSF choosing between {}, {}, and {}".format(max_safe_size, max_ins_size, max_conv_size))
            psf = ins.calcPSF(oversample=self.oversample,fov_pixels=min(max_safe_size, max_ins_size, max_conv_size))
            psf[0].header['VERSION'] = webbpsf.__version__
            if os.path.exists(os.path.join(self.out_path, "psf_cache")):
                dest = os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("NIRCam", self.filter, self.oversample))
                pyfits.writeto(dest, psf[0].data, header=psf[0].header, overwrite=True)
            self.psf = AstroImage(data=psf[0].data,detname="NIRCam %s PSF" % (self.filter),logger=self.logger)
            self.updateState(base_state)
    
    def generateReadnoise(self):
        """
        Readnoise formula that is similar to JWST ETC.
        """
        k = 55.07
        a = -0.26

        if self.exptime > 1000:
            numGroups = np.ceil(self.exptime/1000.0)
            timePerGroup = self.exptime/numGroups
        else:
            numGroups = 1
            timePerGroup = self.exptime

        rdns = np.sqrt(numGroups) * k * timePerGroup**(a)
        return rdns
    
    @classmethod
    def handleDithers(cls,form):
        """
        Handles a dither pattern request with NIRCam (in some form) as an instrument.
        """
        dither_pattern = form['dither_type']
        dither_subpixel = int(form['dither_subpixel'])
        if dither_subpixel <= 8:
            subpixel_dithers = cls.DITHER_OFFSETS['SUBPIXEL'][dither_subpixel-1]
        else:
            subpixel_dithers = cls.DITHER_OFFSETS['SUBPIXEL'][8][:dither_subpixel]
        if dither_pattern == "FULL":
            dither_points = form['dither_points']
            initial_dithers = cls.DITHER_OFFSETS[dither_pattern][dither_points]
        elif dither_pattern == "INTRAMODULE": #NIRCam
            dither_points = int(form['dither_points'])
            initial_dithers = cls.DITHER_OFFSETS[dither_pattern][:dither_points]
        elif dither_pattern == "INTRASCA": #NIRCam
            dither_points = int(form['dither_points'])
            dither_size = form['dither_size']
            initial_dithers = cls.DITHER_OFFSETS[dither_pattern][dither_size][:dither_points]
        else: #no pattern, only sub-pixel dithers
            initial_dithers = [(0.,0.)]
        return cls.doSubpixel(initial_dithers,subpixel_dithers)

    INSTRUMENT = "NIRCam"
    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    DETECTOR_SIZE = (2040,2040) #pixels
    PIXEL_SIZE = 18.0 #um
    DISTORTION = {
                    'A1': { 'DIST_A': [[0.,             0.,             0.,             0.,             0.,             0.],
                                       [3.11450E-2,     0.,             0.,             0.,             0.,             0.],
                                       [3.19296E-9,     -2.09854E-7,    -3.45451E-8,    0.,             0.,             0.],
                                       [9.98433E-12,    2.07312E-12,    1.16144E-11,    1.02269E-12,    0.,             0.],
                                       [-7.28858E-16,   -9.01713E-16,   -9.36343E-16,   -8.89534E-16,   -1.69081E-16,   0.],
                                       [1.40202E-19,    1.07762E-20,    2.87120E-19,    9.20540E-21,    1.43515E-19,    -1.75089E-22]],
                            'DIST_B': [[0.,             0.,             0.,             0.,             0.,             0.],
                                       [5.35493E-5,     3.13335E-2,     0.,             0.,             0.,             0.],
                                       [6.52433E-8,     3.91325E-8,     -1.46390E-7,    0.,             0.,             0.],
                                       [6.20123E-13,    9.58663E-12,    1.63279E-12,    1.13863E-11,    0.,             0.],
                                       [-1.42E-16,      -6.51E-16,      -1.15E-15,      -6.25E-16,      -9.33E-16,      0.],
                                       [4.18516E-21,    1.55935E-19,    1.75359E-20,    2.92096E-19,    1.22758E-20,    1.34127E-19]]},
                    'A2': { 'DIST_A': [[0.,             0.,             0.,             0.,             0.,             0.],
                                       [3.07617E-2,     3.27971E-18,    0.,             0.,             0.,             0.],
                                       [3.78692E-9,     -1.66075E-7,    -3.29971E-8,    0.,             0.,             0.],
                                       [9.39581E-12,    -1.85616E-12,   9.90046E-12,    -4.16255E-13,   0.,             0.],
                                       [-7.05186E-16,   3.50820E-16,    -8.74004E-16,   3.62134E-16,    -1.68652E-16,   0.],
                                       [1.40799E-19,    1.11211E-20,    2.88368E-19,    1.01389E-20,    1.44172E-19,    4.33434E-22]],
                            'DIST_B': [[0.,             0.,             0.,             0.,             0.,             0.],
                                       [-2.73E-5,       3.08613E-2,     0.,             0.,             0.,             0.],
                                       [8.38352E-8,     3.85602E-8,     -8.46336E-8,    0.,             0.,             0.],
                                       [-7.68727E-13,   8.76109E-12,    -2.17489E-12,   9.65818E-12,    0.,             0.],
                                       [2.02700E-16,    -5.78E-16,      7.68233E-16,    -5.22E-16,      5.32993E-16,    0.],
                                       [3.44618E-21,    1.56499E-19,    1.54886E-20,    2.93225E-19,    1.10094E-20,    1.34691E-19]]}#,
                 }
    BACKGROUND = {
                    'none': {   'F070W': 0., 'F090W': 0., 'F115W': 0., 'F140M': 0., 'F150W': 0., 'F162M': 0., 'F164N': 0., 
                                'F182M': 0., 'F187N': 0., 'F200W': 0., 'F210M': 0., 'F212N': 0., 'F250M': 0., 'F277W': 0., 
                                'F300M': 0., 'F322W': 0., 'F323N': 0., 'F335M': 0., 'F356W': 0., 'F360M': 0., 'F405N': 0., 
                                'F410M': 0., 'F430M': 0., 'F444W': 0., 'F460M': 0., 'F466N': 0., 'F470N': 0., 'F480M': 0.},
                    'avg':  {   'F070W': 9.416E-02, 'F090W': 1.523E-01, 'F115W': 1.924E-01, 'F140M': 7.563E-02, 'F150W': 2.325E-01, 
                                'F162M': 7.287E-02, 'F164N': 6.918E-03, 'F182M': 7.709E-02, 'F187N': 5.899E-03, 'F200W': 1.881E-01, 
                                'F210M': 6.148E-02, 'F212N': 5.379E-03, 'F250M': 1.224E-01, 'F277W': 4.446E-01, 'F300M': 1.180E-01, 
                                'F322W': 8.741E-01, 'F323N': 8.456E-03, 'F335M': 1.143E-01, 'F356W': 3.959E-01, 'F360M': 1.418E-01, 
                                'F405N': 1.909E-02, 'F410M': 2.692E-01, 'F430M': 1.627E-01, 'F444W': 1.309E+00, 'F460M': 4.544E-01, 
                                'F466N': 3.750E-02, 'F470N': 3.934E-02, 'F480M': 4.955E-01}
                 }
    BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
    BACKGROUNDS = ['None', 'Average zodiacal background', 'Median zodiacal background', 'Maximum zodiacal background', 'Minimum zodiacal background']
    BGTEXT = {'none': 'None', 'avg': 'Average zodiacal background', 'med': 'Median zodiacal background', 'max': 'Maximum zodiacal background', 'min': 'Minimum zodiacal background'}
    #PHOTFNU has units of Jy
    PHOTFNU = { 'F070W':5.085E-08, 'F090W':3.722E-08, 'F115W':3.171E-08, 'F140M':8.313E-08, 
                'F150W':2.678E-08, 'F162M':8.396E-08, 'F164N':8.724E-07, 'F182M':7.073E-08, 
                'F187N':8.645E-07, 'F200W':2.635E-08, 'F210M':7.511E-08, 'F212N':8.460E-07,
                'F250M':1.067E-07, 'F277W':2.252E-08, 'F300M':7.077E-08, 'F322W':1.166E-08, 
                'F323N':9.419E-07, 'F335M':7.286E-08, 'F356W':2.575E-08, 'F360M':7.258E-08, 
                'F405N':9.920E-07, 'F410M':7.551E-08, 'F430M':1.648E-07, 'F444W':2.535E-08,
                'F460M':8.679E-08, 'F466N':1.143E-06, 'F470N':1.160E-06, 'F480M':1.044E-07
              }
    #PHOTPLAM has units of um
    PHOTPLAM = {    'F070W':0.6967, 'F090W':0.8993, 'F115W':1.1421, 'F140M':1.3999, 'F150W':1.4962, 
                    'F162M':1.6196, 'F164N':1.6440, 'F182M':1.8196, 'F187N':1.8753, 'F200W':1.9836, 
                    'F210M':2.0974, 'F212N':2.1218, 'F250M':2.5009, 'F277W':2.7710, 'F300M':2.9966, 
                    'F322W':3.1723, 'F323N':3.2350, 'F335M':3.3402, 'F356W':3.5290, 'F360M':3.5979, 
                    'F405N':4.0524, 'F410M':4.0949, 'F430M':4.2979, 'F444W':4.4088, 'F460M':4.5847, 
                    'F466N':4.6560, 'F470N':4.7049, 'F480M':4.7900
               }
    DITHERS = ("SUBPIXEL ONLY","FULL","INTRAMODULE","INTRASCA")
    DITHER_POINTS = {
                        "SUBPIXEL ONLY": ["0"],
                        "FULL": ["3TIGHT","3","6","9","15","21","27","36","45"],
                        "INTRAMODULE": [str(x) for x in range(1,17)],
                        "INTRASCA": [str(x) for x in range(1,26)]
                     }
    DITHER_SIZE = {
                    "SUBPIXEL ONLY": ["STANDARD"],
                    "FULL": ["STANDARD"],
                    "INTRAMODULE": ["STANDARD"],
                    "INTRASCA": ["SMALL","MEDIUM","LARGE"]
                  }
    DITHER_SUBPIXEL = {
                        "SUBPIXEL ONLY": [str(x) for x in range(2,65)],
                        "FULL": [str(x) for x in range(1,65)],
                        "INTRAMODULE": [str(x) for x in range(1,65)],
                        "INTRASCA": [str(x) for x in range(1,65)]
                      }
    DITHER_OFFSETS = {
                        "FULL": {
                                    "3TIGHT":   [(-58.,-7.5),(0.,0.),(58.,7.5)],
                                    "3":    [(-58.,-23.5),(0.,0.),(58.,23.5)],
                                    "6":    [(-72.,-30.),(-43.,-18.),(-14.,-6.),
                                             (15.,6.),(44.,18.),(73.,30.)],
                                    "9":    [(-78.,-33.),(-58.,-24.),(-38.,-15.), 
                                             (-20.,-8.),(0.,0.),(20.,8.), 
                                             (38.,15.),(58.,24.),(78.,33.)],
                                    "15":   [(-77.,-32.),(-72.,-33.),(-57.,-23.),
                                             (-43.,-20.),(-37.,-14.),(-19.,-7.),
                                             (-14.,-9.),(1.,0.),(15.,3.),
                                             (21.,9.),(39.,16.),(44.,15.),
                                             (59.,25.),(73.,27.),(79.,32.)],
                                    "21":   [(-77.,-36.),(-75.,-30.),(-70.,-21.),
                                             (-55.,-23.),(-48.,-18.),(-41.,-12.),
                                             (-35.,-12.),(-19.,-5.),(-17.,-7.),
                                             (-12.,2.),(3.,0.),(10.,5.),
                                             (17.,11.),(23.,12.),(39.,18.),
                                             (41.,17.),(46.,27.),(61.,27.),
                                             (68.,34.),(75.,29.),(81.,34.)],
                                    "27":   [(-84.,-41.),(-77.,-37.),(-71.,-20.),
                                             (-64.,-32.),(-57.,-28.),(-51.,-11.),
                                             (-44.,-23.),(-37.,-19.),(-31.,-2.),
                                             (-26.,-16.),(-19.,-12.),(-13.,5.),
                                             (-6.,-9.),(1.,-5.),(7.,12.),
                                             (14.,0.),(21.,4.),(27.,21.),
                                             (32.,7.),(39.,11.),(45.,28.),
                                             (52.,16.),(59.,20.),(65.,37.),
                                             (72.,23.),(79.,27.),(85.,44.)],
                                    "36":   [(-83.,-54.),(-78.,-27.),(-77.,-37.),
                                             (-72.,-10.),(-63.,-45.),(-59.,-18.),
                                             (-57.,-28.),(-52.,-1.),(-43.,-36.),
                                             (-38.,-9.),(-38.,-19.),(-32.,8.),
                                             (-25.,-29.),(-20.,-2.),(-19.,-12.),
                                             (-14.,15.),(-5.,-22.),(0.,5.),
                                             (1.,-5.),(6.,22.),(15.,-13.),
                                             (20.,14.),(21.,4.),(26.,31.),
                                             (33.,-6),(38.,21.),(39.,11.),
                                             (44.,38.),(53.,3.),(58.,30.),
                                             (59.,20.),(64.,47.),(73.,10.),
                                             (78.,37.),(79.,27.),(84.,54.)],
                                    "45":   [(-84.,-61.),(-81.,-47.),(-76.,-20.),
                                             (-75.,-30.),(-70.,-3.),(-64.,-52.),
                                             (-61.,-38.),(-56.,-11.),(-55.,-21.),
                                             (-50.,6.),(-44.,-43.),(-41.,-29.),
                                             (-36.,-2.),(-35.,-12.),(-30.,15.),
                                             (-26.,-36.),(-23.,-22.),(-18.,5.),
                                             (-17.,-5.),(-12.,22.),(-6.,-29.),
                                             (-3.,-15.),(2.,12.),(3.,2.),
                                             (8.,20.),(14.,-20.),(17.,-6.),
                                             (22.,21.),(23.,11.),(28.,38.),
                                             (32.,-13.),(35.,1.),(40.,28.),
                                             (41.,18.),(46.,45.),(52.,-4.),
                                             (55.,10.),(60.,37.),(61.,27.),
                                             (66.,54.),(72.,3.),(75.,17.),
                                             (80.,44.),(81.,34.),(86.,61.)]
                                },
                        "INTRAMODULE":  [(-3.75,3.75),(3.75,3.75),(11.25,-11.25),
                                         (-11.25,11.25),(11.25,11.25),(-11.25,-11.25),
                                         (-3.75,-3.75),(3.75,3.75),(-11.25,3.75),
                                         (11.25,-3.75),(3.75,-11.25),(-3.75,11.25),
                                         (-11.25,-3.75),(11.25,3.75),(-3.75,-11.25),
                                         (3.75,11.25)],
                        "INTRASCA": {
                                        "SMALL":    [(0.00,0.00),(8.19,8.19),(8.19,-8.19),
                                                     (-8.19,-8.19),(-8.19,8.19),(0.00,8.19),
                                                     (8.19,0.00),(0.00,-8.19),(-8.19,0.00),
                                                     (-4.09,4.09),(4.09,4.09),(4.09,-4.09),
                                                     (-4.09,-4.09),(-4.09,0.00),(0.00,4.09),
                                                     (4.09,0.00),(0.00,-4.09),(-4.09,-8.19),
                                                     (-8.19,-4.09),(-8.19,4.09),(-4.09,8.19),
                                                     (4.09,8.19),(8.19,4.09),(8.19,-4.09),
                                                     (4.09,-8.19)],
                                        "MEDIUM":   [(0.00,0.00),(16.38,16.38),(16.38,-16.38),
                                                     (-16.38,-16.38),(-16.38,16.38),(0.00,16.38),
                                                     (16.38,0.00),(0.00,-16.38),(-16.38,0.00),
                                                     (-8.19,8.19),(8.19,8.19),(8.19,-8.19),
                                                     (-8.19,-8.19),(-8.19,0.00),(0.00,8.19),
                                                     (8.19,0.00),(0.00,-8.19),(-8.19,-16.38),
                                                     (-16.38,-8.19),(-16.38,8.19),(-8.19,16.38),
                                                     (8.19,16.38),(16.38,8.19),(16.38,-8.19),(8.19,-16.38)],
                                        "LARGE":    [(0.00,0.00),(24.56,24.56),(24.56,-24.56),(-24.56,-24.56),
                                                     (-24.56,24.56),(0.00,24.56),(24.56,0.00),(0.00,-24.56),(-24.56,0.00),
                                                     (-12.28,12.28),(12.28,12.28),(12.28,-12.28),(-12.28,-12.28),
                                                     (-12.28,0.00),(0.00,12.28),(12.28,0.00),(0.00,-12.28),(-12.28,-24.56),
                                                     (-24.56,-12.28),(-24.56,12.28),(-12.28,24.56),(12.28,24.56),
                                                     (24.56,12.28),(24.56,-12.28),(12.28,-24.56)]
                                    },
                        "SUBPIXEL": [   
                                        [(0.,0.)],
                                        [(0.,0.),(0.24,0.176)],
                                        [(0.,0.),(0.2347,0.0853),(0.1493,0.2347)],
                                        [(0.,0.),(0.16,0.144),(0.08,0.224),(0.24,0.048)],
                                        [(0.,0.),(0.1024,0.2048),(0.2048,0.4096),(0.3072,0.1024),(0.4096,0.3072)],
                                        [(0.0,0.0),(0.0747,0.1653),(0.1493,0.2667),(0.1600,0.0480),(0.2347,0.2133),
                                         (0.3093,0.1227)],
                                        [(0.0,0.0),(0.0366,0.1692),(0.1052,0.2743),(0.1417,0.0913),(0.2103,0.2606),
                                         (0.2468,0.0457),(0.3154,0.1508)],
                                        [(0.0,0.0),(0.0320,0.1440),(0.0800,0.3200),(0.1440,0.0800),(0.1680,0.2320),
                                         (0.2960,0.0560),(0.2800,0.1680),(0.3120,0.3440)],
                                        [(0.0,0.0),(0.0352,0.1387),(0.0640,0.3093),(0.1067,0.0640),(0.1707,0.1707),
                                         (0.2027,0.3413),(0.2777,0.0320),(0.3093,0.1387),(0.3413,0.3093),(0.544,0.728),
                                         (0.720,0.552),(0.496,0.504),(0.296,0.776),(0.072,0.728),(0.248,0.552),
                                         (0.024,0.504),(-0.188,0.772),(-0.412,0.724),(-0.236,0.548),(-0.460,0.500),
                                         (-0.180,0.292),(-0.404,0.244),(-0.228,0.068),(-0.452,0.020),(-0.476,-0.468),
                                         (-0.252,-0.388),(-0.396,-0.244),(-0.172,-0.164),(0.012,-0.468),(0.236,-0.388),
                                         (0.092,-0.244),(0.316,-0.164),(0.484,-0.480),(0.708,-0.400),(0.564,-0.256),
                                         (0.788,-0.176),(0.972,-0.480),(1.196,-0.400),(1.052,-0.256),(1.276,-0.176),
                                         (0.964,0.008),(1.188,0.088),(1.044,0.232),(1.268,0.312),(0.972,0.488),
                                         (1.196,0.568),(1.052,0.712),(1.276,0.792),(1.248,1.252),(1.024,1.204),
                                         (1.200,1.028),(0.976,0.980),(0.776,1.252),(0.552,1.204),(0.728,1.028),
                                         (0.504,0.980),(0.288,1.260),(0.064,1.212),(0.240,1.036),(0.016,0.988),
                                         (-0.184,1.260),(-0.408,1.212),(-0.232,1.036),(-0.456,0.988)]
                                    ]
                     }

class NIRCamShort(NIRCamBase):
    """
    The NircamShort class contains the necessary constants and modifications to run NircamShort
    observations.
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        instrument: string, which instrument this is
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Init is super init.
        """
        self.classname = self.__class__.__name__
        #Initialize superclass
        super(NIRCamShort,self).__init__(**kwargs)
                
    DETECTOR = "NIRCamShort"
    MODE = "sw_imaging"
    APERTURE = "sw"

    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    # Instrument Offset is from JWST V2,V3 centre
    INSTRUMENT_OFFSET = (-0.57,-488.61,0.)
    # Detector Offsets are RELATIVE TO INSTRUMENT OFFSET
    DETECTOR_OFFSETS = ((123.24,-33.69,1.21),(123.24,35.37,1.04),(53.67,-33.69,0.62),(53.67,35.37,0.03),
                        (-123.27,33.51,-1.04),(-123.27,-36.33,-1.21),(-53.49,33.51,-0.03),(-53.49,-36.33,-0.62))
    OFFSET_NAMES = ("A1","A2","A3","A4","B1","B2","B3","B4")
    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
    N_DETECTORS = [1, 4, 8]
    N_OFFSET = {1: (123.24, -33.69, 1.21), 4: (88.17, -1.35, 0.775), 8: (0., 0., 0.)}
    INSTRUMENT_OFFSET = (-0.57,-488.61,0.) #RA,DEC,PA
    SCALE = [0.0311, 0.0311]
    FILTERS = ('F070W','F090W','F115W','F140M','F150W','F162M','F164N','F182M','F187N','F200W','F210M','F212N')
    DEFAULT_FILTER = 'F115W'
    FLATFILE = 'err_flat_nircam.fits'
    DARKFILE = 'err_rdrk_nircam_short.fits'  # ETC short

class NIRCamLong(NIRCamBase):
    __classtype__ = "detector"
    """
    The NIRCamLong class contains the necessary constants and modifications to run NIRCamLong
    observations.
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        instrument: string, which instrument this is
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Init is super init.
        """
        self.classname = self.__class__.__name__
        #Initialize superclass
        super(NIRCamLong,self).__init__(**kwargs)
            
    DETECTOR = "NIRCamLong"
    MODE = "lw_imaging"
    APERTURE = "lw"
    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    # Instrument Offset is from JWST V2,V3 centre
    INSTRUMENT_OFFSET = (-0.57,-488.61,0.)
    # Detector Offsets are RELATIVE TO INSTRUMENT OFFSET
    DETECTOR_OFFSETS = ((88.17,1.35,1.21),(-88.17,-1.35,-1.04))
    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
    N_DETECTORS = [1, 2]
    N_OFFSET = {1: (88.17, 1.35, 1.21), 2: (0., 0., 0.)}
    OFFSET_NAMES = ("A5","B5")
    SCALE = [0.063, 0.063] #arcsec/pixel
    FILTERS = ('F250M','F277W','F300M','F323N','F335M','F356W','F360M','F405N','F410M','F430M',
               'F444W','F460M','F466N','F470N','F480M')
    DEFAULT_FILTER = 'F277W'
    FLATFILE = 'err_flat_nircam.fits'
    DARKFILE = 'err_rdrk_nircam_long.fits'  # ETC short
    DITHERS = ("SUBPIXEL ONLY","FULL","INTRAMODULE","INTRASCA")
    DITHER_OFFSETS = {
                        "FULL": {
                                    "3TIGHT":   [(-58.,-7.5),(0.,0.),(58.,7.5)],
                                    "3":    [(-58.,-23.5),(0.,0.),(58.,23.5)],
                                    "6":    [(-72.,-30.),(-43.,-18.),(-14.,-6.),
                                             (15.,6.),(44.,18.),(73.,30.)],
                                    "9":    [(-78.,-33.),(-58.,-24.),(-38.,-15.), 
                                             (-20.,-8.),(0.,0.),(20.,8.), 
                                             (38.,15.),(58.,24.),(78.,33.)],
                                    "15":   [(-77.,-32.),(-72.,-33.),(-57.,-23.),
                                             (-43.,-20.),(-37.,-14.),(-19.,-7.),
                                             (-14.,-9.),(1.,0.),(15.,3.),
                                             (21.,9.),(39.,16.),(44.,15.),
                                             (59.,25.),(73.,27.),(79.,32.)],
                                    "21":   [(-77.,-36.),(-75.,-30.),(-70.,-21.),
                                             (-55.,-23.),(-48.,-18.),(-41.,-12.),
                                             (-35.,-12.),(-19.,-5.),(-17.,-7.),
                                             (-12.,2.),(3.,0.),(10.,5.),
                                             (17.,11.),(23.,12.),(39.,18.),
                                             (41.,17.),(46.,27.),(61.,27.),
                                             (68.,34.),(75.,29.),(81.,34.)],
                                    "27":   [(-84.,-41.),(-77.,-37.),(-71.,-20.),
                                             (-64.,-32.),(-57.,-28.),(-51.,-11.),
                                             (-44.,-23.),(-37.,-19.),(-31.,-2.),
                                             (-26.,-16.),(-19.,-12.),(-13.,5.),
                                             (-6.,-9.),(1.,-5.),(7.,12.),
                                             (14.,0.),(21.,4.),(27.,21.),
                                             (32.,7.),(39.,11.),(45.,28.),
                                             (52.,16.),(59.,20.),(65.,37.),
                                             (72.,23.),(79.,27.),(85.,44.)],
                                    "36":   [(-83.,-54.),(-78.,-27.),(-77.,-37.),
                                             (-72.,-10.),(-63.,-45.),(-59.,-18.),
                                             (-57.,-28.),(-52.,-1.),(-43.,-36.),
                                             (-38.,-9.),(-38.,-19.),(-32.,8.),
                                             (-25.,-29.),(-20.,-2.),(-19.,-12.),
                                             (-14.,15.),(-5.,-22.),(0.,5.),
                                             (1.,-5.),(6.,22.),(15.,-13.),
                                             (20.,14.),(21.,4.),(26.,31.),
                                             (33.,-6),(38.,21.),(39.,11.),
                                             (44.,38.),(53.,3.),(58.,30.),
                                             (59.,20.),(64.,47.),(73.,10.),
                                             (78.,37.),(79.,27.),(84.,54.)],
                                    "45":   [(-84.,-61.),(-81.,-47.),(-76.,-20.),
                                             (-75.,-30.),(-70.,-3.),(-64.,-52.),
                                             (-61.,-38.),(-56.,-11.),(-55.,-21.),
                                             (-50.,6.),(-44.,-43.),(-41.,-29.),
                                             (-36.,-2.),(-35.,-12.),(-30.,15.),
                                             (-26.,-36.),(-23.,-22.),(-18.,5.),
                                             (-17.,-5.),(-12.,22.),(-6.,-29.),
                                             (-3.,-15.),(2.,12.),(3.,2.),
                                             (8.,20.),(14.,-20.),(17.,-6.),
                                             (22.,21.),(23.,11.),(28.,38.),
                                             (32.,-13.),(35.,1.),(40.,28.),
                                             (41.,18.),(46.,45.),(52.,-4.),
                                             (55.,10.),(60.,37.),(61.,27.),
                                             (66.,54.),(72.,3.),(75.,17.),
                                             (80.,44.),(81.,34.),(86.,61.)]
                                },
                        "INTRAMODULE":  [(-3.75,3.75),(3.75,3.75),(11.25,-11.25),
                                         (-11.25,11.25),(11.25,11.25),(-11.25,-11.25),
                                         (-3.75,-3.75),(3.75,3.75),(-11.25,3.75),
                                         (11.25,-3.75),(3.75,-11.25),(-3.75,11.25),
                                         (-11.25,-3.75),(11.25,3.75),(-3.75,-11.25),
                                         (3.75,11.25)],
                        "INTRASCA": {
                                        "SMALL":    [(0.00,0.00),(16.38,16.38),(16.38,-16.38),
                                                     (-16.38,-16.38),(-16.38,16.38),(0.00,16.38),
                                                     (16.38,0.00),(0.00,-16.38),(-16.38,0.00),
                                                     (-8.19,8.19),(8.19,8.19),(8.19,-8.19),
                                                     (-8.19,-8.19),(-8.19,0.00),(0.00,8.19),
                                                     (8.19,0.00),(0.00,-8.19),(-8.19,-16.38),
                                                     (-16.38,-8.19),(-16.38,8.19),(-8.19,16.38),
                                                     (8.19,16.38),(16.38,8.19),(16.38,-8.19),
                                                     (8.19,-16.38)],
                                        "MEDIUM":   [(0.,0.),(32.76,32.76),(32.76,-32.76),(-32.76,-32.76),
                                                     (-32.76,32.76),(0.,32.76),(32.76,0.),(0.,-32.76),
                                                     (-32.76,0.),(-16.38,16.38),(16.38,16.38),(16.38,-16.38),
                                                     (-16.38,-16.38),(-16.38,0.),(0.,16.38),(16.38,0.),
                                                     (0.,-16.38),(-16.38,-32.76),(-32.76,-16.38),(-32.76,16.38),
                                                     (-16.38,32.76),(16.38,32.76),(32.76,16.38),(32.76,-16.38),
                                                     (16.38,-32.76)],
                                        "LARGE":    [(0.,0.),(49.12,49.12),(49.12,-49.12),(-49.12,-49.12),
                                                     (-49.12,49.12),(0.,49.12),(49.12,0.),(0.,-49.12),
                                                     (-49.12,0.),(-24.56,24.56),(24.56,24.56),(24.56,-24.56),
                                                     (-24.56,-24.56),(-24.56,0.),(0.,24.56),(24.56,0.),
                                                     (0.,-24.56),(-24.56,-49.12),(-49.12,-24.56),(-49.12,24.56),
                                                     (-24.56,49.12),(24.56,49.12),(49.12,24.56),(49.12,-24.56),
                                                     (24.56,-49.12)]
                                    },
                        "SUBPIXEL": [   
                                        [(0.,0.)],
                                        [(0.,0.),(0.24,0.176)],
                                        [(0.,0.),(0.2347,0.0853),(0.1493,0.2347)],
                                        [(0.,0.),(0.16,0.144),(0.08,0.224),(0.24,0.048)],
                                        [(0.,0.),(0.1024,0.2048),(0.2048,0.4096),(0.3072,0.1024),(0.4096,0.3072)],
                                        [(0.0,0.0),(0.0747,0.1653),(0.1493,0.2667),(0.1600,0.0480),(0.2347,0.2133),
                                         (0.3093,0.1227)],
                                        [(0.0,0.0),(0.0366,0.1692),(0.1052,0.2743),(0.1417,0.0913),(0.2103,0.2606),
                                         (0.2468,0.0457),(0.3154,0.1508)],
                                        [(0.0,0.0),(0.0320,0.1440),(0.0800,0.3200),(0.1440,0.0800),(0.1680,0.2320),
                                         (0.2960,0.0560),(0.2800,0.1680),(0.3120,0.3440)],
                                        [(0.0,0.0),(0.0352,0.1387),(0.0640,0.3093),(0.1067,0.0640),(0.1707,0.1707),
                                         (0.2027,0.3413),(0.2777,0.0320),(0.3093,0.1387),(0.3413,0.3093),(0.544,0.728),
                                         (0.720,0.552),(0.496,0.504),(0.296,0.776),(0.072,0.728),(0.248,0.552),
                                         (0.024,0.504),(-0.188,0.772),(-0.412,0.724),(-0.236,0.548),(-0.460,0.500),
                                         (-0.180,0.292),(-0.404,0.244),(-0.228,0.068),(-0.452,0.020),(-0.476,-0.468),
                                         (-0.252,-0.388),(-0.396,-0.244),(-0.172,-0.164),(0.012,-0.468),(0.236,-0.388),
                                         (0.092,-0.244),(0.316,-0.164),(0.484,-0.480),(0.708,-0.400),(0.564,-0.256),
                                         (0.788,-0.176),(0.972,-0.480),(1.196,-0.400),(1.052,-0.256),(1.276,-0.176),
                                         (0.964,0.008),(1.188,0.088),(1.044,0.232),(1.268,0.312),(0.972,0.488),
                                         (1.196,0.568),(1.052,0.712),(1.276,0.792),(1.248,1.252),(1.024,1.204),
                                         (1.200,1.028),(0.976,0.980),(0.776,1.252),(0.552,1.204),(0.728,1.028),
                                         (0.504,0.980),(0.288,1.260),(0.064,1.212),(0.240,1.036),(0.016,0.988),
                                         (-0.184,1.260),(-0.408,1.212),(-0.232,1.036),(-0.456,0.988)]
                                    ]
                     }
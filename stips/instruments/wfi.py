__filetype__ = "detector"

# Local Modules
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS
from .roman_instrument import RomanInstrument
from soc_roman_tools.siaf import siaf
import numpy as np


class WFI(RomanInstrument):
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

        # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
        self.N_DETECTORS = [1]
        self.INSTRUMENT_OFFSET = (0., 0., 0.)  # Presumably there is one, but not determined
        self.DETECTOR_SIZE = (4088, 4088)  # pixels
        self.PIXEL_SIZE = 10.0  # um (Assume for now)
        self.SCALE = [0.11, 0.11]  # Assume for now
        self.FILTERS = ('F062', 'F087', 'F106', 'F129', 'F158', 'F184', 'F146', 'F213')
        self.DEFAULT_FILTER = 'F184'  # Assume for now
        self.SCA_ROTATION = -60  # Rotation of SCA with respect to SIAF

        # Get RA and DEC that will be used for detector offset calculations
        # RA = kwargs.get('ra', 0)
        DEC = kwargs.get('dec', 0)

        # Calculate the detector offsets based on the Roman SIAF file from soc_roman_tools
        # The telescope (tel) frame is defined such that (0, 0) is the center of the boresight
        # path of the WFI in the observatory coordinate system.

        # Get Roman SIAF
        rsiaf = siaf.RomanSiaf()

        # Make array of detectors
        # self.SCA_NAMES = ["WFI01_FULL", "WFI02_FULL", "WFI03_FULL", "WFI04_FULL", "WFI05_FULL",
        #                   "WFI06_FULL", "WFI07_FULL", "WFI08_FULL", "WFI09_FULL", "WFI10_FULL",
        #                   "WFI11_FULL", "WFI12_FULL", "WFI13_FULL", "WFI14_FULL", "WFI15_FULL",
        #                   "WFI16_FULL", "WFI17_FULL", "WFI18_FULL"]
        self.DETECTOR_NAMES = ["WFI01_FULL", "WFI02_FULL", "WFI03_FULL",
                               "WFI04_FULL", "WFI05_FULL", "WFI06_FULL",
                               "WFI07_FULL", "WFI08_FULL", "WFI09_FULL",
                               "WFI10_FULL", "WFI11_FULL", "WFI12_FULL",
                               "WFI13_FULL", "WFI14_FULL", "WFI15_FULL",
                               "WFI16_FULL", "WFI17_FULL", "WFI18_FULL"]
        self.SCA_NAMES = self.DETECTOR_NAMES.copy()  # for legacy compatibility?

        # Default detector offsets in (arcseconds_ra,arcseconds_dec,degrees_angle)
        self.DETECTOR_OFFSETS = [# WFI01           WFI02             WFI03
                                 (0.0, 0.0, 0.0),  (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                 # WFI04           WFI05             WFI06
                                 (0.0, 0.0, 0.0),  (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                 # WFI07           WFI08             WFI09
                                 (0.0, 0.0, 0.0),  (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                 # WFI10           WFI11             WFI12
                                 (0.0, 0.0, 0.0),  (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                 # WFI13           WFI14             WFI15
                                 (0.0, 0.0, 0.0),  (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                 # WFI16           WFI17             WFI18
                                 (0.0, 0.0, 0.0),  (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]

        # Reference point
        s01 = rsiaf['WFI01_FULL'].sci_to_tel(self.DETECTOR_SIZE[0] / 2, self.DETECTOR_SIZE[1] / 2)
        center = SkyCoord(s01[0] * u.arcsec,
                          s01[1] * u.arcsec).skyoffset_frame(self.SCA_ROTATION * u.deg)

        for i, det in enumerate(self.DETECTOR_NAMES):
            # Get V2, V3 coordinate pair in the telescope frame at the center of the SCA
            V2, V3 = rsiaf[det].sci_to_tel(self.DETECTOR_SIZE[0] / 2, self.DETECTOR_SIZE[1] / 2)
            # Transform these to the center reference point
            out = ICRS(V2 * u.arcsec, V3 * u.arcsec).transform_to(center)
            # Get out offsets in X and Y direction of SCA array, convert to arcsec
            X, Y = out.lon.value * 3600, out.lat.value * 3600
            # Offset in arcsec, arcsec, and degrees angle, converting from deg to radians
            self.DETECTOR_OFFSETS[i] = (X / np.cos(DEC / 180 * np.pi), Y, 0)

        # Copy the results into these variables needed by some functions in STIPS
        self.OFFSET_NAMES = ("WFI01", "WFI02", "WFI03", "WFI04", "WFI05", "WFI06",
                             "WFI07", "WFI08", "WFI09", "WFI10", "WFI11", "WFI12",
                             "WFI13", "WFI14", "WFI15", "WFI16", "WFI17", "WFI18")
        self.N_OFFSET = {1: self.DETECTOR_OFFSETS[0], 2: self.DETECTOR_OFFSETS[1],
                         3: self.DETECTOR_OFFSETS[2], 4: self.DETECTOR_OFFSETS[3],
                         5: self.DETECTOR_OFFSETS[4], 6: self.DETECTOR_OFFSETS[5],
                         7: self.DETECTOR_OFFSETS[6], 8: self.DETECTOR_OFFSETS[7],
                         9: self.DETECTOR_OFFSETS[8], 10: self.DETECTOR_OFFSETS[9],
                         11: self.DETECTOR_OFFSETS[10], 12: self.DETECTOR_OFFSETS[11],
                         13: self.DETECTOR_OFFSETS[12], 14: self.DETECTOR_OFFSETS[13],
                         15: self.DETECTOR_OFFSETS[14], 16: self.DETECTOR_OFFSETS[15],
                         17: self.DETECTOR_OFFSETS[16], 18: self.DETECTOR_OFFSETS[17]}

        # This is a set of offsets derived from "WFIRST-STSCI-TR1506A"
        #   Further, it assumes no rotation or imperfection (see ASCII diagrams).
        #       07             16
        #       08 04       13 17
        #       09    01 10    18
        #          05       14
        #          06 02 11 15
        #             03 12
        #
        # The following positions come from the spreadsheet
        # "GRISM off orders position_Zernikes_efficiency data_v4_20200424.xlsx" on
        # the "WSM-GRISM" tab, with values in mm. Note that for all detectors values
        # are taken from Field Position 1 (centre) for a wavelength of 1 micron.
        #
        # These are referred to as FPA positions because they are the detector
        # positions on the focal plane array, albeit with the Y axis inverted to
        # yield WFI-local co-ordinates rather than FPA-local co-ordinates.

        self.FPA_GLOBAL_ROTATION = -120.  # degrees

        # From WFIRST-STScI-TR1506A(3).pdf, 2.5mm=27.5", and 8.564mm=94.2"
        self.DETECTOR_MM_ARCSEC = 11.

        # I can't currently figure out rotation angles (or rather figure out which
        # rotation angle to use), so I'm setting everything to just zero for now.
        self.DETECTOR_ROTATION = (0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0.)

        self.PSF_INSTRUMENT = "WFI"

        # Reference Files
        self.FLATFILE = 'err_flat_wfi.fits'  # Use for the moment
        self.DARKFILE = 'err_rdrk_wfi.fits'  # IREF, IHB (use for the moment)

        # Background Values
        self.BACKGROUND = {'none': {'F062': 0., 'F087': 0., 'F106': 0., 'F129': 0.,
                                    'F158': 0., 'F184': 0., 'F146': 0., 'F213': 0.},
                           'avg':  {'F062': 0.253, 'F087': 0.254, 'F106': 0.280, 'F129': 0.270,
                                    'F158': 0.292, 'F184': 0.296, 'F146': 4.498, 'F213': 1.811}}
        self.BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
        self.BACKGROUNDS = ['None', 'Average zodiacal background', 'Median zodiacal background',
                            'Maximum zodiacal background', 'Minimum zodiacal background']
        self.BGTEXT = {'none': 'None', 'avg': 'Average zodiacal background',
                       'med': 'Median zodiacal background', 'max': 'Maximum zodiacal background',
                       'min': 'Minimum zodiacal background', 'custom': 'Custom thermal background rate',
                       'pandeia': 'Pandeia background rate'}

        # PHOTPLAM values from RDox as of June 10, 2025:
        # https://roman-docs.stsci.edu/roman-instruments-home/wfi-imaging-mode-user-guide/wfi-design/wfi-optical-elements
        self.PHOTPLAM = {'F062': 0.6291, 'F087': 0.8696, 'F106': 1.0567, 'F129': 1.2901,
                         'F158': 1.5749, 'F184': 1.8394, 'F146': 1.4378, 'F213': 2.1230}

        self.ZEROPOINTS_AB = {'F062': 26.77, 'F087': 26.43, 'F106': 26.45, 'F129': 26.47,
                              'F158': 26.50, 'F184': 26.12, 'F146': 27.70, 'F213': 26.06}

        # PHOTFNU has units of Jy
        self.PHOTFNU = {}
        for i in self.ZEROPOINTS_AB:
            self.PHOTFNU[i] = 10 ** (0.4 * (8.9 - self.ZEROPOINTS_AB[i]))

        self.DITHERS = ("SUBPIXEL ONLY", "BOX-UVIS", "BLOB")  # Assume for now
        self.DITHER_POINTS = {
                                 "SUBPIXEL ONLY": ["0"],
                                 "BOX-UVIS": ["4"],
                                 "BLOB": ["1", "2"]
                              }
        self.DITHER_SIZE = {
                             "SUBPIXEL ONLY": ["STABDARD"],
                             "BOX-UVIS": ["Standard"],
                             "BLOB": ["Standard"]
                           }
        self.DITHER_SUBPIXEL = {
                                 "SUBPIXEL ONLY": ["LINE", "LINE-3PT", "BOX-MIN"],
                                 "BOX-UVIS": ["NONE", "LINE", "LINE-3PT", "BOX-MIN"],
                                 "BLOB": ["NONE", "LINE", "LINE-3PT", "BOX-MIN"]
                               }
        self.DITHER_OFFSETS = {
                                 "BOX-UVIS": [(-11.071, -17.744), (11.947, -17.457),
                                                (11.071, 17.744), (-11.947, 17.457)],
                                 "BLOB": [(-1.930, -1.729), (1.930, 1.729)],
                                 "SUBPIXEL": {
                                                 "NONE":     [(0.000, 0.000)],
                                                 "BOX-MIN":  [(0.000, 0.000), (0.542, 0.182), (0.339, 0.485), (-0.203, 0.303)],
                                                 "LINE":     [(0.000, 0.000), (0.474, 0.424)],
                                                 "LINE-3PT": [(0.000, 0.000), (0.451, 0.403), (0.902, 0.806)]
                                             }
                              }

        # Initialize superclass
        super().__init__(**kwargs)

    def generateReadnoise(self):
        """
        Readnoise formula that is similar to HST ETC.
        """
        return 12.  # max ramp, lowest noise

    @classmethod
    def handleDithers(cls, form):
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
        else:  # dither_pattern == "SUBPIXEL ONLY" and instrument potentially matters again
            initial_dithers = [(0., 0.)]
        dither_subpixel = form['dither_subpixel']
        subpixel_dithers = cls.DITHER_OFFSETS['SUBPIXEL'][dither_subpixel]
        return cls.doSubpixel(initial_dithers, subpixel_dithers)

    INSTRUMENT = "WFI"
    DETECTOR = "WFI"

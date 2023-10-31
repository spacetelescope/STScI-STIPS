__filetype__ = "detector"

# Local Modules
from .roman_instrument import RomanInstrument


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

    OFFSET_NAMES = ("SCA01", "SCA02", "SCA03", "SCA04", "SCA05", "SCA06",
                    "SCA07", "SCA08", "SCA09", "SCA10", "SCA11", "SCA12",
                    "SCA13", "SCA14", "SCA15", "SCA16", "SCA17", "SCA18")
    N_OFFSET = {1:      (0.0,    0.0,  0.0),  2:     (0.0,  464.31, 0.0),  3:     (0.0,   995.32, 0.0),
                4:   (-533.06, 131.74, 0.0),  5:  (-533.06, 596.05, 0.0),  6:  (-533.06, 1127.06, 0.0),
                7:  (-1066.12, 338.76, 0.0),  8: (-1066.12, 803.07, 0.0),  9: (-1066.12, 1334.08, 0.0),
                10:   (533.06,   0.0,  0.0), 11:   (533.06, 464.31, 0.0), 12:   (533.06,  995.32, 0.0),
                13:  (1066.12, 131.74, 0.0), 14:  (1066.12, 596.05, 0.0), 15:  (1066.12, 1127.06, 0.0),
                16:  (1599.18, 338.76, 0.0), 17:  (1599.18, 803.07, 0.0), 18:  (1599.18, 1334.08, 0.0)}

    # Each detector is 450.56 x 450.56", and according to the WFIRST-STSCI-TR1506A
    # Document, the offsets are (27.5”) in the long direction and (94.2" and 27.5")
    # in the short direction. Assuming 0 offset for SCA01. Then RA offsets are equal
    # to 478.06", and DEC offsets are 478.06" for the lower SCAs and 544.76 for the
    # upper ones. Assuming no rotation offset.

    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    '''
    DETECTOR_OFFSETS = (# SCA01                                        SCA02                                         SCA03                                         SCA04
                        ((2*27.5+478.06)*( 0),(478.06-27.5/2)*( 0),0.),            ((2*27.5+478.06)*( 0),(478.06-27.5/2)*( 1),0.),            ((2*27.5+478.06)*( 0),(478.06-27.5/2)*( 2)+66.7,0.),       ((2*27.5+478.06)*(-1),(478.06-27.5/2)*( 0)+188.2*0.7,0.),
                        # SCA05                                        SCA06                                         SCA07                                         SCA08
                        ((2*27.5+478.06)*(-1),(478.06-27.5/2)*( 1)+188.2*0.7,0.),      ((2*27.5+478.06)*(-1),(478.06-27.5/2)*( 2)+66.7+188.2*0.7,0.), ((2*27.5+478.06)*(-2),(478.06-27.5/2)*( 0)+376.4*0.9,0.),      ((2*27.5+478.06)*(-2),(478.06-27.5/2)*( 1)+376.4*0.9,0.),
                        # SCA09                                        SCA10                                         SCA11                                         SCA12
                        ((2*27.5+478.06)*(-2),(478.06-27.5/2)*( 2)+66.7+376.4*0.9,0.), ((2*27.5+478.06)*( 1),(478.06-27.5/2)*( 0),0.),            ((2*27.5+478.06)*( 1),(478.06-27.5/2)*( 1),0.),            ((2*27.5+478.06)*( 1),(478.06-27.5/2)*( 2)+66.7,0.),
                        # SCA13                                        SCA14                                         SCA15                                         SCA16
                        ((2*27.5+478.06)*( 2),(478.06-27.5/2)*( 0)+188.2*0.7,0.),      ((2*27.5+478.06)*( 2),(478.06-27.5/2)*( 1)+188.2*0.7,0.),      ((2*27.5+478.06)*( 2),(478.06-27.5/2)*( 2)+66.7+188.2*0.7,0.), ((2*27.5+478.06)*( 3),(478.06-27.5/2)*( 0)+376.4*0.9,0.),
                        # SCA17                                        SCA18
                        ((2*27.5+478.06)*( 3),(478.06-27.5/2)*( 1)+376.4*0.9,0.),      ((2*27.5+478.06)*( 3),(478.06-27.5/2)*( 2)+66.7+376.4*0.9,0.))
    '''

    DETECTOR_OFFSETS = (  # SCA01                  SCA02                    SCA03
                        (0.0, 0.0, 0.0),         (0.0, 464.31, 0.0),      (0.0, 995.32, 0.0),
                        # SCA04                  SCA05                    SCA06
                        (-533.06, 131.74, 0.0),  (-533.06, 596.05, 0.0),  (-533.06, 1127.06, 0.0),
                        # SCA07                  SCA08                    SCA09
                        (-1066.12, 338.76, 0.0), (-1066.12, 803.07, 0.0), (-1066.12, 1334.08, 0.0),
                        # SCA10                  SCA11                    SCA12
                        (533.06, 0.0, 0.0),      (533.06, 464.31, 0.0),   (533.06, 995.32, 0.0),
                        # SCA13                  SCA14                    SCA15
                        (1066.12, 131.74, 0.0),  (1066.12, 596.05, 0.0),  (1066.12, 1127.06, 0.0),
                        # SCA16                  SCA17                    SCA18
                        (1599.18, 338.76, 0.0),  (1599.18, 803.07, 0.0),  (1599.18, 1334.08, 0.0))

    # This is a set of offsets derived from "WFIRST-STSCI-TR1506A"
    #
    # Since that document doesn't actually cover the column offsets, they are
    #   assumed to be 188.2" e.g. between 08 and 07, and 376.4" e.g. between
    #   09 and 08 (and the same for each other column)(see ASCII diagrams).
    #   Further, it assumes no rotation or imperfection.
    #       09             18
    #
    #       08 06       15 17
    #       07    03 12    16
    #          05       14
    #          04 02 11 13
    #             01 10

    # The detector positions WRT to WFI Local FOV are as follows (with the
    # X-axis pointing right and the Y-axis pointing up):
    #
    #             01 10
    #          04 02 11 13
    #          05       14
    #       07    03 12    16
    #       08 06       15 17
    #
    #       09             18
    #
    # with the large gaps (e.g. 09-08) showing the larger offset between the
    # bottom row of detectors and the top two rows, or, in the case of the
    # offsets between columns (e.g. the Y positions of 07 vs. 04 vs. 01) the
    # staggered offsets of the columns.
    #
    # The following positions come from the spreadsheet
    # "GRISM off orders position_Zernikes_efficiency data_v4_20200424.xlsx" on
    # the "WSM-GRISM" tab, with values in mm. Note that for all detectors values
    # are taken from Field Position 1 (centre) for a wavelength of 1 micron.
    #
    # These are referred to as FPA positions because they are the detector
    # positions on the focal plane array, albeit with the Y axis inverted to
    # yield WFI-local co-ordinates rather than FPA-local co-ordinates.

    FPA_GLOBAL_ROTATION = -120.  # degrees

    # Values in mm. From "FPA Position" X and Y columns (F/G in xlsx)
    # Specifically, 1.55 um is supposed to be the "undistorted" wavelength, so
    # I'm using the detector centre at 1.55um.
    DETECTOR_FPA = (  # SCA01                SCA02               SCA03
                      (-22.029, -11.956), (-22.181, 36.421), (-22.331, 80.749),
                    # SCA04                SCA05               SCA06
                      (-65.558, -20.505), (-66.055, 27.856), (-66.543, 71.928),
                    # SCA07                SCA08               SCA09
                      (-109.05,  -41.316), (-109.83,   7.001), (-110.97,  50.358),
                    # SCA10                SCA11               SCA12
                      (21.509, -11.955), (21.65,  36.42),  (21.787, 80.747),
                    # SCA13                SCA14               SCA15
                      (65.03,  -20.503), (65.517, 27.854), (65.993, 71.923),
                    # SCA16                SCA17               SCA18
                      (108.507, -41.309), (109.282,  7.002), (110.409, 50.353),
                   )

    # From WFIRST-STScI-TR1506A(3).pdf, 2.5mm=27.5", and 8.564mm=94.2"
    DETECTOR_MM_ARCSEC = 11.

    # I can't currently figure out rotation angles (or rather figure out which
    # rotation angle to use), so I'm setting everything to just zero for now.
    DETECTOR_ROTATION = (0., 0., 0., 0., 0., 0., 0., 0., 0.,
                         0., 0., 0., 0., 0., 0., 0., 0., 0.)

    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
    N_DETECTORS = [1]
    INSTRUMENT_OFFSET = (0., 0., 0.)  # Presumably there is one, but not determined
    DETECTOR_SIZE = (4088, 4088)  # pixels
    PIXEL_SIZE = 10.0  # um (Assume for now)
    SCALE = [0.11, 0.11]  # Assume for now
    FILTERS = ('F062', 'F087', 'F106', 'F129', 'F158', 'F184', 'F146')
    DEFAULT_FILTER = 'F184'  # Assume for now

    PSF_INSTRUMENT = "WFI"

    # Reference Files
    FLATFILE = 'err_flat_wfi.fits'  # Use for the moment
    DARKFILE = 'err_rdrk_wfi.fits'  # IREF, IHB (use for the moment)

    # Background Values
    BACKGROUND = {'none': {'F062': 0., 'F087': 0., 'F106': 0., 'F129': 0., 'F158': 0., 'F184': 0., 'F146': 0.},
                  'avg':  {'F062': 1.401E+00, 'F087': 1.401E+00, 'F106': 1.401E+00, 'F129': 7.000E-01,
                           'F158': 7.521E-01, 'F184': 8.500E-01, 'F146': 7.000E-01}}
    BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
    BACKGROUNDS = ['None', 'Average zodiacal background', 'Median zodiacal background', 'Maximum zodiacal background', 'Minimum zodiacal background']
    BGTEXT = {'none': 'None', 'avg': 'Average zodiacal background',
              'med': 'Median zodiacal background', 'max': 'Maximum zodiacal background',
              'min': 'Minimum zodiacal background', 'custom': 'Custom thermal background rate'}
    # PHOTFNU has units of Jy
    # For now, just assuming similar PHOTFNU to WFC3IR.
    PHOTFNU = {'F062': 3.25e-08, 'F087': 8.87e-08, 'F106': 3.94e-08, 'F129': 3.51e-08, 'F158': 3.13e-08, 'F184': 1.18e-07, 'F146': 1.63e-08}
    # PHOTPLAM has units of um
    # For now, just put them in the middle
    PHOTPLAM = {'F062': 0.6700, 'F087': 0.8735, 'F106': 1.0595, 'F129': 1.2925, 'F158': 1.577, 'F184': 1.5815, 'F146': 1.4635}
    # For now, just put in HST-style dithers.
    DITHERS = ("SUBPIXEL ONLY", "BOX-UVIS", "BLOB")  # Assume for now
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
                        "SUBPIXEL ONLY": ["LINE", "LINE-3PT", "BOX-MIN"],
                        "BOX-UVIS": ["NONE", "LINE", "LINE-3PT", "BOX-MIN"],
                        "BLOB": ["NONE", "LINE", "LINE-3PT", "BOX-MIN"]
                      }
    DITHER_OFFSETS = {
                        "BOX-UVIS": [(-11.071, -17.744), (11.947, -17.457), (11.071, 17.744), (-11.947, 17.457)],
                        "BLOB": [(-1.930, -1.729), (1.930, 1.729)],
                        "SUBPIXEL": {
                                        "NONE":     [(0.000, 0.000)],
                                        "BOX-MIN":  [(0.000, 0.000), (0.542, 0.182), (0.339, 0.485), (-0.203, 0.303)],
                                        "LINE":     [(0.000, 0.000), (0.474, 0.424)],
                                        "LINE-3PT": [(0.000, 0.000), (0.451, 0.403), (0.902, 0.806)]
                                    }
                     }

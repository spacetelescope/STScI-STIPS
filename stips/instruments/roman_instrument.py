__filetype__ = "base"

# Local Modules
from .instrument import Instrument


class RomanInstrument(Instrument):
    """
    The RomanInstrument class contains the necessary constants and modifications specific to Roman
        but independent of any specific instrument. This class is also not intended to be
        implemented directly, but rather through its children (e.g. WFI).
        It contains the following constants and the following new variables:

        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Still only a base class. Init is super init.
        """
        super().__init__(**kwargs)
        self.REFS = {
                        'comptable': self.COMPFILES[-1],
                        'graphtable': self.GRAPHFILES[-1],
                        'thermtable': self.THERMFILES[-1],
                        'area': 45238.93416,
                        'waveset': (500, 26000, 10000., 'log')
                    }

    TELESCOPE = 'ROMAN'
    # WARNING : The definition of the area is currently under discussion, the 
    # value presented here is simply pi * r^2, where r = 0.6m, half the diameter
    # of the Roman mirror. This is subject to change.
    AREA = 45238.93416
    DBNAME = "IsochroneGrid.db"
    MODE = 'imaging'
    APERTURE = 'any'
    PANDEIA_VERSION = '0.1.dev0'

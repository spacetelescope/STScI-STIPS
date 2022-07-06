from __future__ import absolute_import, division
__filetype__ = "base"

#Local Modules
from .instrument import Instrument

class JwstInstrument(Instrument):
    """
    The JwstInstrument class contains the necessary constants and modifications specific to JWST
        but independent of any specific instrument. This class is also not intended to be 
        implemented directly, but rather through its children (e.g. NirCamShort, Miri, NirCamLong).
        It contains the following constants and the following new variables:
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Still only a base class. Init is super init.
        """
        super(JwstInstrument,self).__init__(**kwargs)
        self.REFS = {   'comptable': self.COMPFILES[-1],
                        'graphtable': self.GRAPHFILES[-1],
                        'thermtable': self.THERMFILES[-1],
                        'area': 253260.0,
                        'waveset': (5000,260000,10000.,'log')}
            
    TELESCOPE = 'JWST'
    AREA = 253260.0
    DBNAME = "IsochroneGrid.db"
    MODE = 'imaging'
    APERTURE = 'imager'
    PANDEIA_VERSION = '1.1.1'

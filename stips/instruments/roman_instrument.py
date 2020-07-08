from __future__ import absolute_import,division
__filetype__ = "base"

#Local Modules
from .instrument import Instrument

class RomanInstrument(Instrument):
    """
    The RomanInstrument class contains the necessary constants and modifications 
    specific to the Roman Space Telescope (Roman) but independent of any 
    specific instrument. This class is also not intended to be implemented 
    directly, but rather through its children (e.g. WFI).
    """

    def __init__(self, **kwargs):
        """
        Still only a base class. Init is super init.
        """
        super(WfirstInstrument,self).__init__(**kwargs)
        self.REFS = {   'comptable': self.COMPFILES[-1],
                        'graphtable': self.GRAPHFILES[-1],
                        'thermtable': self.THERMFILES[-1],
                        'area': 45238.93416,
                        'waveset': (500,26000,10000.,'log')
                    }
    
    TELESCOPE = 'WFIRST'
    AREA = 45238.93416
    DBNAME = "IsochroneGrid.db"
    MODE = 'imaging'
    APERTURE = 'any'
    PANDEIA_VERSION = '0.1.dev0'

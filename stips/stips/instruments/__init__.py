__all__ = ['hst_instrument', 'instrument', 'jwst_instrument', 'miri', 'nircam', 'wfc3ir', 'wfi', 'wfirst_instrument']
# Don't modify the line above, or this line!
import automodinit
automodinit.automodinit(__name__, __file__, globals())
del automodinit
# Local Import
from .instrument import Instrument

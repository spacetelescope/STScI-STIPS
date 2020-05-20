from __future__ import absolute_import

# External modules
import glob, logging, os, sys

# Local modules
from ..utilities import GetStipsData, InstrumentList, OffsetPosition, StipsDataTable
from ..astro_image import AstroImage

#-----------
class ObservationModule(object):

    #-----------
    def __init__(self, obs, **kwargs):
        """
        JWST instrument observation module

        :Author: Brian York

        :Organization: Space Telescope Science Institute

        :History:
            * 2010/10/19 PLL created this module as SceneModule.py
            * 2011/06/14 PLL added single star simulation.
            * 2011/06/28 PLL reorganized functions.
            * 2011/10/28 PLL added galaxies simulation.
            * 2014/03/05 BY split out the Observation components creating ObservationModule.py
            * 2014/12/03 BAQ revised the ObservationModule to use instrument classes

        ----------
        self: obj
            Class instance.
            
        kwargs: dictionary
            Contains the necessary instrument and scene parameters
        """
        # Initialize logger
        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            self.logger = logging.getLogger('__stips__')
            self.logger.setLevel(logging.INFO)
            if not len(self.logger.handlers):
                stream_handler = logging.StreamHandler(sys.stderr)
                stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))# [in %(pathname)s:%(lineno)d]'))
                self.logger.addHandler(stream_handler)
        
        #initialize parameters from the supplied observation
        self.instrument_name = obs.get('instrument', 'NIRCamShort')
        self.filters = obs.get('filters', [])
        self.offsets = obs.get('offsets', [])
        self._log('info', "Got offsets as {}".format(self.offsets))
        self.oversample = int(obs.get('oversample', 1))
        self.psf_commands = obs.get('pupil_mask', "")
        self.id = obs.get('observations_id', '0')
        self.detectors = int(obs.get('detectors', 1))
        self.excludes = obs.get('excludes', [])
        self.exptime = float(obs.get('exptime', 1.))
        self.distortion = obs.get('distortion', False)
        self.background = obs.get('background', 'none')
        self.custom_background = obs.get('custom_background', 0.)
        self.small_subarray = obs.get('small_subarray', False)
        if 'cache' in obs:
            self.use_local_cache = obs['use_local_cache']
        if len(self.filters) == 0 and 'filter' in obs:
            self.filters.append(obs['filter'])
        
        #initialize parameters from the supplied keyword arguments
        self.prefix = kwargs.get('out_prefix', 'sim')
        self.cat_path = kwargs.get('cat_path', os.getcwd())
        self.out_path = kwargs.get('out_path', os.getcwd())
        self.cat_type = kwargs.get('cat_type', 'fits')
        self.convolve_size = kwargs.get('convolve_size', 8192)
        self.parallel = kwargs.get('parallel', False)
        self.cores = kwargs.get('cores', None)
        
        if 'scene_general' in kwargs:
            self.ra = kwargs['scene_general'].get('ra', 0.0)
            self.dec = kwargs['scene_general'].get('dec', 0.0)
            self.pa = kwargs['scene_general'].get('pa', 0.0)
            self.seed = kwargs['scene_general'].get('seed', 0)
            if 'cache' in kwargs['scene_general'] and not hasattr(self, 'use_local_cache'):
                self.use_local_cache = kwargs['scene_general']['cache']
        else:
            self.ra = kwargs.get('ra', 0.0)
            self.dec = kwargs.get('dec', 0.0)
            self.pa = kwargs.get('pa', 0.0)
            self.seed = kwargs.get('seed', 0)
            if not hasattr(self, 'use_local_cache'):
                self.use_local_cache = kwargs.get('cache', False)
        if 'residual' in kwargs:
            self.flat = kwargs['residual'].get('flat', True)
            self.dark = kwargs['residual'].get('dark', True)
            self.cosmic = kwargs['residual'].get('cosmic', True)
            self.poisson = kwargs['residual'].get('poisson', True)
            self.readnoise = kwargs['residual'].get('readnoise', True)
        else:
            self.flat = kwargs.get('flatfield', True)
            self.dark = kwargs.get('dark', True)
            self.cosmic = kwargs.get('cosmic', True)
            self.poisson = kwargs.get('poisson', True)
            self.readnoise = kwargs.get('readnoise', True)
        self.version = kwargs.get('version', '0.0')
        self.memmap = kwargs.get('memmap', True)
        self.set_celery = kwargs.get('set_celery', None)
        self.get_celery = kwargs.get('get_celery', None)

        if sys.version_info[0] >= 3:
            if isinstance(self.instrument_name, bytes):
                self.instrument_name = self.instrument_name.decode('utf8')
        else:
            self.instrument_name = self.instrument_name.encode('ascii')
        if self.instrument_name == 'WFIRST':
            self.instrument_name = 'WFI'
        
        self.observations = []
        for filter in self.filters:
            for offset in self.offsets:
                self._log("info","Adding observation with filter %s and offset (%f,%f,%f)" % (filter,offset['offset_ra'],offset['offset_dec'],offset['offset_pa']))
                self.observations.append({'filter':filter,'offset':offset})
        self._log("info","Added %d observations" % (len(self.observations)))
        self.obs_count = -1 #initially. Will advance to 0 when nextObservation() is first called.
        self.images = {}
        
        self.imgbase = os.path.join(self.out_path, "{}_{}".format(self.prefix, self.id))
        
        self.instruments = InstrumentList(excludes=self.excludes)
        self.instrument = self.instruments[str(self.instrument_name)](**self.__dict__)
    
    #-----------
    def initParams(self):
        """
        Initialize the parsed parameters list
        
        Parameters
        ----------
        self: obj
            Class instance
        """
        self.params = [
            'Instrument: %s' % self.instrument_name,
            'Filters: %s' % self.instrument.filter,
            'Pixel scale: (%.3f,%.3f) arcsec/pix' % (self.instrument.SCALE[0],self.instrument.SCALE[1]),
            'Pivot wavelength: %.3f micron' % self.instrument.photplam,
            'Background Value: %s' % self.instrument.BGTEXT[self.background],
            'Exposure time: %s s' % self.exptime,
            'Input unit: counts/s',
            'Flatfield correction error applied? %s' % self.flat,
            'Dark correction error applied? %s' % self.dark,
            'Cosmic rays correction error applied? %s' % self.cosmic,
            'X size: %i pix' % self.instrument.DETECTOR_SIZE[0],
            'Y size: %i pix' % self.instrument.DETECTOR_SIZE[1],
            'Random seed: %s' % self.seed]
    
    #-----------
    def prepImage(self, img):
        """
        Initialize an image (to be added in later). Saves time to initialize each once.
        """
        imname = img['current_filename']
        imscale = float(img['scale'])
        if img['wcs']:
            self.images[imname] = AstroImage.initFromFits(imname,ext=img['ext'],logger=self.logger)
        else:
            self.images[imname] = AstroImage.initDataFromFits(imname,ext=img['ext'],ra=self.ra,dec=self.dec,pa=self.pa,scale=[imscale,imscale],logger=self.logger)
    
    #-----------
    def nextObservation(self):
        """
        Given that an observation can contain a set of offsets and filters, this moves to the next
        observation in that set, and re-initializes the local instrument. If possible, it keeps the
        same calculated PSF.
        """
        self.obs_count += 1
        self._log("info","Initializing Observation %d of %d" % (self.obs_count,len(self.observations)))
        if self.obs_count < len(self.observations):
            filter = self.observations[self.obs_count]['filter']
            self._log("info","Observation Filter is %s" % (filter))
            offset_ra = float(self.observations[self.obs_count]['offset']['offset_ra'])/3600.
            offset_dec = float(self.observations[self.obs_count]['offset']['offset_dec'])/3600.
            offset_pa = float(self.observations[self.obs_count]['offset']['offset_pa'])
            ra,dec = OffsetPosition(self.ra,self.dec,float(offset_ra),float(offset_dec))
            pa = (self.pa + float(offset_pa))%360.
            if self.observations[self.obs_count]['offset']['offset_centre']:
                offset_ra = self.instrument.INSTRUMENT_OFFSET[0]/3600.
                offset_dec = self.instrument.INSTRUMENT_OFFSET[1]/3600.
                offset_pa = self.instrument.INSTRUMENT_OFFSET[2]
                ra,dec = OffsetPosition(ra,dec,offset_ra,offset_dec)
                pa = (pa + offset_pa)%360.
            self._log("info","Observation (RA,DEC) = (%f,%f) with PA=%f" % (ra,dec,pa))
            self.instrument.reset(ra, dec, pa, filter, self.obs_count)
            self._log("info","Reset Instrument")
            self.initParams()
            return self.obs_count
        else:
            return None
    
    #-----------
    def addCatalogue(self, catalogue, *args, **kwargs):
        """
        Add a catalogue to the internal image.
        
        Parameters
        ----------
        
        self: obj
            Class instance.
        
        catalogue: string
            Name of catalogue file
        """
        self._log("info","Running catalogue %s" % (catalogue))
        if 'parallel' not in kwargs:
            kwargs['parallel'] = self.parallel
        if 'cores' not in kwargs:
            kwargs['cores'] = self.cores
        cats = self.instrument.addCatalogue(catalogue, self.id, *args, **kwargs)
        self._log("info",'Finished catalogue %s' % (catalogue))
        return cats

    #-----------
    def addTable(self, table, table_type, *args, **kwargs):
        """
        Add an in-memory data table to the internal image
        
        Parameters
        ----------
        
        self: obj
            Class instance.
        
        table: astropy.table.Table
            table of observation data
        
        table_type: string
            what type of table it is
        """
        self._log("info","Running {} table".format(table_type))
        if 'parallel' not in kwargs:
            kwargs['parallel'] = self.parallel
        if 'cores' not in kwargs:
            kwargs['cores'] = self.cores
        tables = self.instrument.addTable(table, table_type, *args, **kwargs)
        self._log("info", "Finished {} table".format(table_type))
        return tables

    #-----------
    def addImage(self, img, units, *args, **kwargs):
        """
        Create an internal image from a user-provided file
        
        Parameters
        ----------
        
        self: obj
            Class instance.
        
        file: string
            File to read
        
        units: string
            Unit type of file
        """
        self._log("info",'Adding image %s to observation' % (file))
#        if True:
        if self.images[img].header["ASTROIMAGEVALID"]:
            self.instrument.addImage(self.images[img], units, *args, **kwargs)
        self._log("info",'Image Added')

    #-----------
    def addError(self, *args, **kwargs):
        """
        Add internal sources of error to the image
        
        Parameters
        ----------
        
        self: obj
            Class instance.
        """
        psf_names = []
        psf_path = os.path.join(self.out_path, "psf_cache")
        if os.path.exists(psf_path):
            psf_names = glob.glob(os.path.join(psf_path, "*.fits"))
        self._log("info","Adding Error")
        if 'parallel' not in kwargs:
            kwargs['parallel'] = self.parallel
        if 'cores' not in kwargs:
            kwargs['cores'] = self.cores
        self.instrument.addError(poisson=self.poisson, readnoise=self.readnoise, 
                                 flat=self.flat, dark=self.dark, 
                                 cosmic=self.cosmic, *args, **kwargs)
        self._log("info","Finished Adding Error")
        return psf_names
        
    #-----------
    def finalize(self, mosaic=True, *args, **kwargs):
        """
        Finalize FITS file
        
        Parameters
        ----------
        
        self: obj
            Class instance.
        """
        self.instrument.toFits("%s_%d.fits" % (self.imgbase, self.obs_count), *args, **kwargs)
        mosaics = None
        if mosaic:
            mosaics = self.instrument.toMosaic("%s_%d_mosaic.fits"%(self.imgbase, self.obs_count), *args, **kwargs)
        return "%s_%d.fits"%(self.imgbase, self.obs_count), mosaics, self.params
    
    #-----------
    def totalObservations(self):
        """
        Return the total number of observations
        """
        return len(self.observations)

    def _log(self,mtype,message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self,'logger'):
            getattr(self.logger,mtype)(message)
        else:
            sys.stderr.write("%s: %s\n" % (mtype,message))
    
    def updateState(self, state):
        if self.set_celery is not None:
            self.set_celery(state)
    
    def getState(self):
        if self.get_celery is not None:
            return self.get_celery()
        return ""

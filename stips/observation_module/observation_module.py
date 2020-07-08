from __future__ import absolute_import

# External modules
import glob, logging, os, sys

# Local modules
from ..utilities import GetStipsData
from ..utilities import InstrumentList
from ..utilities import OffsetPosition
from ..utilities import StipsDataTable
from ..utilities import SelectParameter
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
            log_level = SetParameter('log_level', kwargs)
            self.logger.setLevel(getattr(logging, log_level))
            if not len(self.logger.handlers):
                stream_handler = logging.StreamHandler(sys.stderr)
                stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))# [in %(pathname)s:%(lineno)d]'))
                self.logger.addHandler(stream_handler)
        
        #initialize parameters from the supplied observation
        self.instrument_name = obs.get('instrument', 'NIRCamShort')
        self.filters = obs.get('filters', [])
        self.offsets = obs.get('offsets', [])
        self._log('info', "Got offsets as {}".format(self.offsets))
        self.psf_commands = obs.get('pupil_mask', "")
        self.id = obs.get('observations_id', '0')
        self.detectors = int(obs.get('detectors', 1))
        self.excludes = obs.get('excludes', [])
        self.exptime = float(obs.get('exptime', 1.))
        self.background = SelectParameter('background', obs)
        self.custom_background = obs.get('custom_background', 0.)
        self.background_location = SelectParameter('jbt_location', obs)
        self.small_subarray = obs.get('small_subarray', False)
        if len(self.filters) == 0 and 'filter' in obs:
            self.filters.append(obs['filter'])
        
        #initialize parameters from the supplied keyword arguments
        self.oversample = SelectParameter('oversample', kwargs)
        self.distortion = SelectParameter('distortion', kwargs)
        self.prefix = kwargs.get('out_prefix', 'sim')
        self.cat_path = SelectParameter('cat_path', kwargs)
        self.out_path = SelectParameter('out_path', kwargs)
        self.cat_type = SelectParameter('cat_type', kwargs)
        self.convolve_size = SelectParameter('convolve_size', kwargs)
        self.parallel = SelectParameter('parallel', kwargs)
        self.cores = SelectParameter('cores', kwargs)
        self.psf_grid_size = SelectParameter('psf_grid_size', kwargs)
        self.memmap = SelectParameter('memmap', kwargs)
        
        if 'scene_general' in kwargs:
            self.ra = kwargs['scene_general'].get('ra', 0.0)
            self.dec = kwargs['scene_general'].get('dec', 0.0)
            self.pa = kwargs['scene_general'].get('pa', 0.0)
            self.seed = SelectParameter('seed', kwargs['scene_general'])
        else:
            self.ra = kwargs.get('ra', 0.0)
            self.dec = kwargs.get('dec', 0.0)
            self.pa = kwargs.get('pa', 0.0)
            self.seed = SelectParameter('seed', kwargs)
        if 'residual' in kwargs:
            kwdict = kwargs['residual']
            self.residual_flat = SelectParameter('residual_flat', kwdict)
            self.residual_dark = SelectParameter('residual_dark', kwdict)
            self.residual_cosmic = SelectParameter('residual_cosmic', kwdict)
            self.residual_poisson = SelectParameter('residual_poisson', kwdict)
            self.residual_readnoise = SelectParameter('residual_readnoise', 
                                                      kwdict)
        else:
            self.residual_flat = SelectParameter('residual_flat', kwargs)
            self.residual_dark = SelectParameter('residual_dark', kwargs)
            self.residual_cosmic = SelectParameter('residual_cosmic', kwargs)
            self.residual_poisson = SelectParameter('residual_poisson', kwargs)
            self.residual_readnoise = SelectParameter('residual_readnoise', 
                                                      kwargs)

        # Parameters used if called from a celery queue manager
        self.set_celery = kwargs.get('set_celery', None)
        self.get_celery = kwargs.get('get_celery', None)

        if sys.version_info[0] >= 3:
            if isinstance(self.instrument_name, bytes):
                self.instrument_name = self.instrument_name.decode('utf8')
        else:
            self.instrument_name = self.instrument_name.encode('ascii')
        if self.instrument_name == 'WFIRST' or self.instrument_name == 'ROMAN':
            self.instrument_name = 'WFI'
        
        self.observations = []
        for filter in self.filters:
            for offset in self.offsets:
                msg = "Adding observation with filter {} and offset ({},{},{})"
                self._log("info",msg.format(filter, offset['offset_ra'],
                                            offset['offset_dec'], 
                                            offset['offset_pa']))
                self.observations.append({'filter':filter,'offset':offset})
        self._log("info","Added %d observations" % (len(self.observations)))

        # Initially
        self.obs_count = -1
        self.images = {}
        
        self.imgbase = os.path.join(self.out_path, "{}_{}".format(self.prefix, self.id))
        
        self.instruments = InstrumentList(excludes=self.excludes)
        instrument_class = self.instruments[str(self.instrument_name)]
        self.instrument = instrument_class(**self.__dict__)
    

    #-----------
    def initParams(self):
        """
        Initialize the parsed parameters list
        
        Parameters
        ----------
        self: obj
            Class instance
        """
        scale = (self.instrument.SCALE[0],self.instrument.SCALE[1])
        try:
            bg = self.instrument.BGTEXT[self.background]
        except Exception as e:
            self.logger.error("Error Parsing Background: {}".format(e))
            bg = 'None'
        self.params = [
            'Instrument: {}'.format(self.instrument_name),
            'Filters: {}'.format(self.instrument.filter),
            'Pixel Scale: ({:.3f},{:.3f}) arcsec/pix'.format(scale[0], 
                                                             scale[1]),
            'Pivot Wavelength: {:.3f} um'.format(self.instrument.photplam),
            'Background Type: {}'.format(bg),
            'Exposure Time: {}s'.format(self.exptime),
            'Input Unit: counts/s',
            'Added Flatfield Residual: {}'.format(self.residual_flat),
            'Added Dark Current Residual: {}'.format(self.residual_dark),
            'Added Cosmic Ray Residual: {}'.format(self.residual_cosmic),
            'Added Readnoise: {}'.format(self.residual_readnoise),
            'Added Poisson Noise: {}'.format(self.residual_poisson),
            'Detector X size: {}'.format(self.instrument.DETECTOR_SIZE[0]),
            'Detector Y size: {}'.format(self.instrument.DETECTOR_SIZE[1]),
            'Random Seed: {}'.format(self.seed)]

    
    #-----------
    def prepImage(self, img):
        """
        Initialize an image (to be added in later). Saves time to initialize each once.
        """
        imname = img['current_filename']
        imscale = float(img['scale'])
        if img['wcs']:
            img = AstroImage.initFromFits(imname, ext=img['ext'], psf=False,
                                          logger=self.logger)
        else:
            scale = [imscale, imscale]
            img = AstroImage.initDataFromFits(imname, ext=img['ext'], psf=False,
                                              ra=self.ra, dec=self.dec, 
                                              pa=self.pa, scale=scale,
                                              logger=self.logger)
        self.images[imname] = img

    
    #-----------
    def nextObservation(self):
        """
        Given that an observation can contain a set of offsets and filters, this moves to the next
        observation in that set, and re-initializes the local instrument. If possible, it keeps the
        same calculated PSF.
        """
        self.obs_count += 1
        msg = "Initializing Observation {} of {}"
        self._log("info",  msg.format(self.obs_count, len(self.observations)))
        if self.obs_count < len(self.observations):
            filter = self.observations[self.obs_count]['filter']
            self._log("info","Observation Filter is {}".format(filter))
            offset = self.observations[self.obs_count]['offset']
            offset_ra = float(offset['offset_ra'])/3600.
            offset_dec = float(offset['offset_dec'])/3600.
            offset_pa = float(offset['offset_pa'])
            ra, dec = OffsetPosition(self.ra, self.dec, float(offset_ra), 
                                     float(offset_dec))
            pa = (self.pa + float(offset_pa))%360.
            if self.observations[self.obs_count]['offset']['offset_centre']:
                offset_ra = self.instrument.INSTRUMENT_OFFSET[0]/3600.
                offset_dec = self.instrument.INSTRUMENT_OFFSET[1]/3600.
                offset_pa = self.instrument.INSTRUMENT_OFFSET[2]
                ra,dec = OffsetPosition(ra,dec,offset_ra,offset_dec)
                pa = (pa + offset_pa)%360.
            msg = "Observation (RA,DEC) = ({:.3f},{:.3f}) with PA={:.3f}"
            self._log("info", msg.format(ra,dec,pa))
            self.instrument.reset(ra, dec, pa, filter, self.obs_count)
            self._log("info", "Reset Instrument")
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
        self._log("info", "Running catalogue %s" % (catalogue))
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
        self._log("info", 'Adding image {} to observation'.format(file))
        if self.images[img].header["ASTROIMAGEVALID"]:
            self.instrument.addImage(self.images[img], units, *args, **kwargs)
        self._log("info", 'Image Added')


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
        psf_path = os.path.join(SelectParameter('psf_cache_location'),
                                SelectParameter('psf_cache_directory'))
        if os.path.exists(psf_path):
            psf_names = glob.glob(os.path.join(psf_path, "*.fits"))
        self._log("info", "Adding Error")
        if 'parallel' not in kwargs:
            kwargs['parallel'] = self.parallel
        if 'cores' not in kwargs:
            kwargs['cores'] = self.cores
        self.instrument.addError(residual_poisson=self.residual_poisson, 
                                 residual_readnoise=self.residual_readnoise, 
                                 residual_flat=self.residual_flat, 
                                 residual_dark=self.residual_dark, 
                                 residual_cosmic=self.residual_cosmic, 
                                 *args, **kwargs)
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
        fits_name = "{}_{:03d}.fits".format(self.imgbase, self.obs_count)
        fits_file = os.path.join(self.out_path, fits_name)
        self.instrument.toFits(fits_file, *args, **kwargs)
        mosaics = None
        if mosaic:
            mosaic_name = "{}_{:03d}_mosaic.fits"
            mosaic_name = mosaic_name.format(self.imgbase, self.obs_count)
            mosaic_file = os.path.join(self.out_path, mosaic_name)
            mosaics = self.instrument.toMosaic(mosaic_file, *args, **kwargs)
        return fits_file, mosaics, self.params
    

    #-----------
    def totalObservations(self):
        """
        Return the total number of observations
        """
        return len(self.observations)


    def _log(self, mtype, message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self,'logger'):
            getattr(self.logger, mtype)(message)
        else:
            sys.stderr.write("{}: {}\n".format(mtype, message))

    
    def updateState(self, state):
        if self.set_celery is not None:
            self.set_celery(state)

    
    def getState(self):
        if self.get_celery is not None:
            return self.get_celery()
        return ""

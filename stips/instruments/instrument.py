from __future__ import absolute_import,division
__filetype__ = "base"

#External Modules
import glob, logging, os, shutil, sys, types, uuid

import numpy as np
import synphot as syn
import stsynphot as stsyn

from astropy import units as u
from astropy.io import fits as pyfits
from astropy.table import Table, Column
from functools import wraps

#Local Modules
from ..stellar_module import StarGenerator
from ..astro_image import AstroImage
from ..utilities import GetStipsData
from ..utilities import internet
from ..utilities import OffsetPosition
from ..utilities import read_metadata
from ..utilities import read_table
from ..utilities import SelectParameter
from ..utilities import StipsDataTable

if sys.version_info[0] >= 3:
    import builtins
    from io import StringIO
else:
    import __builtin__
    from cStringIO import StringIO


class Instrument(object):
    """
    The Instrument class represents a virtual base class which will be implemented as a variety of
        JWST, HST, and Roman actual instruments. The Instrument class contains:
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        filter    : string, what filter of the instrument is being observed
        out_path  : place to put temporary files and the like
    """

    def __init__(self, **kwargs):
        """
        Instrument. The __init__ function creates a (potentially) empty instrument.
        """
        self.COMPFILES =  sorted(glob.glob(os.path.join(os.environ["PYSYN_CDBS"],"mtab","*tmc.fits")))
        self.GRAPHFILES = sorted(glob.glob(os.path.join(os.environ["PYSYN_CDBS"],"mtab","*tmg.fits")))
        self.THERMFILES = sorted(glob.glob(os.path.join(os.environ["PYSYN_CDBS"],"mtab","*tmt.fits")))

        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            self.logger = logging.getLogger('__stips__')
            log_level = SelectParameter('log_level', kwargs)
            self.logger.setLevel(getattr(logging, log_level))
            if not len(self.logger.handlers):
                stream_handler = logging.StreamHandler(sys.stderr)
                format = '%(asctime)s %(levelname)s: %(message)s'
                stream_handler.setFormatter(logging.Formatter(format))
                self.logger.addHandler(stream_handler)

        self.out_path = SelectParameter('out_path', kwargs)
        self.prefix = kwargs.get('prefix', '')
        self.cat_type = SelectParameter('cat_type', kwargs)
        self.flatfile = GetStipsData(os.path.join("residual_files", self.FLATFILE))
        self.darkfile = GetStipsData(os.path.join("residual_files", self.DARKFILE))
        self.oversample = SelectParameter('oversample', kwargs)
        self.psf_grid_size = SelectParameter('psf_grid_size', kwargs)
        self.seed = SelectParameter('seed', kwargs)
        self.imgbase = kwargs.get('imgbase', '')
        self.ra = kwargs.get('ra', 0.)
        self.dec = kwargs.get('dec', 0.)
        self.pa = kwargs.get('pa', 0.)
        self.distortion = SelectParameter('distortion', kwargs)
        self.exptime = kwargs.get('exptime', 1.)
        self.small_subarray = kwargs.get('small_subarray', False)
        self.filter = None
        self.detectors = None
        self.psf_commands = kwargs.get('psf_commands', None)
        self.instrument = kwargs.get('instrument', "")
        self.background_value = SelectParameter('background', kwargs)
        self.background_location = SelectParameter('jbt_location', kwargs)
        self.custom_background = kwargs.get('custom_background', 0.)
        self.CENTRAL_OFFSET = (0., 0., 0.)
        self.convolve_size = SelectParameter('convolve_size', kwargs)
        self.memmap = SelectParameter('memmap', kwargs)
        self.set_celery = kwargs.get('set_celery', None)
        self.get_celery = kwargs.get('get_celery', None)

        #Adjust # of detectors based on keyword:
        n_detectors = int(kwargs.get('detectors', len(self.DETECTOR_OFFSETS)))
        self.DETECTOR_OFFSETS = self.DETECTOR_OFFSETS[:n_detectors]
        self.OFFSET_NAMES = self.OFFSET_NAMES[:n_detectors]
        self.CENTRAL_OFFSET = self.N_OFFSET[n_detectors]
        msg = "{} with {} detectors. Central offset {}"
        self._log('info', msg.format(self.DETECTOR, n_detectors, 
                                     self.CENTRAL_OFFSET))
    
    @classmethod
    def initFromImage(cls, image, **kwargs):
        """
            Takes an input AstroImage, and does an add-with-align for every detector present. 
            
            If units are present, does a unit conversion.
        """
        units = kwargs.get('unit', 'c')
        self._log("info","Initializing from with units %s" % (units))
        ins = cls(**kwargs)
        img = image * ins.convertToCounts(unit,scalex=image.scale[0],scaley=image.scale[1])
        self._log("info","Converted image units")
        for detector in ins.detectors:
            self._log("info","Adding image to detector %s" % (detector.name))
            detector.addWithAlignment(img)
        self._log("info","Finished initialization")
        return ins

    @classmethod
    def initFromCatalogue(cls,catalogue,**kwargs):
        """
            Takes an input catalogue, and observes that catalogue with all detectors.
            
            It currently assumes that the input catalogue has the following columns:
                RA: RA of source
                DEC: DEC of source
                FLUX: flux of source
                TYPE: type of source (point, sersic)
                N: sersic index
                Re: radius containing half of the light of the sersic profile
                Phi: angle of the major axis of the sersic profile
                Ratio: axial ratio of the Sersic profile
            Obtaining the correct values for FLUX (if not done before initialization) is a job for
            the subclasses.
        """
        self._log("info","Initializing with catalogue %s" % (catalogue))
        ins = cls(**kwargs)
        self._log("info","Converting catalogue to internal format")
        cat = ins.convertCatalogue(catalogue)
        for detector in ins.detectors:
            self._log("info","Adding image to detector %s" % (detector.name))
            detector.addCatalogue(cat, dist=self.distortion)
        self._log("info","Finished initialization")
        return ins
    
    def reset(self, ra, dec, pa, filter, obs_count, psf=True, detectors=True, celery=None):
        """
        Reset instrument parameters.
        """
        self._log("info","Resetting")
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.obs_count = obs_count
        if filter != self.filter:
            if filter not in self.FILTERS:
                msg = "Filter {} is not a valid {} filter"
                raise ValueError(msg.format(filter, self.instrument))
            self.filter = filter
            self.background = self.pixel_background
            self.photfnu = self.PHOTFNU[self.filter]
            self.photplam = self.PHOTPLAM[self.filter]
            if hasattr(self, "_bp"):
                del self._bp
        if detectors:
            self.resetDetectors(psf=psf)
    

    def resetDetectors(self, psf=True):
        if self.detectors is not None:
            del self.detectors
        #Create Detectors
        self.detectors = []
        for offset, name in zip(self.DETECTOR_OFFSETS, self.OFFSET_NAMES):
            distortion = None
            if self.distortion and hasattr(self, 'DISTORTION'):
                distortion = self.DISTORTION[name]
            (delta_ra, delta_dec, delta_pa) = offset
            delta_ra = (delta_ra - self.CENTRAL_OFFSET[0])/3600.
            delta_dec = (delta_dec - self.CENTRAL_OFFSET[1])/3600.
            delta_pa = delta_pa - self.CENTRAL_OFFSET[2]
            ra,dec = OffsetPosition(self.ra, self.dec, delta_ra, delta_dec)
            pa = (self.pa + delta_pa)%360.
            hdr = {"DETECTOR":name, "FILTER":self.filter}
            msg = "Initialized {} Detector {} with filter {}"
            hist = [msg.format(self.instrument, name, self.filter)]
            msg = "Creating Detector {} with (RA,DEC,PA) = ({},{},{})"
            self._log("info", msg.format(name, ra, dec, pa))
            msg = "Creating Detector {} with offset ({},{})"
            self._log("info", msg.format(name, delta_ra, delta_dec))
            detector = AstroImage(parent=self, ra=ra, dec=dec, pa=pa, psf=psf,
                                  header=hdr, history=hist, detname=name,
                                  distortion=distortion)
            self._log("info", "Detector {} created".format(name))
            self.detectors.append(detector)
        
    def toFits(self,outfile):
        """
        Takes the detectors and turns them into a multi-extension FITS file.
        """
        self._log("info","Converting to FITS file")
        hdus = [pyfits.PrimaryHDU()]
        for detector in self.detectors:
            self._log("info","Converting detector %s to FITS extension" % (detector.name))
            hdus.append(detector.imageHdu)
        hdulist = pyfits.HDUList(hdus)
        hdulist.writeto(outfile, overwrite=True)
        self._log("info","Created FITS file %s" % (outfile))
    
    def toMosaic(self,outfile):
        """
        Creates a single FITS file from each detector, and then uses Montage to create a mosaic
        of all of the said files.
        """
        if len(self.detectors) > 1:
            import montage_wrapper as montage
            self._log("info","Converting to FITS mosaic")
            tmp_dir = os.path.join(self.out_path,"tmp-"+str(uuid.uuid4()))
            os.makedirs(tmp_dir)
            tmp_out_dir = os.path.join(self.out_path,"tmp-out-"+str(uuid.uuid4()))
            tmp_work_dir = os.path.join(self.out_path,"tmp-work-"+str(uuid.uuid4()))
            self.toFits(os.path.join(tmp_dir,"image.fits"))
            montage.mosaic(tmp_dir,tmp_out_dir,background_match=True,work_dir=tmp_work_dir)
            self._log("info","Mosaic finished running")
            shutil.copy(os.path.join(tmp_out_dir,"mosaic.fits"),outfile)
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            if os.path.exists(tmp_out_dir):
                shutil.rmtree(tmp_out_dir)
            if os.path.exists(tmp_work_dir):
                shutil.rmtree(tmp_work_dir)
            self._log("info","Created FITS file %s and cleaned up" % (outfile))
            return [outfile]
        else:
            self._log("info","Not creating single-detector mosaic")
            return []

    def addImage(self, image, unit='c'):
        """
            Takes an input AstroImage, and does an add-with-align for every detector present. 
            
            If units are present, does a unit conversion.
        """
        self._log("info","Adding image with units %s" % (unit))
        img = image * self.convertToCounts(unit,scalex=image.scale[0],scaley=image.scale[1])
        self._log("info","Converted image count rate")
        for detector in self.detectors:
            self._log("info","Adding image to detector %s" % (detector.name))
            detector.addWithAlignment(img)
        self._log("info","Finished adding image")
    
    def addCatalogue(self, catalogue, obs_num, *args, **kwargs):
        """
            Takes an input catalogue, and observes that catalogue with all detectors.
            
            It currently assumes that the input catalogue has the following columns:
                RA: RA of source
                DEC: DEC of source
                FLUX: flux of source
                TYPE: type of source (point, sersic)
                N: sersic index
                Re: radius containing half of the light of the sersic profile
                Phi: angle of the major axis of the sersic profile
                Ratio: axial ratio of the Sersic profile
            Obtaining the correct values for FLUX (if not done before initialization) is a job for
            the subclasses.
        """
        base_state = self.getState()
        self._log("info","Adding catalogue %s" % (catalogue))
        self.updateState(base_state + "<br /><span class='indented'>Converting Catalogue to Internal Format</span>")
        cat = self.convertCatalogue(catalogue, obs_num)
        self._log("info","Finished converting catalogue to internal format")
        cats = [cat]
        for detector in self.detectors:
            self.updateState(base_state + "<br /><span class='indented'>Detector {}</span>".format(detector.name))
            self._log("info","Adding catalogue to detector %s" % (detector.name))
            cats.append(detector.addCatalogue(cat, dist=self.distortion, *args, **kwargs))
            self.updateState(base_state)
        return cats
        self._log("info","Finished Adding Catalogue")
    
    def addTable(self, table, table_type, *args, **kwargs):
        """
            Takes an input table (still in memory), and observes that table with all detectors.
            
            It assumes that the table type is a supported type.

            Obtaining the correct values for FLUX (if not done before initialization) is a job for
            the subclasses.
        """
        self._log("info","Adding {} table".format(table_type))
        conversion_fn = self.getTableFunction(table_type)
        internal_table, cache = conversion_fn(table, self.bandpass)
        self._log("info","Finished converting table to internal format")
        tables = [internal_table]
        for detector in self.detectors:
            self._log("info","Adding table to detector %s" % (detector.name))
            tables.append(detector.addTable(internal_table, dist=self.distortion))
        self._log("info","Finished Adding Catalogue")        
        return tables
    
    def convertToCounts(self,unit,scalex=None,scaley=None):
        """
        Convert input to Counts.
    
        unit: one of 'p' (photons/s), 'e' (erg/s), 'c' (counts/s), 'j' (Jansky), 's' (W/m/m^2/Sr)
        
        scale: needed for surface brightness conversions. Arcseconds/pixel
    
        returns: factor (multiplicative conversion factor)
        """
        units = ('p','e','c','j','s')
        if unit not in units: raise ValueError("Unit %s is not one of %s" % (unit,units))
        
        if unit == 'c':
            return 1.
        elif unit == 'j':
            return 1./self.photfnu
        elif unit =='p':
            freq = 299792458000000. / self.photplam #c in um/s
            energy = 6.6260755e-27 * freq # h in erg*s to give energy in ergs
            return 1.e23 * energy / (self.photfnu * self.AREA * freq)
        elif unit == 'e':
            freq = 299792458000000. / self.photplam #c in um/s
            return 1.e23 / (self.photfnu * self.AREA * freq)
        else: #unit == 's'
            #W/m/m^2/Sr -> Jy = 1.e14 * photplam^2 * (scalex*scaley/206265.^2) / 3.e8
            # Jy -> counts = 1./photfnu
            # Combined as follows
            return 1.e14 * self.photplam**2 * (scalex*scaley/42545250225.) / (3.e8 * self.photfnu)

    def convertCatalogue(self, catalogue, obs_num):
        """
        Converts a catalogue to the expected format for AstroImage, including doing unit conversions
        of columns if necessary. Acceptable formats are:
            - Phoenix (models from the Phoenix stellar grid)
            - BC95 (galaxy models from BC95)
            - Internal (has columns RA/DEC/FLUX/TYPE/N/Re/Phi/Ratio/ID/Notes)
            - Mixed (has columns RA/DEC/FLUX/UNITS/TYPE/N/Re/Phi/Ratio/ID/Notes)
            - Generic (has columns RA/DEC/FILTER where FILTER == self.filter (and possibly also
                       other filters)
        
        catalogue: catalogue name of input catalogue
        
        returns: cat: new catalogue in Internal format
        """
        (in_cat_path, in_cat_name) = os.path.split(catalogue)
        (in_cat_base, ext) = os.path.splitext(in_cat_name)
        obs_cat_name = "{}_{:02d}_conv_{}.{}".format(in_cat_base, obs_num, self.filter, self.cat_type)
        obsname = os.path.join(self.out_path, obs_cat_name)
        in_data_table = StipsDataTable.dataTableFromFile(catalogue)
        cols = in_data_table.columns
        if "keywords" in in_data_table.meta:
            meta = {k.lower():v['value'] for k,v in in_data_table.meta['keywords'].items()}
        else:
            meta = {k.lower():v for k,v in in_data_table.meta.items()}
        #Check for built-in metadata
        table_type = ""
#         if 'keywords' in in_meta:
#             if 'type' in in_meta['keywords']:
#                 table_type = in_meta['keywords']['type']['value']
        if 'type' in meta:
            table_type = meta['type']
        if table_type in ['phoenix', 'phoenix_realtime', 'bc95']:
            pass
        elif table_type == 'internal':
            filter = meta['filter'].lower()
#             filter = t.meta['keywords']['filter']['value'].lower()
            if filter != self.filter.lower(): 
                raise ValueError("Adding catalogue with filter {} to {} {}".format(filter, self.DETECTOR, self.filter))
            return catalogue
        elif table_type == 'mixed':
            filter = meta['filter'].lower()
#             filter = t.meta['keywords']['filter']['value'].lower()
            if filter != self.filter.lower(): 
                raise ValueError("Adding catalogue with filter {} to {} {}".format(filter, self.DETECTOR, self.filter))
        elif table_type == 'multifilter':
            if self.filter.lower() not in [c.lower() for c in in_data_table.columns]:
                raise ValueError("Adding catalogue with filters {} to {} {}".format(in_data_table.columns, self.DETECTOR, self.filter))
        else: #check for necessary columns
            #We want RA, DEC, and count rate in the appropriate filter
            if 'ra' not in cols or 'dec' not in cols or self.filter.lower() not in [c.lower() for c in cols]:
                raise ValueError("Can't parse catalogue without proper columns")
        return self.handleConversion(catalogue, table_type, obsname)
    
    def getTableFunction(self, table_type):
        if table_type == 'phoenix':
            return self.readPhoenixTable
        elif table_type == 'phoenix_realtime':
            return self.readPhoenixRealtimeTable
        elif table_type == 'pandeia':
            return self.readPandeiaTable
        elif table_type == 'bc95':
            return self.readBC95Table
        elif table_type == 'mixed':
            return self.readMixedTable
        elif table_type == 'multifilter':
            return self.readMultiTable
        return self.readGenericTable
    
    def getTableFormat(self, table_type):
        if table_type == 'phoenix':
            return {'ra': ' %10g', 'dec': ' %10g', 'flux': ' %12g', 'type': '%6s', 'n': '%4s', 're': '%4s', 'phi': '%4s', 'ratio': '%6s', 'id': '%8d', 'notes': '%-25s'}
        elif table_type == 'phoenix_realtime':
            return {'ra': ' %10g', 'dec': ' %10g', 'flux': ' %12g', 'type': '%6s', 'n': '%4s', 're': '%4s', 'phi': '%4s', 'ratio': '%6s', 'id': '%8d', 'notes': '%-25s'}
        elif table_type == 'pandeia':
            return {'ra': ' %10g', 'dec': ' %10g', 'flux': ' %12g', 'type': '%6s', 'n': '%4s', 're': '%4s', 'phi': '%4s', 'ratio': '%6s', 'id': '%8d', 'notes': '%-25s'}
        elif table_type == 'bc95':
            return {'ra': ' %10g', 'dec': ' %10g', 'flux': ' %12g', 'type': '%6s', 'n': '%6.3g', 're': '%10g', 'phi': ' %10g', 'ratio': '%10g', 'id': '%8d', 'notes': '%-25s'}
        elif table_type == 'mixed':
            return {}
        elif table_type == 'multifilter':
            return {}
        return {}
    
    def handleConversion(self, catalogue, table_type, obsname):
        """
        Converts an input catalogue into an internal format catalogue. Needs the appropriate 
        function (e.g. readPhoenixTable) in order to call it.
        """
        self._log("info","Converting {} catalogue".format(table_type))
        self._log("info", "Preparing output table")
        if os.path.isfile(obsname): 
            os.remove(obsname)
        cat_function = self.getTableFunction(table_type)
        
        in_data_table = StipsDataTable.dataTableFromFile(catalogue)
        out_data_table = StipsDataTable.dataTableFromFile(obsname)
        out_data_table.meta = {'name': 'Internal Format Catalogue', 'type': 'internal', 
                               'filter': self.filter}
        
        
        bp = self.bandpass
        cached = -1
        current_chunk = in_data_table.read_chunk()
        while current_chunk is not None:
            self._log("info", "Converting chunk {}".format(in_data_table.chunk))
            out_chunk, cached = cat_function(current_chunk, bp, cached)
            out_data_table.write_chunk(out_chunk)
            current_chunk = in_data_table.read_chunk()
        return obsname
    
    def readPhoenixTable(self, table, bp, cached=-1):
        """
        Takes a table (or fraction of a table) with the data in the Phoenix catalogue format, and 
        return an output table (not yet in ascii form) in the internal format, with those sources
        in it.
        """
        if table is None:
            return None
        self._log("info","Converting Phoenix Table to Internal format")
        ids = table['id']
        datasets = table['dataset']
        ras = table['ra']
        decs = table['dec']
        ages = table['age']
        metallicities = table['metallicity']
        masses = table['mass']
        distances = table['distance']
        binaries = table['binary']
        rates = np.zeros_like(ras)
        all_datasets = np.unique(datasets)
        self._log("info","{} datasets".format(len(all_datasets)))
        for dataset in all_datasets:
            idx = np.where(datasets == dataset)
            if len(idx[0]) > 0:
                stargen = StarGenerator(ages[idx][0], metallicities[idx][0], seed=self.seed, logger=self.logger)
                rates[idx] = stargen.make_cluster_rates(masses[idx], self, bp)
                del stargen
        rates = rates * 100 / (distances**2) #convert absolute to apparent rates
        if cached > 0:
            rates[0] += cached
            cached = -1
        #Now, deal with binaries. Remember that if Binary == 1, then the star below is the binary companion.
        idx = np.where(binaries==1.0)[0] #Stars with binary companions
        idxp = (idx+1) #binary companions
        if len(idxp) > 0 and idxp[-1] >= len(rates): #last one is a binary. Cache it.
            cached = rates[-1]
            idx, idxp = idx[:-1], idxp[:-1]
            ids, datasets, ras, decs, ages, metallicities, = ids[:-1], datasets[:-1], ras[:-1], decs[:-1], ages[:-1], metallicities[:-1]
            masses, distances, binaries, rates = masses[:-1], distances[:-1], binaries[:-1], rates[:-1]
        rates[idx] += rates[idxp] #add count rates together
        # Now that we've added rates together, remove the binary companions
        ras = np.delete(ras,idxp)
        decs = np.delete(decs,idxp)
        rates = np.delete(rates,idxp)
        ids = np.delete(ids,idxp)
        binaries = np.delete(binaries,idxp)
        notes = np.empty_like(ras,dtype="S6")
        notes[np.where(binaries==1)] = 'Binary'
        notes[np.where(binaries!=1)] = 'None'
        t = Table()
        t['ra'] = Column(data=ras)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=np.full_like(ras,"point",dtype="S6"))
        t['n'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['re'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['phi'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['ratio'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=notes)
        
        return t, cached

    def readPhoenixRealtimeTable(self, table, bp, cached=-1):
        """
        Converts a set of Phoenix sources specified as (id, ra, dec, T, Z, log(g), apparent) into individual stars, observes them,
        and then produces an output catalogue
        """
        if table is None:
            return None
        ids = table['id']
        ras = table['ra']
        decs = table['dec']
        temps = table['teff']
        gravs = table['log_g']
        metallicities = table['metallicity']
        apparents = table['apparent']
        norm_bp = '{}'.format(apparents.unit)
        if norm_bp == '':
            norm_bp = 'johnson,i'
        rates = np.zeros_like(ras)
        for index in range(len(ids)):
            t, g, Z, a = temps[index], gravs[index], metallicities[index], apparents[index]
            sp = stsyn.grid_to_spec('phoenix', t, Z, g)
            sp = self.normalize(sp, a, norm_bp)
            obs = syn.Observation(sp, bp, binset=sp.waveset)
            rates[index] = obs.countrate(area=self.AREA)
        t = Table()
        t['ra'] = Column(data=ras)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=np.full_like(ras,"point",dtype="S6"))
        t['n'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['re'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['phi'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['ratio'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=np.full_like(ras,"None",dtype="S6"))
        
        return t, cached
    
    def readPandeiaTable(self, table, bp, cached=-1):
        """
        Converts a set of Pandeia phoenix sources specified as (id, ra, dec, key, apparent) into individual stars, observes them,
        and then produces an output catalogue
        """
        if table is None:
            return None
        from pandeia.engine.sed import SEDFactory
        ids = table['id']
        ras = table['ra']
        decs = table['dec']
        keys = table['key']
        apparents = table['apparent']
        norm_bp = '{}'.format(apparents.unit)
        if norm_bp == '':
            norm_bp = 'johnson,i'
        rates = np.array((), dtype='float32')
        for a, key in zip(apparents, keys):
            config = {'sed_type': 'phoenix', 'key': key}
            spectrum = SEDFactory(config=config)
            wave, flux = spectrum.get_spectrum()
            sp = self.normalize((wave, flux), a, norm_bp)
            obs = syn.Observation(sp, bp, binset=sp.wave)
            rates = np.append(rates, obs.countrate(area=self.AREA))
        t = Table()
        t['ra'] = Column(data=ras)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=np.full_like(ras,"point",dtype="S6"))
        t['n'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['re'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['phi'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['ratio'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=np.full_like(ras,"None",dtype="S6"))
        
        return t, cached
    
    def readBC95Table(self, table, bp, cached=-1):
        """
        Converts a BC95 galaxy grid of sources into the internal source table standard
        """
        from pandeia.engine.custom_exceptions import SynphotError as pSynError
        if table is None:
            return None
        # This function is needed because I can't get python not to read '50E8' as a number, or to output it as a correctly formatted string
        def stringify(num):
            num = float(num)
            exponent = int(np.floor(np.log10(num)))
            value = int(10*(num/(10**np.floor(np.log10(num)))))
            return "{}E{}".format(value, exponent-1)
        
        self._log("info", "Converting BC95 Catalogue")
        proflist = {"expdisk":1,"devauc":4}
        distance_type = "redshift"
        ids = table['id']
        ras = table['ra']
        decs = table['dec']
        try:
            zs = table['redshift']
        except KeyError: #distances instead
            zs = table['distance']
            distance_type = "pc"
        models = table['model']
        ages = table['age']
        profiles = table['profile']
        radii = table['radius']/self.SCALE[0] #at some point, may want to figure out multiple scales.
        ratios = table['axial_ratio']
        pas = (table['pa'] + (self.pa*180./np.pi) )%360.
        vmags = table['apparent_surface_brightness']
        norm_bp = '{}'.format(vmags.unit)
        self._log("info", "Normalization Bandpass is {} ({})".format(norm_bp, type(norm_bp)))
        if norm_bp == '' or norm_bp is None or norm_bp == 'None':
            norm_bp = 'johnson,v'
        self._log("info", "Normalization Bandpass is {}".format(norm_bp))
        rates = np.array(())
        indices = np.array(())
        notes = np.array((),dtype='object')
        total = len(models)
        for i, (z,model,age,profile,radius,ratio,pa,mag) in enumerate(zip(zs,models,ages,profiles,radii,ratios,pas,vmags)):
#             self._log("info", "{} of {}: {} {} {} {} {} {} {} {}".format(i, total, z, model, age, profile, radius, ratio, pa, mag))
            fname = "bc95_{}_{}.fits".format(model, stringify(age))
            try:
                sp = syn.SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'],"grid","bc95","templates",fname))
                if distance_type == "redshift":
                    sp = syn.SourceSpectrum(sp.model, z=z)
                sp = self.normalize(sp, mag, norm_bp)
                obs = syn.Observation(sp, self.bandpass, force='taper', 
                                      binset=sp.waveset)
                rate = obs.countrate(area=self.AREA).value
            except pSynError as e:
                msg = 'Source {} of {}: Pysynphot Error {} encountered'
                self._log('warning', msg.format(i, total, e))
                rate = 0.
            rates = np.append(rates, rate)
            indices = np.append(indices, proflist[profile])
            note = "BC95_{}_{}_{}".format(model, stringify(age), mag)
            notes = np.append(notes, note)
        t = Table()
        t['ra'] = Column(data=ras, dtype=np.float)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=np.full_like(ras, 'sersic', dtype='S7'))
        t['n'] = Column(data=indices)
        t['re'] = Column(data=radii)
        t['phi'] = Column(data=pas)
        t['ratio'] = Column(data=ratios)
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=notes, dtype='S25')
        
        return t, cached

    def readMixedTable(self, table, bp, cached=-1):
        """
        Converts a mixed internal list of sources into the internal source table standard
        """
        if table is None:
            return None
        ras = table['ra']
        decs = table['dec']
        types = table['type']
        indices = table['n']
        radii = table['re']
        pas = table['phi']
        ratios = table['ratio']
        ids = table['id']
        notes = table['notes']

        rates = table['flux']
        units = table['units']
        idxp = np.where(units == 'p')
        rates[idxp] *= convertToCounts('p')
        idxe = np.where(units == 'e')
        rates[idxe] *= convertToCounts('e')
        idxj = np.where(units == 'j')
        rates[idxj] *= convertToCounts('j')
        t = Table()
        t['ra'] = Column(data=ras)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=types)
        t['n'] = Column(data=indices)
        t['re'] = Column(data=radii)
        t['phi'] = Column(data=pas)
        t['ratio'] = Column(data=ratios)
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=notes)
        
        return t, cached
    
    def readMultiTable(self, table, bp, cached=-1):
        """
        Converts an internal multifilter list of sources into the internal source table standard
        """
        if table is None:
            return None
        ras = table['ra']
        decs = table['dec']
        types = table['type']
        indices = table['n']
        radii = table['re']/self.SCALE[0] #at some point, may want to figure out multiple scales.
        pas = table['phi']
        ratios = table['ratio']
        ids = table['id']
        notes = table['notes']
        rates = table[self.filter]
        t = Table()
        t['ra'] = Column(data=ras)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=types)
        t['n'] = Column(data=indices)
        t['re'] = Column(data=radii)
        t['phi'] = Column(data=pas)
        t['ratio'] = Column(data=ratios)
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=notes)
        
        return t, cached

    def readGenericTable(self, table, bp, cached=-1):
        """
        Converts a generic list of point sources into the internal source table standard
        """
        if table is None:
            return None
        ras = table['ra']
        decs = table['dec']
        rates = table[self.filter.lower()]
        if 'id' in table:
            ids = table['id']
        else:
            ids = np.arange(len(ras),dtype=int)
        t = Table()
        t['ra'] = Column(data=ras)
        t['dec'] = Column(data=decs)
        t['flux'] = Column(data=rates)
        t['type'] = Column(data=np.full_like(ras,"point",dtype="S6"))
        t['n'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['re'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['phi'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['ratio'] = Column(data=np.full_like(ras,"N/A",dtype="S3"))
        t['id'] = Column(data=ids)
        t['notes'] = Column(data=np.full_like(ras, "N/A", dtype="S3"))
        
        return t, cached

    def generateReadnoise(self):
        """Base function for adding read noise"""
        pass

    @classmethod
    def handleDithers(cls,form):
        """Base function for handling dither patterns"""
        pass

    @classmethod
    def doSubpixel(cls,dithers,subpixels):
        """For each (x,y) in dithers, dither around that point for each point in subpixels. Return the full set"""
        my_dithers = []
        for (x,y) in dithers:
            for (i,j) in subpixels:
                my_dithers.append((x+i,y+j))
        return my_dithers
    
    def addError(self, *args, **kwargs):
        """Base function for adding in residual error"""
        self._log("info","Adding residual error")
        cores = SelectParameter('cores', kwargs)
        convolve = SelectParameter('convolve', kwargs)
        poisson = SelectParameter('residual_poisson', kwargs)
        readnoise = SelectParameter('residual_readnoise', kwargs)
        flat = SelectParameter('residual_flat', kwargs)
        dark = SelectParameter('residual_dark', kwargs)
        cosmic = SelectParameter('residual_cosmic', kwargs)
        parallel = SelectParameter('parallel_enable', kwargs)
        snapshots = kwargs.get("snapshots", {})

        base_state = self.getState()
        if flat:
            flat = AstroImage.initDataFromFits(self.flatfile,ext='COMPRESSED_IMAGE', psf=False, logger=self.logger)
        if dark:
            dark = AstroImage.initDataFromFits(self.darkfile,ext='COMPRESSED_IMAGE', psf=False, logger=self.logger)
            dark *= self.exptime
        if readnoise:
            rn = self.generateReadnoise()
        for detector in self.detectors:
            if 'initial' in snapshots:
                detector.toFits(self.imgbase+"_{}_{}_snapshot_initial.fits".format(self.obs_count, detector.name))
            self.updateState(base_state + "<br /><span class='indented'>Detector {}: Adding Background</span>".format(detector.name))
            self._log("info","Adding error to detector %s" % (detector.name))
            self._log("info","Adding background")
            self._log("info","Background is {} counts/s/pixel".format(self.pixel_background))
            detector.addBackground(self.pixel_background)
            if 'background' in snapshots or 'all' in snapshots:
                detector.toFits(self.imgbase+"_{}_{}_snapshot_background.fits".format(self.obs_count, detector.name))
            self._log("info","Inserting correct exposure time")
            self.updateState(base_state + "<br /><span class='indented'>Detector {}: Applying Exposure Time</span>".format(detector.name))
            detector.setExptime(self.exptime)
            if 'exptime' in snapshots or 'all' in snapshots:
                detector.toFits(self.imgbase+"_{}_{}_snapshot_exptime.fits".format(self.obs_count, detector.name))
            self._log("info","Convolving with PSF")
            convolve_state = base_state + "<br /><span class='indented'>Detector {}: Convolving PSF</span>".format(detector.name)
            self.updateState(convolve_state)
            detector.convolve_psf(max_size=self.convolve_size-1, parallel=parallel, cores=cores)
            if 'convolve' in snapshots or 'all' in snapshots:
                detector.toFits(self.imgbase+"_{}_{}_snapshot_convolve.fits".format(self.obs_count, detector.name))
            if self.oversample != 1:
                self._log("info","Binning oversampled image")
                self.updateState(base_state + "<br /><span class='indented'>Detector {}: Binning oversampled image</span>".format(detector.name))
                detector.bin(self.oversample)
                if 'bin' in snapshots or 'all' in snapshots:
                    detector.toFits(self.imgbase+"_{}_{}_snapshot_bin.fits".format(self.obs_count, detector.name))
            if poisson: 
                self._log("info","Adding poisson noise")
                self.updateState(base_state + "<br /><span class='indented'>Detector {}: Adding Poisson Noise</span>".format(detector.name))
                detector.introducePoissonNoise()
                if 'poisson' in snapshots or 'all' in snapshots:
                    detector.toFits(self.imgbase+"_{}_{}_snapshot_poisson.fits".format(self.obs_count, detector.name))
            if readnoise:
                self._log("info","Adding readnoise")
                self.updateState(base_state + "<br /><span class='indented'>Detector {}: Adding Readnoise</span>".format(detector.name))
                detector.introduceReadnoise(rn)
                if 'readnoise' in snapshots or 'all' in snapshots:
                    detector.toFits(self.imgbase+"_{}_{}_snapshot_readnoise.fits".format(self.obs_count, detector.name))
            if flat:
                self._log("info","Adding flatfield residual")
                self.updateState(base_state + "<br /><span class='indented'>Detector {}: Adding Flatfield Residual</span>".format(detector.name))
                detector.introduceFlatfieldResidual(flat)
                if 'flat' in snapshots or 'all' in snapshots:
                    detector.toFits(self.imgbase+"_{}_{}_snapshot_flat.fits".format(self.obs_count, detector.name))
            if dark:
                self._log("info","Adding dark residual")
                self.updateState(base_state + "<br /><span class='indented'>Detector {}: Adding Dark Residual</span>".format(detector.name))
                detector.introduceDarkResidual(dark)
                if 'dark' in snapshots or 'all' in snapshots:
                    detector.toFits(self.imgbase+"_{}_{}_snapshot_dark.fits".format(self.obs_count, detector.name))
            if cosmic:
                self._log("info","Adding cosmic ray residual")
                self.updateState(base_state + "<br /><span class='indented'>Detector {}: Adding Cosmic Ray Residual</span>".format(detector.name))
                detector.introduceCosmicRayResidual(self.PIXEL_SIZE)
                if 'cr' in snapshots or 'all' in snapshots:
                    detector.toFits(self.imgbase+"_{}_{}_snapshot_cr.fits".format(self.obs_count, detector.name))
            self.updateState(base_state)
        self._log("info","Finished adding error")

    def normalize(self, source_spectrum_or_wave_flux, norm_flux, bandpass):
        from pandeia.engine.normalization import NormalizationFactory
        norm_type = self.get_type(bandpass)
        
        norm = NormalizationFactory(type=norm_type, bandpass=bandpass, norm_fluxunit='abmag', norm_flux=norm_flux)
        if isinstance(source_spectrum_or_wave_flux, tuple):
            wave, flux = source_spectrum_or_wave_flux
        else:
            wave = source_spectrum_or_wave_flux.waveset
            flux = source_spectrum_or_wave_flux(wave)
        norm_wave, norm_flux = norm.normalize(wave, flux)
        sp = syn.SourceSpectrum(syn.Empirical1D, points=norm_wave, lookup_table=norm_flux)
        return sp
    
    def get_type(self, bandpass_str):
        if 'miri' in bandpass_str or 'nircam' in bandpass_str:
            return 'jwst'
        #**WFIRST_REMNANT**
        elif 'wfi' in bandpass_str or 'wfirst' in bandpass_str or 'roman' in bandpass_str:
            return 'roman'
        elif 'wfc3' in bandpass_str:
            return 'hst'
        return 'photsys'

    @property
    def bandpass(self):
        if hasattr(self, "_bp"):
            return self._bp
        i = self.pandeia_instrument
        det_params = i.get_detector_pars()
        # 'rn_fudge': multiplied in to match IDT results.
        # 'var_fudge': chromatic fudge factor. quantum yield squared.
        # 'fullwell': 
        # 'ff_electrons':
        # 'pix_size':
        # 
        wr = i.get_wave_range()
        wave = np.linspace(wr['wmin'], wr['wmax'], num=500)
        pce = i.get_total_eff(wave)
        
        if pce[0] != 0.:
            wave = np.insert(wave, 0, wave[0]-(wave[1]-wave[0]))
            pce = np.insert(pce, 0, 0.)
        if pce[-1] != 0.:
            wave = np.append(wave, wave[-1]+(wave[-1]-wave[-2]))
            pce = np.append(pce, 0.)
        
        self._bp = syn.SpectralElement(syn.Empirical1D, points=wave*u.micron, 
                                       lookup_table=pce)
        return self._bp
    
    @property
    def pandeia_instrument(self):
        if hasattr(self, "_instrument"):
            return self._instrument
        from pandeia.engine.calc_utils import build_default_calc
        from pandeia.engine.instrument_factory import InstrumentFactory
    
        translate_instrument = {
                                    #**WFIRST_REMNANT**
                                    'wfi': 'wfirstimager',
                                    'wfi': 'romanimager',
                                    'nircamlong': 'nircam',
                                    'nircamshort': 'nircam',
                                    'miri': 'miri'
                                }

        conf = build_default_calc(self.TELESCOPE.lower(), translate_instrument.get(self.INSTRUMENT.lower(), self.INSTRUMENT.lower()), self.MODE)['configuration']
        conf['instrument']['filter'] = self.filter.lower()
        
        self.logger.info("Creating Instrument with Configuration {}".format(conf['instrument']))
        
        self._instrument = InstrumentFactory(config=conf)
        return self._instrument        
    
    @property
    def zeropoint(self):
        sp = syn.SourceSpectrum.from_vega()
        bp = self.bandpass
        obs = syn.Observation(sp, bp, binset=sp.waveset)
        zeropoint = obs.effstim(flux_unit=u.ABmag).value
        return zeropoint
    
    @property
    def photflam(self):
        return self.photflam_unit.value
    
    @property
    def photflam_unit(self):
        sp = syn.SourceSpectrum(syn.ConstFlux1D, amplitude=(0.*u.STmag))
        bp = self.bandpass
        obs = syn.Observation(sp, bp, binset=sp.waveset)
        pf = (obs.effstim(flux_unit='flam') / obs.countrate(area=self.AREA))
        return pf
    
    @property
    def pixel_background(self):
        return self.pixel_background_unit.value
    
    @property
    def pixel_background_unit(self):
        if isinstance(self.background_value, (int, float)):
            return self.background_value
        elif self.background_value == 'none':
            self._log("info", "Returning background 0.0 for 'none'")
            return 0.
        elif self.background_value == 'custom':
            self._log("info", "Returning background {} for 'custom'".format(self.custom_background))
            return self.custom_background

        bg = None
        if internet() and self.background_location == '$WEB':
            from jwst_backgrounds import jbt
            try:
                bg = jbt.background(self.ra, self.dec, self.PHOTPLAM[self.filter])
            except Exception as e:
                self._log("error", "Accessing JBT background produced error {}".format(e))
                self._log("warning", "Unable to connect to the JBT server")
        
        if bg is None:
            self._log("info", "Using local cache of JBT background data")
            from ..utilities import CachedJbtBackground
            try:
                bg = CachedJbtBackground(self.ra, self.dec, self.PHOTPLAM[self.filter])
            except Exception as e:
                self._log("error", "Retrieving local cache produced error {}".format(e))
                self._log("info", "More complete error: {}".format(repr(e)))
                message = "Unable to retrieve local cache. Returning background 0.0 for '{}'"
                self._log("warning", message.format(self.background_value))
                return 0.

        wave_array = bg.bkg_data['wave_array']
        combined_bg_array = bg.bkg_data['total_bg']
    
        if self.background_value == 'avg':
            flux_array = np.mean(combined_bg_array, axis=0)
        elif self.background_value == 'med':
            flux_array = np.median(combined_bg_array, axis=0)
        elif self.background_value == 'max':
            flux_array = np.max(combined_bg_array, axis=0)
        elif self.background_value == 'min':
            flux_array = np.min(combined_bg_array, axis=0)
        else:
            flux_array = combined_bg_array[0]
    
        # Convert background flux from MJy/sr to mJy/pixel.
        #   Conversion: * 1e9 for MJy -> mJy
        #   Conversion: * 2.3504e-11 for sr^-2 -> arcsec^-2
        #   Conversion: * self.SCALE[0] * self.SCALE[1] for arcsec^-2 -> pixel^-2
        flux_array_pixels = 1e9 * flux_array * 2.3504e-11 * self.SCALE[0] * self.SCALE[1]
    
        sp = syn.SourceSpectrum(syn.Empirical1D, points=wave_array*u.micron, 
                                lookup_table=flux_array_pixels)
        obs = syn.Observation(sp, self.bandpass, binset=sp.waveset, 
                              force='taper')
        bg = obs.countrate(area=self.AREA)
        
        self._log("info", "Returning background {} for '{}'".format(bg, self.background_value))
        return bg

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


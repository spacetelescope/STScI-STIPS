from __future__ import absolute_import,division
__filetype__ = "base"

#External Modules
import glob, logging, numpy, os, shutil, sys, uuid

import pysynphot as ps

from astropy.io import fits as pyfits
from astropy.table import Table,Column
from cStringIO import StringIO
from pandeia.engine.instrument_factory import InstrumentFactory
import montage_wrapper as montage

#Local Modules
from ..stellar_module import StarGenerator
from ..astro_image import AstroImage
from ..utilities import datadir, OffsetPosition, read_metadata, read_table

class Instrument(object):
    """
    The Instrument class represents a virtual base class which will be implemented as a variety of
        JWST, HST, and WFIRST actual instruments. The Instrument class contains:
        
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
        self.in_path = os.environ.get("stips_data", datadir())
        self.out_path = kwargs.get('out_path', os.getcwd())
        self.flatfile = os.path.join(self.in_path, "residual_files", self.FLATFILE)
        self.darkfile = os.path.join(self.in_path, "residual_files", self.DARKFILE)
        self.oversample = 1
        self.ra = kwargs.get('ra', 0.)
        self.dec = kwargs.get('dec', 0.)
        self.pa = kwargs.get('pa', 0.)
        self.distortion = kwargs.get('distortion', False)
        self.filter = None
        self.detectors = None
        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            stream_handler = logging.StreamHandler(sys.stderr)
            stream_handler.setLevel(logging.DEBUG)
            stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
            self.logger = logging.getLogger()
            self.logger.addHandler(stream_handler)
        self.psf_commands = kwargs.get('psf_commands', None)
        self.instrument = kwargs.get('instrument', "")
        self.background_value = kwargs.get('background', 'none')
        self.CENTRAL_OFFSET = (0., 0., 0.)
    
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
    
    def reset(self,ra,dec,pa,filter):
        """
        Reset instrument parameters.
        """
        self._log("info","Resetting")
        self.ra = ra
        self.dec = dec
        self.pa = pa
        if filter != self.filter:
            self.filter = filter
            self.resetPSF()
            self.background = self.BACKGROUND[self.background_value][self.filter]/(self.oversample*self.oversample)
            self.photfnu = self.PHOTFNU[self.filter]
            self.photplam = self.PHOTPLAM[self.filter]
        self.resetDetectors()
    
    def resetPSF(self):
        pass
    
    def resetDetectors(self):
        if self.detectors is not None:
            del self.detectors
        #Create Detectors
        self.detectors = []
        for offset,name in zip(self.DETECTOR_OFFSETS,self.OFFSET_NAMES):
            distortion = None
            if self.distortion:
                distortion = self.DISTORTION[name]
            (delta_ra,delta_dec,delta_pa) = offset
            delta_ra -= self.CENTRAL_OFFSET[0]
            delta_dec -= self.CENTRAL_OFFSET[1]
            delta_pa -= self.CENTRAL_OFFSET[2]
            ra,dec = OffsetPosition(self.ra,self.dec,delta_ra/3600.,delta_dec/3600.)
            pa = (self.pa + delta_pa)%360.
            hdr = {"DETECTOR":name,"FILTER":self.filter}
            hist = ["Initialized %s Detector %s, filter %s" % (self.instrument,name,self.filter)]
            xsize = self.DETECTOR_SIZE[0]*self.oversample
            ysize = self.DETECTOR_SIZE[1]*self.oversample
            scale = [self.SCALE[0]/self.oversample,self.SCALE[1]/self.oversample]
            self._log("info","Creating Detector with (RA,DEC,PA) = (%f,%f,%f)" % (ra,dec,pa))
            self._log("info","Creating Detector with pixel offset ({},{})".format(delta_ra/scale[0], delta_dec/scale[1]))
            detector = AstroImage(out_path=self.out_path, shape=(ysize, xsize), scale=scale, ra=ra, dec=dec, pa=pa, 
                                  header=hdr, history=hist, psf_shape=self.psf.shape, zeropoint=self.zeropoint, 
                                  photflam=self.photflam, detname=name, logger=self.logger, oversample=self.oversample,
                                  distortion=distortion)
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
        hdulist.writeto(outfile, clobber=True)
        self._log("info","Created FITS file %s" % (outfile))
    
    def toMosaic(self,outfile):
        """
        Creates a single FITS file from each detector, and then uses Montage to create a mosaic
        of all of the said files.
        """
        if len(self.detectors) > 1:
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
    
    def addCatalogue(self, catalogue, obs_num):
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
        self._log("info","Adding catalogue %s" % (catalogue))
        cat = self.convertCatalogue(catalogue, obs_num)
        self._log("info","Finished converting catalogue to internal format")
        cats = [cat]
        for detector in self.detectors:
            self._log("info","Adding catalogue to detector %s" % (detector.name))
            cats.append(detector.addCatalogue(cat, dist=self.distortion))
        return cats
        self._log("info","Finished Adding Catalogue")
    
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

    def convertCatalogue(self,catalogue,obs_num):
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
        (catpath,catname) = os.path.split(catalogue)
        (catbase,ext) = os.path.splitext(catname)
        obsname = os.path.join(self.out_path,catbase+"_{:2d}_conv_{}.txt".format(obs_num, self.filter))
        t = read_metadata(catalogue, n_lines=1000)
        #Check for built-in metadata
        if 'keywords' in t.meta:
            if 'type' in t.meta['keywords']:
                table_type = t.meta['keywords']['type']['value']
        if table_type == 'phoenix':
            cat = self.convertPhoenixCatalogue(catalogue,obsname)
        elif table_type == 'phoenix_realtime':
            cat = self.convertPhoenixRealtime(catalogue, obsname)
        elif table_type == 'bc95':
            cat = self.convertBc95Catalogue(catalogue,obsname)
        elif table_type == 'internal':
            filter = t.meta['keywords']['filter']['value']
            if filter != self.filter: raise ValueError("Adding catalogue with filter {} to {} {}".format(filter, self.DETECTOR, self.filter))
            cat = catalogue
        elif table_type == 'mixed':
            filter = t.meta['keywords']['filter']['value']
            if filter != self.filter: raise ValueError("Adding catalogue with filter %s to NIRCamShort %s" % (filter,self.filter))
            cat = self.convertMixedCatalogue(catalogue,obsname)
        elif table_type == 'multifilter':
            if self.filter in t.columns:
                cat = self.convertMultiCatalogue(catalogue, obsname)
        else: #check for necessary columns
            #We want RA, DEC, and count rate in the appropriate filter
            if 'ra' in t.columns and 'dec' in t.columns and self.filter.lower() in t.columns:
                cat = self.convertGenericCatalogue(catalogue,obsname)
            else:
                raise ValueError("Can't parse catalogue without proper columns")
        del t
        return cat
    
    def convertPhoenixCatalogue(self, catalogue, obsname):
        """
        Converts a Phoenix grid of sources into the internal source table standard
        """
        self._log("info","Converting Phoenix catalogues")
        self._log("info", "Preparing output table")
        if os.path.isfile(obsname): 
            os.remove(obsname)
        
        with open(obsname, 'w') as f:
            f.write('\\ Internal Format Catalogue\n')
            f.write('\\ \n')
            f.write('\\ Parameters:\n')
            f.write('\\ \n')
            f.write('\\type="internal"\n')
            f.write('\\filter="%s"\n'%(self.filter))
            f.write('\\ \n')

        bp = self.bandpass
        cached = -1.
        for i, table in enumerate(read_table(catalogue)):
            ids = table['id']
            datasets = table['dataset']
            ras = table['ra']
            decs = table['dec']
            ages = table['age']
            metallicities = table['metallicity']
            masses = table['mass']
            distances = table['distance']
            binaries = table['binary']
            rates = numpy.zeros_like(ras)
            all_datasets = numpy.unique(datasets)
            self._log("info","Chunk {}: {} datasets".format(i+1, len(all_datasets)))
            for dataset in all_datasets:
                idx = numpy.where(datasets == dataset)
                if len(idx[0]) > 0:
                    stargen = StarGenerator(ages[idx][0], metallicities[idx][0], logger=self.logger)
                    rates[idx] = stargen.make_cluster_rates(masses[idx], self.INSTRUMENT, self.filter, bp, self.REFS)
                    del stargen
            rates = rates * 100 / (distances**2) #convert absolute to apparent rates
            if cached > 0:
                rates[0] += cached
                cached = -1.
            #Now, deal with binaries. Remember that if Binary == 1, then the star below is the binary companion.
            idx = numpy.where(binaries==1.0)[0] #Stars with binary companions
            idxp = (idx+1) #binary companions
            if len(idxp) > 0 and idxp[-1] >= len(rates): #last one is a binary. Cache it.
                cached = rates[-1]
                idx, idxp = idx[:-1], idxp[:-1]
                ids, datasets, ras, decs, ages, metallicities, = ids[:-1], datasets[:-1], ras[:-1], decs[:-1], ages[:-1], metallicities[:-1]
                masses, distances, binaries, rates = masses[:-1], distances[:-1], binaries[:-1], rates[:-1]
            rates[idx] += rates[idxp] #add count rates together
            # Now that we've added rates together, remove the binary companions
            ras = numpy.delete(ras,idxp)
            decs = numpy.delete(decs,idxp)
            rates = numpy.delete(rates,idxp)
            ids = numpy.delete(ids,idxp)
            binaries = numpy.delete(binaries,idxp)
            notes = numpy.empty_like(ras,dtype="S6")
            notes[numpy.where(binaries==1)] = 'Binary'
            notes[numpy.where(binaries!=1)] = 'None'
            types = numpy.full_like(ras,"point",dtype="S6")
        
            t = Table()
            t['ra'] = Column(data=ras)
            t['dec'] = Column(data=decs)
            t['flux'] = Column(data=rates)
            t['type'] = Column(data=types)
            t['n'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['re'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['phi'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['ratio'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['id'] = Column(data=ids)
            t['notes'] = Column(data=notes)

            with open(obsname, 'a') as outf:
                data = StringIO()
                t.write(data, format='ascii.ipac')
                data.seek(0)
                if i != 0: # get rid of the header lines
                    data.readline()
                    data.readline()
                    data.readline()
                    data.readline()
                outf.write(data.read())
                del data
            del t
        self._log("info","Finished output table")
        return obsname

    def convertPhoenixRealtime(self, catalogue, obsname):
        """
        Converts a set of Phoenix sources specified as (id, ra, dec, T, Z, log(g), apparent) into individual stars, observes them,
        and then produces an output catalogue
        """
        self._log("info","Converting Phoenix Realtime catalogues")
        self._log("info", "Preparing output table")
        if os.path.isfile(obsname): 
            os.remove(obsname)
        
        with open(obsname, 'w') as f:
            f.write('\\ Internal Format Catalogue\n')
            f.write('\\ \n')
            f.write('\\ Parameters:\n')
            f.write('\\ \n')
            f.write('\\type="internal"\n')
            f.write('\\filter="%s"\n'%(self.filter))
            f.write('\\ \n')

        ps.setref(**self.REFS)
        renorm_bp = ps.ObsBandpass('johnson,i')
        bp = self.bandpass
        cached = -1.
        for i, table in enumerate(read_table(catalogue)):
            ids = table['id']
            ras = table['ra']
            decs = table['dec']
            temps = table['teff']
            gravs = table['log_g']
            metallicities = table['metallicity']
            apparents = table['apparent']
            rates = numpy.zeros_like(ras)
            for index in range(len(ids)):
                t, g, Z, a = temps[index], gravs[index], metallicities[index], apparents[index]
                sp = ps.Icat('phoenix', t, Z, g)
                sp = sp.renorm(a, 'abmag', renorm_bp)
                obs = ps.Observation(sp, self.bandpass)
                rates[index] = obs.countrate()
            notes = numpy.full_like(ras,"None",dtype="S6")
            types = numpy.full_like(ras,"point",dtype="S6")
        
            t = Table()
            t['ra'] = Column(data=ras)
            t['dec'] = Column(data=decs)
            t['flux'] = Column(data=rates)
            t['type'] = Column(data=types)
            t['n'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['re'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['phi'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['ratio'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['id'] = Column(data=ids)
            t['notes'] = Column(data=notes)

            with open(obsname, 'a') as outf:
                data = StringIO()
                t.write(data, format='ascii.ipac')
                data.seek(0)
                if i != 0: # get rid of the header lines
                    data.readline()
                    data.readline()
                    data.readline()
                    data.readline()
                outf.write(data.read())
                del data
            del t
        self._log("info","Finished output table")
        return obsname

    def convertBc95Catalogue(self, catalogue, obsname):
        """
        Converts a BC95 galaxy grid of sources into the internal source table standard
        """
        self._log("info", "Converting BC95 Catalogue")
        proflist = {"expdisk":1,"devauc":4}
        if os.path.isfile(obsname): 
            os.remove(obsname) # No append
        
        with open(obsname, 'w') as f:
            f.write('\\ Internal Format Catalogue\n')
            f.write('\\ \n')
            f.write('\\ Parameters:\n')
            f.write('\\ \n')
            f.write('\\type="internal"\n')
            f.write('\\filter="{}"\n'.format(self.filter))
            f.write('\\ \n')

        ps.setref(**self.REFS)
        renorm_bp = ps.ObsBandpass('johnson,v')
        distance_type = "redshift"
        for i, table in enumerate(read_table(catalogue)):
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
            pas = (table['pa'] + (self.pa*180./numpy.pi) )%360.
            vmags = table['apparent']
            rates = numpy.array(())
            indices = numpy.array(())
            notes = numpy.array((),dtype='object')
            for (z,model,age,profile,radius,ratio,pa,mag) in zip(zs,models,ages,profiles,radii,ratios,pas,vmags):
                fname = "bc95_%s_%s.fits" % (model,age)
                try:
                    sp = ps.FileSpectrum(os.path.join(os.environ['PYSYN_CDBS'],"grid","bc95","templates",fname))
                    if distance_type == "redshift":
                        sp = sp.redshift(z)
                    sp = sp.renorm(mag,'abmag',renorm_bp)
                    obs = ps.Observation(sp,self.bandpass,force='taper')
                    rate = obs.countrate()
                except ValueError, v:
                    rate = 0.
                rates = numpy.append(rates,rate)
                indices = numpy.append(indices,proflist[profile])
                notes = numpy.append(notes,"BC95 %s %s %f" % (model,age,mag))
            types = numpy.full_like(ras,'sersic',dtype='S6')
        
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

            with open(obsname, 'a') as outf:
                data = StringIO()
                t.write(data, format='ascii.ipac')
                data.seek(0)
                if i != 0: # get rid of the header lines
                    data.readline()
                    data.readline()
                    data.readline()
                    data.readline()
                outf.write(data.read())
                del data
            del t

        return obsname

    def convertMixedCatalogue(self, catalogue, obsname):
        """
        Converts a mixed internal list of sources into the internal source table standard
        """

        if os.path.isfile(obsname): 
            os.remove(obsname) # No append
        
        with open(obsname, 'w') as f:
            f.write('\\ Internal Format Catalogue\n')
            f.write('\\ \n')
            f.write('\\ Parameters:\n')
            f.write('\\ \n')
            f.write('\\type="internal"\n')
            f.write('\\filter="{}"\n'.format(self.filter))
            f.write('\\ \n')

        for i, table in enumerate(read_table(catalogue)):
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
            idxp = numpy.where(units == 'p')
            rates[idxp] *= convertToCounts('p')
            idxe = numpy.where(units == 'e')
            rates[idxe] *= convertToCounts('e')
            idxj = numpy.where(units == 'j')
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

            with open(obsname, 'a') as outf:
                data = StringIO()
                t.write(data, format='ascii.ipac')
                data.seek(0)
                if i != 0: # get rid of the header lines
                    data.readline()
                    data.readline()
                    data.readline()
                    data.readline()
                outf.write(data.read())
                del data
            del t

        return obsname

    def convertMultiCatalogue(self, catalogue, obsname):
        """
        Converts an internal multifilter list of sources into the internal source table standard
        """

        if os.path.isfile(obsname): 
            os.remove(obsname) # No append
        
        with open(obsname, 'w') as f:
            f.write('\\ Internal Format Catalogue\n')
            f.write('\\ \n')
            f.write('\\ Parameters:\n')
            f.write('\\ \n')
            f.write('\\type="internal"\n')
            f.write('\\filter="%s"\n'%(self.filter))
            f.write('\\ \n')

        for i, table in enumerate(read_table(catalogue)):
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

            with open(obsname, 'a') as outf:
                data = StringIO()
                t.write(data, format='ascii.ipac')
                data.seek(0)
                if i != 0: # get rid of the header lines
                    data.readline()
                    data.readline()
                    data.readline()
                    data.readline()
                outf.write(data.read())
                del data
            del t

        return obsname

    def convertGenericCatalogue(self, catalogue, obsname):
        """
        Converts a generic list of point sources into the internal source table standard
        """

        if os.path.isfile(obsname): 
            os.remove(obsname) # No append

        with open(obsname, 'w') as f:
            f.write('\\ Internal Format Catalogue\n')
            f.write('\\ \n')
            f.write('\\ Parameters:\n')
            f.write('\\ \n')
            f.write('\\type="internal"\n')
            f.write('\\filter="%s"\n'%(self.filter))
            f.write('\\ \n')

        for i, table in enumerate(read_table(catalogue)):
            ras = table['ra']
            decs = table['dec']
            rates = table[self.filter.lower()]
            types = numpy.full_like(ras,"point",dtype="S6")
            if 'id' in table:
                ids = table['id']
            else:
                ids = numpy.arange(len(ras),dtype=int)
        
            t = Table()
            t['ra'] = Column(data=ras)
            t['dec'] = Column(data=decs)
            t['flux'] = Column(data=rates)
            t['type'] = Column(data=types)
            t['n'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['re'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['phi'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['ratio'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))
            t['id'] = Column(data=ids)
            t['notes'] = Column(data=numpy.full_like(ras,"N/A",dtype="S3"))

            with open(obsname, 'a') as outf:
                data = StringIO()
                t.write(data, format='ascii.ipac')
                data.seek(0)
                if i != 0: # get rid of the header lines
                    data.readline()
                    data.readline()
                    data.readline()
                    data.readline()
                outf.write(data.read())
                del data
            del t

        return obsname

    def generateReadnoise(self,exptime):
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
    
    def addError(self,poisson=True,readnoise=True,flat=True,dark=True,cosmic=True,exptime=1.,coadd=1):
        """Base function for adding in residual error"""
        self._log("info","Adding residual error")
        if flat:
            flat = AstroImage.initDataFromFits(self.flatfile,ext='COMPRESSED_IMAGE',logger=self.logger)
        if dark:
            dark = AstroImage.initDataFromFits(self.darkfile,ext='COMPRESSED_IMAGE',logger=self.logger)
            dark *= exptime
        if readnoise:
            rn = self.generateReadnoise(exptime)
        for detector in self.detectors:
            self._log("info","Adding error to detector %s" % (detector.name))
            self._log("info","Adding background")
            bkg = float(self.background/numpy.sqrt(coadd))
            self._log("info","Background is {}".format(bkg))
            detector += bkg
            self._log("info","Inserting correct exposure time")
            detector *= exptime
            self._log("info","Convolving with PSF")
            detector.convolve(self.psf)
            if self.oversample != 1:
                self._log("info","Binning oversampled image")
                detector.bin(self.oversample)
            if poisson: 
                self._log("info","Adding poisson noise")
                detector.introducePoissonNoise()
            if readnoise:
                self._log("info","Adding readnoise")
                detector.introduceReadnoise(rn)
            if flat:
                self._log("info","Adding flatfield residual")
                detector.introduceFlatfieldResidual(flat)
            if dark:
                self._log("info","Adding dark residual")
                detector.introduceDarkResidual(dark)
            if cosmic:
                self._log("info","Adding cosmic ray residual")
                detector.introduceCosmicRayResidual(self.PIXEL_SIZE,exptime)
        self._log("info","Finished adding error")

    @property
    def bandpass(self):
    
        translate = {
                        'wfi': 'wfirstimager'
                    }

        obsmode = {
                   'instrument': translate.get(self.INSTRUMENT.lower(), self.INSTRUMENT.lower()),
                   'mode': self.MODE,
                   'filter': self.filter.lower(),
                   'aperture': self.APERTURE,
                   'disperser': None
                   }

        conf = {'instrument': obsmode}
        
        self.logger.info("Creating Instrument with Configuration {}".format(obsmode))

        i = InstrumentFactory(config=conf)
        wr = i.get_wave_range()
        wave = numpy.linspace(wr['wmin'], wr['wmax'], num=500)
        pce = i.get_total_eff(wave)
        
        return ps.ArrayBandpass(wave=wave*1.e4, throughput=pce, waveunits='angstroms', name='bp_{}_{}'.format(self.instrument, self.filter))
    
    @property
    def zeropoint(self):
        ps.setref(**self.REFS)
        sp = ps.FileSpectrum(os.path.join(os.environ["PYSYN_CDBS"],"standards","alpha_lyr_stis_007.fits"))
        bp = self.bandpass
        sp = sp.renorm(0.0,"VEGAMAG",bp)
        obs = ps.Observation(sp,bp)
        zeropoint = obs.effstim("OBMAG")
        return zeropoint
    
    @property
    def photflam(self):
        ps.setref(**self.REFS)
        sp = ps.FlatSpectrum(0, fluxunits='stmag')
        bp = self.bandpass
        obs = ps.Observation(sp, bp)
        return obs.effstim('flam') / obs.countrate()

    def _log(self,mtype,message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self,'logger'):
            getattr(self.logger,mtype)(message)
        else:
            sys.stderr.write("%s: %s\n" % (mtype,message))


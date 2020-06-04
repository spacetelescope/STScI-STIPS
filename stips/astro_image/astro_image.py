from __future__ import absolute_import,division
__filetype__ = "base"

#External Modules
import logging, os, shutil, sys, time, uuid

import numpy as np

from astropy import wcs
from astropy.io import fits
from astropy.table import Table, Column
from copy import deepcopy
from photutils import CircularAperture, aperture_photometry
from photutils.psf.models import GriddedPSFModel

from scipy.ndimage.interpolation import zoom, rotate

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from cStringIO import StringIO

#Local Modules
from .. import __version__ as stips_version
from ..utilities import OffsetPosition
from ..utilities import overlapadd2
from ..utilities import overlapaddparallel
from ..utilities import read_table
from ..utilities import ImageData
from ..utilities import Percenter
from ..utilities import StipsDataTable
from ..errors import GetCrProbs, GetCrTemplate, MakeCosmicRay


class AstroImage(object):
    """
    The AstroImage class represents a generic astronomical image. The image has the following
    data associated with it:
        _file   : string of file name (including path) containing mem-mapped numpy array.
        data    : mem-mapped numpy double-precision 2D array of image data, in counts
        scale   : array of 2 double-precision floating point values, forming X and Y scale in
                  arcseconds/pixel
        wcs     : astropy WCS object containing image WCS information.
        header  : key/value array. Contains FITS header information and metadata
        history : array of strings holding the FITS HISTORY section
    """

    def __init__(self, **kwargs):
        """
        Astronomical image. The __init__ function creates an empty image with all other data values
        set to zero.
        """
        default = self.INSTRUMENT_DEFAULT

        if 'parent' in kwargs:
            self.parent = kwargs['parent']
            self.logger = self.parent.logger
            self.out_path = self.parent.out_path
            self.prefix = self.parent.prefix
            self.seed = self.parent.seed
            self.telescope = self.parent.TELESCOPE.lower()
            self.instrument = self.parent.PSF_INSTRUMENT
            self.filter = self.parent.filter
            self.oversample = self.parent.oversample
            self.shape = np.array(self.parent.DETECTOR_SIZE)*self.oversample
            self._scale = np.array(self.parent.SCALE)/self.oversample
            self.zeropoint = self.parent.zeropoint
            self.photflam = self.parent.photflam
            self.photplam = self.parent.PHOTPLAM[self.filter]
            background = self.parent.background
            self.psf_grid_size = self.parent.psf_grid_size
            self.psf_commands = self.parent.psf_commands
            small_subarray = self.parent.small_subarray
            self.cat_type = self.parent.cat_type
            self.set_celery = self.parent.set_celery
            self.get_celery = self.parent.get_celery
            self.convolve_size = self.parent.convolve_size
            self.memmap = self.parent.memmap
        else:
            self.parent = None
            if 'logger' in kwargs:
                self.logger = kwargs['logger']
            else:
                self.logger = logging.getLogger('__stips__')
                self.logger.setLevel(logging.INFO)
                if not len(self.logger.handlers):
                    stream_handler = logging.StreamHandler(sys.stderr)
                    format = '%(asctime)s %(levelname)s: %(message)s'
                    stream_handler.setFormatter(logging.Formatter(format))
                    self.logger.addHandler(stream_handler)
            self.out_path = kwargs.get('out_path', os.getcwd())
            self.oversample = kwargs.get('oversample', default['oversample'])
            self.shape = kwargs.get('shape', default['shape'])
            self.shape = np.array(self.shape) * self.oversample
            self._scale = kwargs.get('scale', np.array(default['scale']))
            self.prefix = kwargs.get('prefix', '')
            self.cat_type = kwargs.get('cat_type', 'fits')
            self.set_celery = kwargs.get('set_celery', None)
            self.get_celery = kwargs.get('get_celery', None)
            self.seed = kwargs.get('seed', 1234)
            small_subarray = kwargs.get('small_subarray', False)
            self.zeropoint = kwargs.get('zeropoint', default['zeropoint'])
            self.photflam = kwargs.get('photflam', default['photflam'])
            self.photplam = kwargs.get('photplam', default['photplam'])
            background = kwargs.get('background', default['background'])
            self.telescope = kwargs.get('telescope', default['telescope'])
            self.instrument = kwargs.get('instrument', default['instrument'])
            self.filter = kwargs.get('filter', default['filter'])
            self.psf_grid_size = kwargs.get('psf_grid_size', 
                                            default['psf_grid_size'])
            self.psf_commands = kwargs.get('psf_commands', '')
            self.convolve_size = kwargs.get('convolve_size', 8192)
            self.memmap = kwargs.get('memmap', True)
        
        if self.get_celery is None:
            self.get_celery = lambda: ""
        if self.set_celery is None:
            self.set_celery = lambda x: None

        #Set unique ID and figure out where the numpy memmap will be stored
        self.name = kwargs.get('detname', default['detector'][self.instrument])
        self.detector = self.name
        if self.memmap:
            fname = self.prefix+"_"+uuid.uuid4().hex+"_"+self.name+".tmp"
            self.fname = os.path.join(self.out_path, fname)

        if self.psf_commands is None:
            self.psf_commands = ''
        psf = kwargs.get('psf', True)
        if psf:
            self.make_psf()
        
        data = kwargs.get('data', None)
        if data is not None:
            base_shape = np.array(data.shape)
        else:
            #restrict data size to PSF size
            if small_subarray:
                if not hasattr(self, 'psf'):
                    msg = "{}: Unable to set image size to PSF size when image "
                    msg += "has no valid PSF."
                    raise ValueError(msg.format(self.name))
                base_shape = self.psf_shape
            else:
                base_shape = self.shape
        self._init_dat(base_shape, self.psf_shape, data)
        
        #Get WCS values if present, or set up a default
        self.wcs = self._getWcs(**kwargs)
        self._prepRaDec()
        
        #Header
        self.header = kwargs.get('header', {})
        if self.header == {}:
            self.header = kwargs.get('imh', {})
        self._prepHeader()
        if 'exptime' in self.header:
            self.exptime = self.header['exptime']
        else:
            self.exptime = kwargs.get('exptime', 1.)
        self.updateHeader('exptime', self.exptime)
        
        #History
        self.history = kwargs.get('history', [])
        
        #Special values for Sersic profile generation
        self.profile_multiplier = kwargs.get('profile_multiplier', 100.)
        self.noise_floor = max(background, kwargs.get('noise_floor', 1.))        

    
    def __del__(self):
        if os.path.exists(self.fname):
            os.remove(self.fname)

    def copy(self):
        other = AstroImage(out_path=self.out_path, detname=self.name, wcs=self.wcs, header=self.header, history=self.history,
                           xsize=self.xsize, ysize=self.ysize, zeropoint=self.zeropoint, photflam=self.photflam, 
                           logger=self.logger)
        try:
            if os.path.exists(other.fname):
                os.remove(other.fname)
            shutil.copy(self.fname, other.fname)
        except Exception as e:
            if os.path.exists(other.fname):
                remove(other.fname)
            raise e
        return other
    
    @classmethod
    def initFromFits(cls, file, **kwargs):
        """
            Takes the entire image from the FITS file (i.e. the FITS header overrides everything else).
            Looks for optional ext keyword to indicate which FITS extension to use.
        """
        img = cls(**kwargs)
        if file != '':
            try:
                with fits.open(file) as inf:
                    ext = kwargs.get('ext', 0)
                    dat = inf[ext].data
                    img._init_dat(base_shape=np.array(dat.shape), data=dat)
                    my_wcs = wcs.WCS(inf[ext].header)
                    for k,v in inf[ext].header.items():
                        if k != '' and k not in my_wcs.wcs.to_header():
                            img.header[k] = v
                img.wcs = img._normalizeWCS(my_wcs)
                img._prepRaDec()
                img._prepHeader()
                img.updateHeader("ASTROIMAGEVALID", True)
                img.addHistory("Created from FITS file %s" % (os.path.split(file)[1]))
                img._log("info","Created AstroImage %s from FITS file %s" % (img.name,os.path.split(file)[1]))
            except IOError as e:
                img.updateHeader("ASTROIMAGEVALID", False)
                img.addHistory("Attempted to create from invalid FITS file %s" % (os.path.split(file)[1]))
                img._log("warning","Attempted to create AstroImage %s from invalid FITS file %s" % (img.name,os.path.split(file)[1]))
                img._log("warning","Error is {}".format(e))
        return img
    
    @classmethod
    def initDataFromFits(cls, file, **kwargs):
        """Takes *only* data from the FITS, and everything else from kwargs (or sets to default)"""
        img = cls(**kwargs)
        if file != '':
            try:
                with fits.open(file) as inf:
                    ext = kwargs.get('ext', 0)
                    dat = inf[ext].data
                    img._init_dat(base_shape=np.array(dat.shape), data=dat)
                img.wcs = img._getWcs(**kwargs)
                img._prepRaDec()
                img._prepHeader()
                img.updateHeader("ASTROIMAGEVALID", True)
                img.addHistory("Data imported from FITS file %s" % (os.path.split(file)[1]))
                img._log("info","Created AstroImage %s and imported data from FITS file %s" % (img.name,os.path.split(file)[1]))
            except IOError as e:
                img.updateHeader("ASTROIMAGEVALID", False)
                img.addHistory("Attempted to create from invalid FITS file %s" % (os.path.split(file)[1]))
                img._log("warning","Attempted to create AstroImage %s from invalid FITS file %s" % (img.name,os.path.split(file)[1]))
        return img
    
    @classmethod
    def initFromArray(cls, array, **kwargs):
        """Convenience to initialize an image from a numpy array."""
        img = cls(data=array, **kwargs)
        img._log("info", "Creating AstroImage {} from array".format(img.name))
        return img
    
    @classmethod
    def initFromPoints(cls, xs, ys, rates, **kwargs):
        """Convenience to initialize a blank image and add a set of points to it"""
        img = cls(**kwargs)
        img._log("info","Creating AstroImage %s from points" % (img.name))
        if img.xsize < np.ceil(np.max(xs)) + 1:
            img.xsize = np.ceil(np.max(xs)) + 1
        if img.ysize < np.ceil(np.max(ys)) + 1:
            img.ysize = np.ceil(np.max(ys)) + 1
        img.addPoints(xs, ys, rates)
        return img
    
    @classmethod
    def initFromProfile(cls,posX,posY,flux,n,re,phi,axialRatio,**kwargs):
        """Convenience to initialize a blank image and add a sersic profile to it"""
        img = cls(**kwargs)
        img._log("info","Creating AstroImage %s from Sersic Profile" % (img.name))
        img.addSersicProfile(posX, posY, flux, n, re, phi, axialRatio)
        return img
    
    @property
    def xsize(self):
        return self.shape[1]
        
    @xsize.setter
    def xsize(self, size, offset=0):
        """Change the horizontal size. The offset will be applied to the new image before adding"""
        self.crop(size, self.ysize, offset, 0)
        self.shape = np.array((self.ysize, size))
    
    @property
    def ysize(self):
        return self.shape[0]
    
    @ysize.setter
    def ysize(self, size, offset=0):
        """Change the vertical size. The offset will be applied to the new image before adding"""
        self.crop(self.xsize, size, 0, offset)
        self.shape = np.array((size, self.xsize))
    
    @property
    def xscale(self):
        return abs(self.scale[0])*3600.
    
    @property
    def yscale(self):
        return abs(self.scale[1])*3600.
    
    @property
    def scale(self):
        return self._scale
    
    @property
    def rascale(self):
        return abs(self.scale[self.ranum])*3600.

    @property
    def decscale(self):
        return abs(self.scale[self.decnum])
    
    @property
    def distorted(self):
        return self.wcs.sip is not None
    
    @property
    def ra(self):
        if self.wcs.wcs.lngtyp == 'RA':
            return self.wcs.wcs.crval[self.wcs.wcs.lng]
        elif self.wcs.wcs.lattyp == 'RA':
            return self.wcs.wcs.crval[self.wcs.wcs.lat]
        else:
            raise ValueError("WCS has longitude %s and latitude %s. Can't get RA" % (self.wcs.wcs.lngtyp,self.wcs.wcs.lattyp))
    
    @ra.setter
    def ra(self,ra):
        if self.wcs.wcs.lngtyp == 'RA':
            self.wcs.wcs.crval[self.wcs.wcs.lng] = ra%360.
            self.addHistory("Set RA to %f" % (ra))
        elif self.wcs.wcs.lattyp == 'RA':
            self.wcs.wcs.crval[self.wcs.wcs.lat] = ra%360.
            self.addHistory("Set RA to %f" % (ra))
        else:
            raise ValueError("WCS has longitude %s and latitude %s. Can't set RA" % (self.wcs.wcs.lngtyp,self.wcs.wcs.lattyp))

    @property
    def dec(self):
        if self.wcs.wcs.lngtyp == 'DEC':
            return self.wcs.wcs.crval[self.wcs.wcs.lng]
        elif self.wcs.wcs.lattyp == 'DEC':
            return self.wcs.wcs.crval[self.wcs.wcs.lat]
        else:
            raise ValueError("WCS has longitude %s and latitude %s. Can't get DEC" % (self.wcs.wcs.lngtyp,self.wcs.wcs.lattyp))
    
    @dec.setter
    def dec(self,dec):
        if self.wcs.wcs.lngtyp == 'DEC':
            self.wcs.wcs.crval[self.wcs.wcs.lng] = dec
            self.addHistory("Set DEC to %f" % (dec))
        elif self.wcs.wcs.lattyp == 'DEC':
            self.wcs.wcs.crval[self.wcs.wcs.lat] = dec
            self.addHistory("Set DEC to %f" % (dec))
        else:
            raise ValueError("WCS has longitude %s and latitude %s. Can't set DEC" % (self.wcs.wcs.lngtyp,self.wcs.wcs.lattyp))
 
    @property
    def pa(self):
        """WCS is normalized, so PA from angle will work"""
        return self._getPA(self.wcs,self.scale,self.decnum)
    
    @pa.setter
    def pa(self,pa):
        cpa = np.cos(np.radians(pa%360.))
        spa = np.sin(np.radians(pa%360.))
        self.wcs.wcs.pc = np.array([[cpa, -spa], [spa, cpa]])
        self.addHistory("Set PA to %f" % (pa))

    
    @property
    def celery_state(self):
        if self.get_celery is not None:
            return self.get_celery()
        return ""
    
    @celery_state.setter
    def celery_state(self, state):
        if self.set_celery is not None:
            self.set_celery(state)

    
    @property
    def hdu(self):
        """Output AstroImage as a FITS Primary HDU"""
        with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
            hdu = fits.PrimaryHDU(dat, header=self.wcs.to_header(relax=True))
        hdu.header['CDELT1'] = self.scale[0]/3600.
        hdu.header['CDELT2'] = self.scale[0]/3600.
        if sys.version_info[0] >= 3:
            for k,v in self.header.items():
                if k != "ASTROIMAGEVALID":
                    hdu.header[k] = v
        else:
            for k,v in self.header.iteritems():
                if k != "ASTROIMAGEVALID":
                    hdu.header[k] = v
        for item in self.history:
            hdu.header.add_history(item)
        self._log("info","Created Primary HDU from AstroImage %s" % (self.name))
        return hdu
    
    @property
    def imageHdu(self):
        """Output AstroImage as a FITS Extension HDU"""
        self._log("info","Creating Extension HDU from AstroImage %s" % (self.name))
        with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
            hdu = fits.ImageHDU(dat, header=self.wcs.to_header(relax=True), name=self.name)
        if sys.version_info[0] >= 3:
            for k,v in self.header.items():
                hdu.header[k] = v
        else:
            for k,v in self.header.iteritems():
                hdu.header[k] = v
        for item in self.history:
            hdu.header.add_history(item)
        self._log("info","Created Extension HDU from AstroImage %s" % (self.name))
        return hdu

    @property
    def psf_constructor(self):
        import webbpsf
        return getattr(getattr(webbpsf, self.telescope), self.instrument)()
    
    
    @property
    def psf_shape(self):
        if hasattr(self, 'psf'):
            sampled_shape = self.psf.data.shape
            return np.array([sampled_shape[1], sampled_shape[2]])
        return (0, 0)

    
    def toFits(self, outFile):
        """Create a FITS file from the current state of the AstroImage data."""
        self._log("info","Writing AstroImage %s to FITS" % (self.name))
        hdulist = fits.HDUList([self.hdu])
        hdulist.writeto(outFile, overwrite=True)
    
    def updateHeader(self,k,v):
        """
        Updates a single keyword in the header dictionary, replacing the current value if there is
        one, otherwise adding a new value
        """
        self.header[k] = v
    
    def addHistory(self,v):
        """Adds an entry to the header history list."""
        self.history.append(v)
    
    def addTable(self, t, dist=False, *args, **kwargs):
        """
        Add a catalogue table to the Image. The Table must have the following columns:
            RA: RA of source
            DEC: DEC of source
            FLUX: flux of source
            TYPE: type of source (point, sersic)
            N: sersic index
            Re: radius containing half of the light of the sersic profile
            Phi: angle of the major axis of the sersic profile
            Ratio: axial ratio of the Sersic profile
            ID: id of source in catalogue
            Notes: any notes of important
        The following will then be done:
            - the table will be shifted from RA,DEC to X,Y (and items outside the FOV will be omitted)
            - the table will be split into point sources and sersic profiles
            - the point sources will be added via addPoints
            - the Sersic Profiles will be iteratively added via addSersicProfile
        The converted table (with X,Y instead of ra,dec and non-visible points removed) will be returned.
        """
        ras = t['ra']
        decs = t['dec']
        self._log("info", "Determining pixel co-ordinates")
        if dist and self.distorted:
            xs, ys = self.wcs.all_world2pix(t['ra'], t['dec'],1, quiet=True, adaptive=True, detect_divergence=True)
        else:
            xs, ys = self.wcs.wcs_world2pix(t['ra'], t['dec'],1)
        to_keep = np.where((xs >= 0) & (xs <= self.xsize) & (ys >= 0) & (ys <= self.ysize))
        self._log("info", "Keeping {} items".format(len(xs[to_keep])))
        ot = None
        base_state = self.getState()
        min_scale = min(self.xscale, self.yscale)
        if len(xs[to_keep]) > 0:
            xs = xs[to_keep]
            ys = ys[to_keep]
            xfs, yfs = self._remap(xs, ys)
            fluxes = t['flux'][to_keep]
            fluxes_observed = np.empty_like(fluxes)
            types = t['type'][to_keep]
            ns = t['n'][to_keep]
            res = t['re'][to_keep]
            phis = t['phi'][to_keep]
            ratios = t['ratio'][to_keep]
            ids = t['id'][to_keep]
            old_notes = t['notes'][to_keep]
            notes = np.empty_like(xs, dtype="S150")
            notes[:] = old_notes[:]
            vegamags = -2.512 * np.log10(fluxes) - self.zeropoint
            stmags = -2.5 * np.log10(fluxes * self.photflam) - 21.10
            stars_idx = np.where(types == 'point')
            if len(xs[stars_idx]) > 0:
                self.updateState(base_state + "<br /><span class='indented'>Adding {} stars</span>".format(len(xs[stars_idx])))
                self._log("info", "Writing {} stars".format(len(xs[stars_idx])))
                self.addPoints(xs[stars_idx], ys[stars_idx], fluxes[stars_idx], *args, **kwargs)
                fluxes_observed[stars_idx] = fluxes[stars_idx]
            gals_idx = np.where(types == 'sersic')
            if len(xs[gals_idx]) > 0:
                self._log("info","Writing {} galaxies".format(len(xs[gals_idx])))
                gxs = xs[gals_idx]
                gys = ys[gals_idx]
                gfluxes = fluxes[gals_idx]
                gtypes = types[gals_idx]
                gns = ns[gals_idx]
                gres = res[gals_idx]
                gphis = phis[gals_idx]
                gratios = ratios[gals_idx]
                gids = ids[gals_idx]
                counter = 1
                total = len(gxs)
                gfluxes_observed = np.empty_like(gxs, dtype='float32')
                gnotes = np.empty_like(gxs, dtype='S150')
                self._log('info', 'Starting Sersic Profiles at {}'.format(time.ctime()))
                for (x, y, flux, n, re, phi, ratio, id) in zip(gxs, gys, gfluxes, gns, gres, gphis, gratios, gids):
                    item_index = np.where(ids==id)[0][0]
                    self._log("info", "Index is {}".format(item_index))
                    self.updateState(base_state + "<br /><span class='indented'>Adding galaxy {} of {}</span>".format(counter, len(xs[gals_idx])))
                    central_flux = self.addSersicProfile(x, y, flux, n, re, phi, ratio, *args, **kwargs)
                    fluxes_observed[item_index] = central_flux
                    notes[item_index] = "{}: surface brightness {:.3f} yielded flux {:.3f}".format(notes[item_index], flux, central_flux)
                    self._log("info", "Finished Galaxy {} of {}".format(counter, total))
                    counter += 1
                self._log('info', 'Finishing Sersic Profiles at {}'.format(time.ctime()))
            ot = Table()
            ot['x'] = Column(data=xfs, unit='pixel')
            ot['y'] = Column(data=yfs, unit='pixel')
            ot['type'] = Column(data=types)
            ot['vegamag'] = Column(data=vegamags)
            ot['stmag'] = Column(data=stmags)
            ot['countrate'] = Column(data=fluxes_observed, unit='counts/s')
            ot['id'] = Column(data=ids)
            ot['notes'] = Column(data=notes)
        return ot
    
    def addCatalogue(self, cat, dist=False, *args, **kwargs):
        """
        Add a catalogue to the Image. The Catalogue must have the following columns:
            RA: RA of source
            DEC: DEC of source
            FLUX: flux of source
            TYPE: type of source (point, sersic)
            N: sersic index
            Re: radius containing half of the light of the sersic profile
            Phi: angle of the major axis of the sersic profile
            Ratio: axial ratio of the Sersic profile
            ID: id of source in catalogue
            Notes: any notes of important
        The following will then be done:
            - the catalogue will be shifted from RA,DEC to X,Y (and items outside the FOV will be
                omitted)
            - the catalogue will be split into point sources and sersic profiles
            - the point sources will be added via addPoints
            - the Sersic Profiles will be iteratively added via addSersicProfile
        """
        (path, catname) = os.path.split(cat)
        (catbase, catext) = os.path.splitext(catname)
        self._log("info","Adding catalogue %s to AstroImage %s" % (catname, self.name))
        obs_file_name = "{}_observed_{}.{}".format(catbase, self.name, self.cat_type)
        obsname = os.path.join(self.out_path, obs_file_name)
        self.addHistory("Adding items from catalogue %s" % (catname))
        data = None
        in_data = StipsDataTable.dataTableFromFile(cat)
        out_data = StipsDataTable.dataTableFromFile(obsname)
        out_data.meta = {'name': 'Observed Catalogue', 'type': 'observed', 'detector': self.name, 'input_catalogue': catname}
        base_state = self.getState()
        counter = 0
        current_chunk = in_data.read_chunk()
        while current_chunk is not None:
            table_length = len(current_chunk['id'])
            self.updateState(base_state + "<br /><span class='indented'>Adding sources {} to {}</span>".format(counter, counter+table_length))
            out_chunk = self.addTable(current_chunk, dist, *args, **kwargs)
            if out_chunk is not None:
                out_data.write_chunk(out_chunk)
            counter += table_length
            current_chunk = in_data.read_chunk()
        self.updateState(base_state)
        self._log("info","Added catalogue %s to AstroImage %s" % (catname, self.name))
        return obsname            
    
    def addPoints(self, xs, ys, rates, *args, **kwargs):
        """Adds a set of point sources to the image given their co-ordinates and count rates."""
        self.addHistory("Adding %d point sources" % (len(xs)))
        self._log("info","Adding %d point sources to AstroImage %s" % (len(xs),self.name))
        xs = np.floor(xs).astype(int)
        ys = np.floor(ys).astype(int)
        with ImageData(self.fname, self.shape, memmap=self.memmap) as dat:
            dat[ys, xs] += rates
    
    def addSersicProfile(self, posX, posY, flux, n, re, phi, axialRatio, *args, **kwargs):
        """
        Adds a single sersic profile to the image given its co-ordinates, count rate, and source 
        type.
        
        (posX,posY) are the co-ordinates of the centre of the profile (pixels).
        flux is the total number of counts to add to the AstroImage.
        n is the Sersic profile index.
        re is the radius enclosing half of the total light (pixels).
        phi is the angle of the major axis (degrees east of north).
        axialRatio is the ratio of major axis to minor axis.
        """
        if flux == 0.:
            return 0.
        self.addHistory("Adding Sersic profile at (%f,%f) with flux %f, index %f, Re %f, Phi %f, and axial ratio %f" % (posX,posY,flux,n,re,phi,axialRatio))
        self._log("info","Adding Sersic: re={}, n={}, flux={}, phi={:.1f}, ratio={}".format(re, n, flux, phi, axialRatio))
        
        # Determine necessary parameters for the Sersic model -- the input radius, surface 
        # brightness and noise floor are all in *detector* pixels.
        pixel_radius = re * self.oversample
        pixel_brightness = flux / (self.oversample*self.oversample)
        noise_floor = self.noise_floor / (self.oversample*self.oversample)

        # Determine the pixel offset of the profile from the centre of the AstroImage
        # Centre of profile is (xc, yc)
        # Centre of image is (self.xsize//2, self.ysize//2)
        # Centre of profile on image is (posX, posY)
        # Let's say we have a 20X20 profile, centre is 10,10
        # Let's say our image is 1024X1024, centre is 512,512
        # Let's say posX=39, posY=27
        # So we want to go from 512,512 to 39,27, so offset is -473,-485
        # Offset is posX-image_centre, posY-image_centre
        offset_x, offset_y = np.floor(posX) - self.xsize//2, np.floor(posY) - self.ysize//2
        fractional_x, fractional_y = posX - np.floor(posX), posY - np.floor(posY)

        from astropy.modeling.models import Sersic2D
        # Figure out an appropriate radius. Start at 5X pixel radius, and continue until the highest value on the outer edge is below the noise floor.
        max_outer_value = 2*noise_floor
        filled = False
        radius_multiplier = 2.5
        full_frame = False
        model_size = int(np.ceil(pixel_radius*radius_multiplier))
        while max_outer_value > noise_floor:
            radius_multiplier *= 2
            model_size = int(np.ceil(pixel_radius*radius_multiplier))
            if not self._filled(offset_x, offset_y, model_size, model_size):
#                 self._log("info", "Creating a {}x{} array for the Sersic model at ({},{})".format(model_size, model_size, posX, posY))
                x, y, = np.meshgrid(np.arange(model_size), np.arange(model_size))
                # In order to get fractional flux per pixel correct, carry the non-integer portion of the model centre through.
                xc, yc = model_size//2 + fractional_x, model_size//2 + fractional_y
                mod = Sersic2D(amplitude=pixel_brightness, r_eff=pixel_radius, n=n, x_0=xc, y_0=yc, ellip=(1.-axialRatio), theta=(np.radians(phi) + 0.5*np.pi))
                img = mod(x, y)
                max_outer_value = max(np.max(img[0,:]), np.max(img[-1,:]), np.max(img[:,0]), np.max(img[:,-1]))
#                 self._log('info', "Max outer value is {}, noise floor is {}".format(max_outer_value, noise_floor))
            else:
                full_frame = True
#                 self._log("info", "Creating full-frame Sersic model at ({},{})".format(posX, posY))
                x, y = np.meshgrid(np.arange(self.ysize), np.arange(self.xsize))
                xc, yc = posX, posY
                mod = Sersic2D(amplitude=pixel_brightness, r_eff=pixel_radius, n=n, x_0=xc, y_0=yc, ellip=(1.-axialRatio), theta=(np.radians(phi) + 0.5*np.pi))
                img = mod(x, y)
                max_outer_value = 0.
        img = np.where(img >= noise_floor, img, 0.)
        aperture = CircularAperture((xc, yc), pixel_radius)
        flux_table = aperture_photometry(img, aperture)
        central_flux = flux_table['aperture_sum'][0]
        self._log("info", "Sersic profile has final size {}x{}, maximum value {}, sum {}".format(model_size, model_size, np.max(img), np.sum(img)))
        self._addArrayWithOffset(img, offset_x, offset_y)
        return central_flux

    def make_psf(self):
        from webbpsf import __version__ as psf_version
        have_psf = False
        psf_name = "psf_{}_{}_{}_{}_{}_{}.fits".format(self.instrument,
                                                       stips_version, 
                                                       self.filter,
                                                       self.oversample,
                                                       self.psf_grid_size,
                                                       self.detector)
        if os.path.exists(os.path.join(self.out_path, "psf_cache")):
            psf_file = os.path.join(self.out_path, "psf_cache", psf_name)
            if os.path.exists(psf_file):
                from webbpsf.utils import to_griddedpsfmodel
                if (self.psf_commands is None or self.psf_commands == ''):
                    self.psf = to_griddedpsfmodel(psf_file)
                    have_psf = True
        if not have_psf:
            base_state = self.celery_state
            update_state = "<br /><span class='indented'>Generating PSF</span>"
            self.celery_state = base_state + update_state
            ins = self.psf_constructor
            if self.psf_commands != '':
                for attribute,value in self.psf_commands.iteritems():
                    setattr(ins,attribute,value)
            ins.filter = self.filter
            ins.detector = self.detector
            scale = self.scale[0]
            # First limit -- PSF no larger than detector
            ins_size = max(self.xsize, self.ysize) * self.oversample
            # Second limit -- PSF no larger than half of max convolution area.
            conv_size = self.convolve_size // (2*self.oversample)
            # Third limit -- prevent aliasing
            safe_size = int(np.floor(30. * self.photplam / (2 * self.scale[0])))
            if safe_size <= 0:
                safe_size = max(self.xsize, self.ysize)
            msg = "PSF choosing between {}, {} and {}"
            self._log("info", msg.format(ins_size, conv_size, safe_size))
            fov_pix = min(ins_size, conv_size, safe_size)
            if fov_pix%2 != 0:
                fov_pix += 1
            num_psfs = self.psf_grid_size*self.psf_grid_size
            if os.path.exists(os.path.join(self.out_path, "psf_cache")):
                save = True
                overwrite = True
                psf_dir = os.path.join(self.out_path, "psf_cache")
                psf_file = "psf_{}_{}_{}_{}_{}".format(self.instrument,
                                                       stips_version, 
                                                       self.filter,
                                                       self.oversample,
                                                       self.psf_grid_size)
            else:
                save = False
                overwrite = False
                psf_dir = None
                psf_file = None
            
            msg = "{}: Starting {}x{} PSF Grid creation at {}"
            self._log("info", msg.format(self.name, self.psf_grid_size, 
                                         self.psf_grid_size, time.ctime()))
            self.psf = ins.psf_grid(all_detectors=False, num_psfs=num_psfs,
                                    fov_pixels=fov_pix, normalize='last',
                                    oversample=self.oversample, save=save,
                                    outdir=psf_dir, outfile=psf_file,
                                    overwrite=overwrite)
            msg = "{}: Finished PSF Grid creation at {}"
            self._log("info", msg.format(self.name, time.ctime()))
            self.celery_state = base_state
    
    def convolve_psf(self, max_size=None, parallel=False, cores=None):
        """
        Convolve the current AstroImage state with the generated PSF (if there
        is a generated PSF). Otherwise, do nothing.
        
        Parameters
        ----------
        max_size : int, default=4095
            The maximum size chunk to use in convolution chunks.
        parallel : bool, default=False
            Whether to perform convolution chunks in parallel
        cores : int, default=None
            How many CPU cores to use for parallel computation (used only if
            parallel=True)
        """
        if hasattr(self, 'psf'):
            if max_size is None:
                max_size = self.convolve_size
            self.convolve(self.psf, max_size=max_size, parallel=parallel,
                          cores=cores, crop=True)


    def convolve(self, other, max_size=None, parallel=False, cores=None,
                 crop=True):
        """
        Convolve the AstroImage with another image. This image can be
            - another AstroImage
            - a numpy NDData array
            - a FITS file
            - a webbpsf PSF grid
        
        Parameters
        ----------
        other : object
            The other image to convolve. See above for possible formats.
        max_size : int, default=4095
            The maximum size chunk to use in convolution chunks.
        parallel : bool, default=False
            Whether to perform convolution chunks in parallel
        cores : int, default=None
            How many CPU cores to use for parallel computation (used only if
            parallel=True)
        crop : bool, default=True
            After convolving, should the AstroImage be cropped down to no
            longer include the PSF overlap region.
        """
        if max_size is None:
            max_size = self.convolve_size
        
        other_type = ""
        if isinstance(other, AstroImage):
            other_type = "astro_image"
            other_data = ImageData(other.fname, other.shape, mode='r', 
                                   memmap=other.memmap)
            other_img = other_data.data
            other_shape = other_img.shape
            other_name = other
        elif isinstance(other, np.ndarray):
            other_type = "ndarray"
            other_img = other
            other_shape = other_img.shape
            other_name = "{}x{} array".format(other.shape[0], other.shape[1])
        elif isinstance(other, GriddedPSFModel):
            other_type = "psfgrid"
            other_img = other
            other_shape = self.psf_shape
            other_name = "PSF grid"
        else:
            # There are so darn many fits HDU classes that I can't
            #   figure out how to do this with them, so FITS file is now
            #   a diagnosis of last resort.
            other_type = 'fits_file'
            if 'primary' in other:
                other_img = other['primary'].data
            elif len(other) > 1:
                other_img = other[1].data
            elif hasattr(other, 'data'):
                other_img = other.data
            else:
                other_img = other[0].data
            other_shape = other_img.shape
            other_name = "FITS file"

        self.addHistory("Convolving {} with {}".format(self.name, other_name))
        self._log("info", "Convolving {} with {}".format(self.name, other_name))
        
        f = os.path.join(self.out_path, uuid.uuid4().hex+"_convolve_01.tmp")
        g = os.path.join(self.out_path, uuid.uuid4().hex+"_convolve_02.tmp")
        try:
            self_y, self_x = self.shape
            other_y, other_x = other_shape
            max_y = min(max_size - other_y, self_y + other_y - 1)
            max_x = min(max_size - other_x, self_x + other_x - 1)
            sub_shape = (max_y, max_x)
            with ImageData(self.fname, self.shape, mode='r', memmap=self.memmap) as dat:
                shape = (self_y + other_y - 1, self_x + other_x - 1)
                if self.memmap:
                    fp_res = np.memmap(f, dtype='float32', mode='w+', shape=shape)
                else:
                    fp_res = np.zeros(shape, dtype='float32')
                centre = (fp_res.shape[0]//2, fp_res.shape[1]//2)
                max_y = min(max_size - other_y, self_y + other_y - 1)
                max_x = min(max_size - other_x, self_x + other_x - 1)
                sub_shape = (max_y, max_x)
                msg = "PSF Shape: {}; Current Shape: {}"
                self._log('info', msg.format(other_shape, self.shape))
                msg = "Choosing between {}-{}={} and {}+{}-1={}"
                self._log('info', msg.format(max_size, other_y, 
                                             max_size-other_y, other_y, 
                                             self_y, other_y+self_y-1))
                msg = "Using overlapping arrays of size {}"
                self._log('info', msg.format(sub_shape))
                msg = "{}: Starting Convolution at {}"
                self._log('info', msg.format(self.name, time.ctime()))
                if parallel:
                    self._log('info', 'Convolving in parallel')
                    if self.memmap:
                        del fp_res
                    else:
                        f = fp_res
                    overlapaddparallel(dat, dat.shape, other_img, other_shape, 
                                       sub_shape, y=f, verbose=True, 
                                       logger=self.logger, 
                                       base_state=self.get_celery(), 
                                       state_setter=self.set_celery, 
                                       path=self.out_path, cores=cores,
                                       memmap=self.memmap)
                    if self.memmap:
                        fp_res = np.memmap(f, dtype='float32', mode='r+', 
                                           shape=shape)
                else:
                    overlapadd2(dat, dat.shape, other_img, other_shape, 
                                sub_shape, y=fp_res, verbose=True, 
                                logger=self.logger, 
                                base_state=self.get_celery(), 
                                state_setter=self.set_celery)
                msg = "{}: Finished Convolution at {}"
                self._log('info', msg.format(self.name, time.ctime()))

            if crop:
                msg = "Cropping convolved image down to detector size"
                self._log('info', msg)
                half = (self.base_shape[0]//2, self.base_shape[1]//2)
                msg = "Image Centre: {}; Image Half-size: {}"
                self._log('info', msg.format(centre, half))
                ly, hy = centre[0]-half[0], centre[0]+half[0]
                lx, hx = centre[1]-half[1], centre[1]+half[1]
                if hx-lx < self.base_shape[1]:
                    hx += 1
                elif hx-lx > self.base_shape[1]:
                    hx -= 1
                if hy-ly < self.base_shape[0]:
                    hy += 1
                elif hy-ly > self.base_shape[0]:
                    hy -= 1
                msg = "Taking [{}:{}, {}:{}]"
                self._log('info', msg.format(ly, hy, lx, hx))
                if self.memmap:
                    fp_crop = np.memmap(g, dtype='float32', mode='w+', 
                                        shape=tuple(self.base_shape))
                else:
                    fp_crop = np.zeros(tuple(self.base_shape), dtype='float32')
                fp_crop[:,:] = fp_res[ly:hy, lx:hx]
                crpix = [half[0], half[1]]
                if self.wcs.sip is not None:
                    sip = wcs.Sip(self.wcs.sip.a, self.wcs.sip.b, None, None, 
                                  crpix)
                else:
                    sip = None
                self.wcs = self._wcs(self.ra, self.dec, self.pa, self.scale, 
                                     crpix=crpix, sip=sip)
                del fp_res
                if self.memmap:
                    del fp_crop
                    if os.path.exists(self.fname):
                        os.remove(self.fname)
                    self.fname = g
                else:
                    del self.fname
                    self.fname = fp_crop
                self.shape = self.base_shape
            if os.path.exists(f):
                os.remove(f)
        except Exception as e:
            if os.path.exists(f):
                os.remove(f)
            if os.path.exists(g):
                os.remove(g)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e

    def rotate(self,angle,reshape=False):
        """
        Rotate the image a number of radians as specified
        
        ..warning:: This function is not necessarily flux-conserving
        """
        self.addHistory("Rotating by %f degrees" % (angle))
        self._log("info","Rotating AstroImage %s by %f degrees" % (self.name,angle))
        self.pa = (self.pa + angle)%360.%360.
        f = os.path.join(self.out_path, uuid.uuid4().hex+"_rotate.tmp")
        try:
            t_shape = tuple(self.shape)
            if self.memmap:
                fp_result = np.memmap(f, dtype='float32', mode='w+', shape=t_shape)
            else:
                fp_result = np.zeros(t_shape, dtype='float32')
            with ImageData(self.fname, self.shape, mode='r', memmap=self.memmap) as dat:
                rotate(dat, angle, order=5, reshape=reshape, output=fp_result)
            if self.memmap:
                del fp_result
                if os.path.exists(self.fname):
                    os.remove(self.fname)
                self.fname = f
            else:
                del self.fname
                self.fname = fp_result
        except Exception as e:
            if os.path.exists(f):
                os.remove(f)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e

    def addWithOffset(self,other,offset_x,offset_y):
        """
        Adds other to self, but with an offset.
        
        offset_n have units of pixels-of-self
        """
        self.addHistory("Adding image with offset of (%f,%f) pixels" % (offset_x,offset_y))
        self._log("info","Adding image %s with offset (%d,%d)" % (other.name,offset_x,offset_y))
        other_copy = other.copy()
        if abs(other.scale[0] - self.scale[0]) > 1.e-5 or abs(other.scale[1]-self.scale[1]) > 1.e-5:
            other_copy.rescale(self.scale)
        self._addWithOffset(other_copy, offset_x, offset_y)
    
    def _findCentre(self, xs, ys, offset_x, offset_y):
        """
        Determines the actual pixel offset for adding another image or array with shape (ys,xs) and
        offset (offset_x, offset_y) to the current image.
        """
        yc, xc = ys//2, xs//2
        
        ycenter, xcenter = self.ysize//2, self.xsize//2
        
        y_overlay = ycenter + offset_y
        x_overlay = xcenter + offset_x
        
        low_x = x_overlay - xc #0 on other array
        low_y = y_overlay - yc #0 on other array
        high_x = low_x + xs #xs on other array
        high_y = low_y + ys #ys on other array
        
        lx = 0
        ly = 0
        hx = xs
        hy = ys
        
        if low_x < 0:
            lx -= low_x
            low_x = 0
        if low_y < 0:
            ly -= low_y
            low_y = 0
        if high_x > self.xsize:
            hx = hx - high_x + self.xsize
            high_x = self.xsize
        if high_y > self.ysize:
            hy = hy - high_y + self.ysize
            high_y = self.ysize
        
        self._log("info","Adding Other[%d:%d,%d:%d] to AstroImage %s[%d:%d,%d:%d]" % (ly,hy,lx,hx,self.name,low_y,high_y,low_x,high_x))
        return int(lx), int(ly), int(hx), int(hy), int(low_x), int(low_y), int(high_x), int(high_y)

    def _filled(self, offset_x, offset_y, size_x, size_y):
        """
        If an image with size (size_x, size_y) is added to the current image with the pixel offset
            (offset_x, offset_y), will that image fill the frame of the current image?
        """
        xs, ys = size_x, size_y
        lx, ly, hx, hy, low_x, low_y, high_x, high_y = self._findCentre(xs, ys, offset_x, offset_y)
        
        if low_y==0 and low_x==0 and high_y==self.shape[0] and high_x==self.shape[1]:
            return True
        return False
        
    
    def _addArrayWithOffset(self, other, offset_x, offset_y):
        """
        Add an AstroImage to the current image with a given pixel offset. This allows for a smaller
        memory use for addWithAlignment.
        """
        ys, xs = other.shape
        lx, ly, hx, hy, low_x, low_y, high_x, high_y = self._findCentre(xs, ys, offset_x, offset_y)
        
        #If any of these are false, the images are disjoint.
        if low_x < self.xsize and low_y < self.ysize and high_x > 0 and high_y > 0:
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
                d = dat[low_y:high_y, low_x:high_x]
                d = other[ly:hy, lx:hx]
                dat[low_y:high_y, low_x:high_x] += other[ly:hy, lx:hx]
        else:
            self.addHistory("Added image is disjoint")
            self._log("warning","%s: Image is disjoint" % (self.name))
        
    def _addWithOffset(self, other, offset_x, offset_y):
        """
        Add an AstroImage to the current image with a given pixel offset. This allows for a smaller
        memory use for addWithAlignment.
        """
        ys, xs = other.shape
        lx, ly, hx, hy, low_x, low_y, high_x, high_y = self._findCentre(xs, ys, offset_x, offset_y)
        
        #If any of these are false, the images are disjoint.
        if low_x < self.xsize and low_y < self.ysize and high_x > 0 and high_y > 0:
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat, ImageData(other.fname, other.shape, mode='r', memmap=other.memmap) as other_data:
                dat[low_y:high_y, low_x:high_x] += other_data[ly:hy, lx:hx]
        else:
            self.addHistory("Added image is disjoint")
            self._log("warning","%s: Image is disjoint" % (self.name))
        
    def addWithAlignment(self,other):
        """Adds other to self after aligning co-ordinates"""
        self.addHistory("Adding image aligned by WCS")
        self._log("info","Adding image %s (RA,DEC,PA)=(%f,%f,%f) aligned by WCS" % (other.name,other.ra,other.dec,other.pa))
        other_copy = other.copy()
        if abs(other.pa - self.pa) > 1.e-5:
            self._log("info","Rotating other image by %f degrees" % ((self.pa-other.pa)))
            other_copy.rotate((self.pa-other.pa), reshape=True)
        if abs(other.scale[0] - self.scale[0]) > 1.e-5 or abs(other.scale[1]-self.scale[1]) > 1.e-5:
            self._log("info","Rescaling other from (%f,%f) to (%f,%f)" % (other.scale[0],other.scale[1],self.scale[0],self.scale[1]))
            other_copy.rescale(self.scale)
            self._log("info","Finished rescaling")
        ra = other.ra
        dec = other.dec
        pix = self.wcs.wcs_world2pix([ra],[dec],1,ra_dec_order=True)
        px = pix[0][0]
        py = pix[1][0]
        offset_x = int(np.round(px - self.wcs.wcs.crpix[0]))
        offset_y = int(np.round(py - self.wcs.wcs.crpix[1]))
        self._log("info","Other has centre co-ordinates (%f,%f) for offset (%d,%d)" % (px,py,offset_x,offset_y))
        self._addWithOffset(other_copy, offset_x, offset_y)

    def rescale(self,scale):
        """Rescales the image to the provided plate scale."""
        self.addHistory("Rescaling to (%f,%f) arcsec/pixel" % (scale[0],scale[1]))
        self._log("info","Rescaling to (%f,%f) arcsec/pixel" % (scale[0],scale[1]))
        f = os.path.join(self.out_path, uuid.uuid4().hex+"_scale.tmp")
        try:
            shape_x = int(round(self.shape[1] * self.scale[0] / scale[0]))
            shape_y = int(round(self.shape[0] * self.scale[1] / scale[1]))
            new_shape = (shape_y, shape_x)
            self._log("info","New shape will be {}".format(new_shape))
            if self.memmap:
                fp_result = np.memmap(f, dtype='float32', mode='w+', shape=new_shape)
            else:
                fp_result = np.zeros(new_shape, dtype='float32')
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
                flux = dat.sum()
                self._log("info", "Max flux is {}, sum is {}".format(np.max(dat), flux))
                zoom(dat, (np.array(self.scale)/np.array(scale)), fp_result)
            factor = flux / fp_result.sum()
            fp_result *= factor
            self._log("info", "Max flux is {}, sum is {}".format(np.max(fp_result), np.sum(fp_result)))
            if self.memmap:
                del fp_result
                if os.path.exists(self.fname):
                    os.remove(self.fname)
                self.fname = f
            else:
                del self.fname
                self.fname = fp_result
            self.shape = new_shape
        except Exception as e:
            if os.path.exists(f):
                os.remove(f)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e
        self._scale = scale
        self.wcs = self._wcs(self.ra, self.dec, self.pa, scale)
        self._prepHeader()
        
    def crop(self,xs,ys,offset_x,offset_y):
        """Crops the image. If the offset < 0, or offset+size > current size, pads with blank pixels"""
        self._addWithOffset(self.copy(), offset_x, offset_y)
        offset_value = np.array([offset_x,offset_y])
        pix_coords = np.array(offset_value+self.wcs.wcs.crpix)
        if self.distorted:
            world_coords = self.wcs.all_pix2world(pix_coords[0],pix_coords[1],1,ra_dec_order=True)
        else:
            world_coords = self.wcs.wcs_pix2world(pix_coords[0],pix_coords[1],1,ra_dec_order=True)
        self.wcs = self._wcs(world_coords[0],world_coords[1],self.pa,self.scale,self.wcs.sip)
        self._prepHeader()
    
    def bin(self,binx,biny=None):
        """Bin xXy pixels into a single pixel. Adjust the scale and size accordingly"""
        if biny is None: biny = binx
        f = os.path.join(self.out_path, uuid.uuid4().hex+"_bin.tmp")
        try:
            shape_x, shape_y = int(self.shape[1] // binx), int(self.shape[0] // biny)
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
                binned = dat.reshape(shape_y, biny, shape_x, binx).sum(axis=(1, 3))
                if self.memmap:
                    mapped = np.memmap(f, dtype='float32', mode='w+', shape=binned.shape)
                    mapped[:] = binned[:]
                    if os.path.exists(self.fname):
                        os.remove(self.fname)
                    self.fname = f
                else:
                    del self.fname
                    self.fname = binned
            self.shape = np.array((shape_y, shape_x))
        except Exception as e:
            if os.path.exists(f):
                os.remove(f)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e
        self.wcs = self._wcs(self.ra, self.dec, self.pa, [self.scale[0]*binx, self.scale[1]*biny])
        self._prepHeader()
        self.oversample = 1
    
    def setExptime(self, exptime):
        """
        Set the exposure time. Multiply data by new_exptime / old_exptime.
        """
        factor = exptime / self.exptime
        with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
            dat *= factor
        self.exptime = exptime
        self.updateHeader('exptime', self.exptime)
    
    def addBackground(self, background):
        """
        Add a constant value per-pixel. Automatically scale that value based on oversampling.
        """
        per_pixel_background = background / (self.oversample*self.oversample)
        self.addHistory("Added background of {} counts/s/detector pixel ({} counts/s/oversampled pixel)".format(background, per_pixel_background))
        self._log("info", "Added background of {} counts/s/detector pixel ({} counts/s/oversampled pixel)".format(background, per_pixel_background))
        with ImageData(self.fname, self.shape, memmap=self.memmap) as dat:
            dat += per_pixel_background

    def introducePoissonNoise(self,absVal=False):
        """
        Generate Poisson noise and add it to the internal image.

        Based on `CALC_SHOT_NOISE` in `instrument__init.pro`
        from JWST IDL Simulator package.

        Parameters
        ----------
        absVal: bool, optional
            Take absolute value of `noiseData`.
        """
        a, n = os.path.join(self.out_path, uuid.uuid4().hex+"_poisson_a.tmp"), os.path.join(self.out_path, uuid.uuid4().hex+"_poisson_n.tmp")
        try:
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
                t_shape = tuple(self.shape)
                if self.memmap:
                    abs_data = np.memmap(a, dtype='float32', mode='w+', shape=t_shape)
                else:
                    abs_data = np.zeros(t_shape, dtype='float32')
                np.absolute(dat, abs_data)

                if self.memmap:
                    noise_data = np.memmap(n, dtype='float32', mode='w+', shape=t_shape)
                else:
                    noise_data = np.zeros(t_shape, dtype='float32')
                noise_data[:,:] = np.random.RandomState(seed=self.seed).normal(size=self.shape) * np.sqrt(abs_data)
                del abs_data
                if absVal:
                    noise_data[:,:] = np.abs(noise_data)
                mean, std = noise_data.mean(), noise_data.std()
                dat += noise_data
                del noise_data
            if os.path.exists(a):
                os.remove(a)
            if os.path.exists(n):
                os.remove(n)
        except Exception as e:
            if os.path.exists(a):
                os.remove(a)
            if os.path.exists(n):
                os.remove(n)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e
        self.addHistory("Adding Poisson Noise with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding Poisson Noise with mean {} and standard deviation {}".format(mean, std))            

    def introduceReadnoise(self,readnoise):
        """
        Simulate kTC or read noise for detector.

        Parameters
        ----------
        readnoise: constant representing average read noise per pixel.
        """
        n = os.path.join(self.out_path, uuid.uuid4().hex+"_readnoise.tmp")
        try:
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
                t_shape = tuple(self.shape)
                if self.memmap:
                    noise_data = np.memmap(n, dtype='float32', mode='w+', shape=t_shape)
                else:
                    noise_data = np.zeros(t_shape, dtype='float32')
                noise_data[:,:] = readnoise * np.random.RandomState(seed=self.seed).randn(self.ysize,self.xsize)            
                mean, std = noise_data.mean(), noise_data.std()
                dat += noise_data
                del noise_data
            if os.path.exists(n):
                os.remove(n)
        except Exception as e:
            if os.path.exists(n):
                os.remove(n)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e
        self.addHistory("Adding Read noise with mean %f and standard deviation %f" % (mean, std))
        self._log("info","Adding readnoise with mean %f and STDEV %f" % (mean, std))

    def introduceFlatfieldResidual(self,flat):
        """
        Simulate flatfield correction error.
        
        flat: AstroImage containing flatfield error values.
        
        returns: mean, std. Mean and standard deviation of used portion of error array.
        """
        with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat, ImageData(flat.fname, flat.shape, mode='r', memmap=flat.memmap) as flat_data:
            err = flat_data[:self.ysize,:self.xsize]
            mean, std = err.mean(), err.std()
            dat *= err
        self.addHistory("Adding Flatfield residual with mean %f and standard deviation %f" % (mean, std))
        self._log("info","Adding Flatfield residual with mean %f and standard deviation %f" % (mean, std))

    def introduceDarkResidual(self,dark):
        """
        Simulate dark correction error.
        
        dark: AstroImage containing dark residual.
        
        returns: mean,std: mean and standard deviation of dark error array.
        """
        with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat, ImageData(dark.fname, dark.shape, mode='r', memmap=dark.memmap) as dark_data:
            err = dark_data[:self.ysize,:self.xsize]
            mean, std = err.mean(), err.std()
            dat += err
        self.addHistory("Adding Dark residual with mean %f and standard deviation %f" % (mean, std))
        self._log("info","Adding Dark residual with mean %f and standard deviation %f" % (mean, std))
    
    def introduceCosmicRayResidual(self, pixel_size):
        """
        Simulate CR correction error.
        
        pixel_size: pixel size in microns (on detector)
        
        returns: mean,std: float. Mean and standard deviation of cosmic ray image that was added.
        """
        energies = (600.0,5000.0) # e- (gal, SS)
        rates = (5.0,5.0) # hits/cm^2/s (gal, SS)
        pixarea = pixel_size**2 / 1E8 # cm^2
        probs = GetCrProbs(rates, pixarea, self.exptime) # hits
        cr_size, cr_psf = GetCrTemplate()

        n = os.path.join(self.out_path, uuid.uuid4().hex+"_cosmic.tmp")
        try:
            with ImageData(self.fname, self.shape, mode='r+', memmap=self.memmap) as dat:
                t_shape = tuple(self.shape)
                if self.memmap:
                    noise_data = np.memmap(n, dtype='float32', mode='w+', shape=t_shape)
                    noise_data.fill(0.)
                else:
                    noise_data = np.zeros(t_shape, dtype='float32')
                for i in range(len(energies)):
                    noise_data += MakeCosmicRay(self.shape[1], self.shape[0], probs[i], energies[i], cr_size, cr_psf, self.seed, verbose=False)
                noise_data *= 0.01
                mean, std = noise_data.mean(), noise_data.std()
                dat += noise_data
                del noise_data
            if os.path.exists(n):
                os.remove(n)
        except Exception as e:
            if os.path.exists(n):
                os.remove(n)
            if os.path.exists(self.fname):
                os.remove(self.fname)
            raise e
        self.addHistory("Adding Cosmic Ray residual with mean %f and standard deviation %f" % (mean, std))
        self._log("info","Adding Cosmic Ray residual with mean %f and standard deviation %f" % (mean, std))

    def _prepHeader(self):
        """
        Prepare the header WCS based on instrument and filter. For now, set RA and DEC to be
        identically zero. 
        """
        self.header['DATE-OBS'] = time.strftime("%Y-%m-%d")
        self.header['TIME-OBS'] = time.strftime("%h:%m:%s")
        self.header['EQUINOX'] = 2000.0
        self.header['PA_APER'] = self.pa
        self.header['VAFACTOR'] = 0.
        self.header['ORIENTAT'] = self.pa
        self.header['RA_APER'] = self.ra
        self.header['DEC_APER'] = self.dec
        self.header['NAXIS1'] = self.xsize
        self.header['NAXIS2'] = self.ysize

    def _prepRaDec(self):
        """Figure out the index of RA and DEC in the FITS image"""
        if self.wcs.wcs.lngtyp == 'RA':
            self.ranum = self.wcs.wcs.lng
            self.decnum = self.wcs.wcs.lat
        else:
            self.ranum = self.wcs.wcs.lat
            self.decnum = self.wcs.wcs.lng
    
    def _getSip(self,**kwargs):
        """Determines if a SIP polynomial is in the input arguments"""
        distortion = None
        sip = None
        distortion = kwargs.get('distortion', None)
        if distortion is not None:
            if 'sip' in distortion:
                sip = distortion['sip']
            else:
                if 'DIST_A' not in distortion:
                    raise ValueError("Can't do co-ordinate distortion without DIST_A array")
                dist_a = distortion['DIST_A']
                if 'DIST_B' not in distortion:
                    raise ValueError("Can't do co-ordinate distortion without DIST_B array")
                dist_b = distortion['DIST_B']
                dist_ap = None
                if 'DIST_AP' in kwargs:
                    dist_ap = distortion['DIST_AP']
                dist_bp = None
                if 'DIST_BP' in kwargs:
                    dist_bp = distortion['DIST_BP']
                sip = self._sip(dist_a,dist_b,dist_ap,dist_bp)
                self._log('info', "SIP A Array: {}".format(sip.a))
                self._log('info', "SIP B Array: {}".format(sip.b))
        return sip

    def _getWcs(self,**kwargs):
        """Makes a WCS from the input arguments"""
        sip = self._getSip(**kwargs)
        if 'wcs' in kwargs:
            wcs = kwargs['wcs']
            if sip is not None and wcs.sip is None:
                wcs.sip = sip
            wcs = self._normalizeWCS(wcs)
        else:
            #get pixel scale (if available)
            scale = self._scale
            ra = kwargs.get('ra', 0.)
            offset_ra = kwargs.get('offset_ra', 0.)
            dec = kwargs.get('dec', 0.)
            offset_dec = kwargs.get('offset_dec', 0.)
            ra,dec = OffsetPosition(ra, dec, offset_ra, offset_dec)
            pa = kwargs.get('pa', 0.)
            wcs = self._wcs(ra, dec, pa, scale, sip=sip)
        return wcs
       
    def _sip(self, da, db, dap, dbp):
        """Create a SIP distortion model from the distortion arrays"""
        crpix = [self.xsize//2, self.ysize//2]
        sip = wcs.Sip(da, db, dap, dbp, crpix)
        return sip
    
    def _wcs(self, ra, dec, pa, scale, crpix=None, sip=None, ranum=0, decnum=1):
        """
        Create a WCS object given the scene information.
        
        ra is right ascension in decimal degrees
        dec is declination in decimal degrees
        pa is the angle between north and east on the tangent plane, with the
            quadrant between northand east being positive
        scale is the pixel scale in arcseconds/pixel, for the (x, y) axes of the 
            image
        crpix is the (x, y) location on the image of the reference pixel 
            (default is centre)
        sip is the simple imaging polynomial distortion object (if present)
        ranum indicates which image axis represents RA (default 0)
        decnum indicates which image axis represents DEC (default 1)
        """
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["",""]
        w.wcs.ctype[ranum] = "RA---TAN"
        w.wcs.ctype[decnum] = "DEC--TAN"
        if crpix is None:
            w.wcs.crpix = [self.xsize//2, self.ysize//2]
        else:
            w.wcs.crpix = crpix
        w.wcs.crval = [0.,0.]
        w.wcs.crval[ranum] = ra
        w.wcs.crval[decnum] = dec
        w.wcs.cdelt = [abs(scale[0])/3600., abs(scale[1])/3600.]
        # For whatever reason, astropy hates CD matrices, and if you give
        #   astropy.wcs a CD matrix, it turns it into a PC matrix *without*
        #   rescaling it, and sets the CDELT keywords to 1.0, which is not
        #   actually a valid FITS WCS format. As such, in order to make
        #   a WCS that astropy won't turn invalid, we need to make a PC matrix
        #   rather than a CD matrix (as was previously done).
        pc_11 = np.cos(np.radians(pa))
        pc_12 = -np.sin(np.radians(pa))
        pc_21 = np.sin(np.radians(pa))
        pc_22 = np.cos(np.radians(pa))
        w.wcs.pc = [[pc_11, pc_12], [pc_21, pc_22]]
        if sip is not None:
            w.wcs.ctype[ranum] = "RA---TAN-SIP"
            w.wcs.ctype[decnum] = "DEC--TAN-SIP"
            w.sip = sip
        self._scale = scale
        msg = "{}: (RA, DEC, PA) := ({}, {}, {}), detected as ({}, {}, {})"
        msg = msg.format(self.name, ra, dec, pa, w.wcs.crval[ranum], 
                         w.wcs.crval[decnum], self._getPA(w, scale, decnum))
        self._log("info", msg)
        return w
    
    def _normalizeWCS(self,w):
        if w.wcs.lngtyp == 'RA' and w.wcs.lattyp == 'DEC':
            ranum = w.wcs.lng
            decnum = w.wcs.lat
        elif w.wcs.lattyp == 'RA' and w.wcs.lngtyp == 'DEC':
            ranum = w.wcs.lat
            decnum = w.wcs.lng
        elif w.wcs.lngtyp.strip() == '' and w.wcs.lattyp.strip() == '':
            ranum = 0
            decnum = 1
        else:
            raise ValueError("Lattype = %s and Lngtype = %s: can't get RA, DEC" % (w.wcs.lngtyp,w.wcs.lattyp))
        ra = w.wcs.crval[ranum]
        dec = w.wcs.crval[decnum]
        if w.wcs.has_cd():
            scale = wcs.utils.proj_plane_pixel_scales(w)*3600.
        else:
            scale = w.wcs.cdelt*3600.
        pa = self._getPA(w,scale,decnum)
        return self._wcs(ra,dec,pa,scale,w.wcs.crpix,w.sip,ranum,decnum)
    
    def _getPA(self,wcs,scale,decnum=1):
        """With non-normalized WCS, find the angle of the pixel offset of increasing DEC by 1 arcsec"""
        offset_value = np.array([0.,0.])
        offset_value[decnum] += 1./3600.
        world_coords = np.array((wcs.wcs.crval,wcs.wcs.crval+offset_value))
        pix_coords = wcs.wcs_world2pix(world_coords, 1)
        offset = (pix_coords[1] - pix_coords[0])
#         offset = offset / scale # Divide by scale to correct for non-square pixels
        pa_north = np.degrees(np.arctan2(offset[0],offset[1])) % 360. % 360.
        return pa_north

    def __add__(self,other):
        """Adds one AstroImage to another."""
        result = self.copy()
        result.addWithAlignment(other)
        for k,v in other.header.iteritems():
            if k not in result.header:
                result.header[k] = v
        for item in other.history:
            result.history.append("%s: %s" (other.name,item))
        result.history.append("Taken from %s and added %s" (self.name,other.name))
        return result

    def __iadd__(self,other):
        """Adds an AstroImage to the current one. (i.e. the '+=' operator)"""
        if isinstance(other, int) or isinstance(other, float):
            self.addHistory("Added constant %f/pixel" % (other))
            with ImageData(self.fname, self.shape, memmap=self.memmap) as dat:
                dat += other
            return self
        else:
            self.addWithAlignment(other)
            for k,v in other.header.iteritems():
                if k not in self.header:
                    self.header[k] = v
            for item in other.history:
                self.history.append("%s: %s" % (other.name,item))
            self.history.append("Added %s" (other.name))
            return self

    def __radd__(self,other):
        """
        Adds an integer or floating-point constant to the AstroImage
        
        .. warning:: Assumes constant value per-pixel
        """
        self.addHistory("Added %f/pixel" % (other))
        with ImageData(self.fname, self.shape, memmap=self.memmap) as dat:
            dat += other
        return self
    
    def __mul__(self, other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        result = self.copy()
        with ImageData(result.fname, result.shape, memmap=result.memmap) as dat:
            dat *= other
        result.addHistory("Multiplied by %f" % (other))
        return result

    def __imul__(self,other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        self.addHistory("Multiplied by %f" % (other))
        with ImageData(self.fname, self.shape, memmap=self.memmap) as dat:
            dat *= other
        return self

    def __rmul__(self,other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        result = self.copy()
        with ImageData(result.fname, result.shape, memmap=result.memmap) as dat:
            dat *= other
        result.addHistory("Multiplied by %f" % (other))
        return result

    @property
    def sum(self):
        """Returns the sum of the flux in the current image"""
        with ImageData(self.fname, self.shape, memmap=self.memmap) as dat:
            sum = np.sum(dat)
        return sum

    def _log(self,mtype,message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self,'logger'):
            getattr(self.logger,mtype)(message)
        else:
            sys.stderr.write("%s: %s\n" % (mtype,message))
    
    def _init_dat(self, base_shape, psf_shape=(0,0), data=None):
        self.base_shape = base_shape
        self.shape = np.array(base_shape) + np.array(psf_shape)
        t_shape = tuple(self.shape)
        if self.memmap:
            if os.path.exists(self.fname):
                os.remove(self.fname)
            fp = np.memmap(self.fname, dtype='float32', mode='w+', 
                           shape=t_shape)
            fp.fill(0.)
        else:
            fp = np.zeros(t_shape, dtype='float32')
            self.fname = fp
        if data is not None:
            centre = self.shape//2
            cy, cx = centre
            half = base_shape//2
            hy, hx = half
            fp[cy-hy:cy+base_shape[0]-hy, cx-hx:cx+base_shape[1]-hx] = data
        if self.memmap:
            del fp
    
    def _remap(self, xs, ys):
        # Step 1 -- compensate for PSF adjustments
        adj = ((self.shape - self.base_shape)/2).astype(np.int32)
        out_y = (ys - adj[0])/self.oversample
        out_x = (xs - adj[1])/self.oversample
        return out_x, out_y

    def updateState(self, state):
        if self.set_celery is not None:
            self.set_celery(state)
    
    def getState(self):
        if self.get_celery is not None:
            return self.get_celery()
        return ""
    
    
    INSTRUMENT_DEFAULT = {
                            'telescope': 'wfirst',
                            'instrument': 'WFI',
                            'filter': 'F062',
                            'detector': {
                                            'WFI': 'SCA01',
                                            'NIRCam': 'A1',
                                            'MIRI': 'MIRI'
                                        },
                            'shape': (4096, 4096),
                            'scale': [0.11,0.11],
                            'zeropoint': 21.0,
                            'photflam': 0.,
                            'photplam': 0.6700,
                            'background': 0.,
                            'oversample': 1,
                            'psf_grid_size': 1
                         }

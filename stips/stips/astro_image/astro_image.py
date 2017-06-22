from __future__ import absolute_import,division
__filetype__ = "base"

#External Modules
import inspect, logging, os, pytest, scipy, shutil, sys, time, tempfile, uuid

import numpy as np

from astropy import wcs
from astropy.convolution import convolve_fft
from astropy.io import fits as pyfits
from astropy.table import Table,Column
from photutils import CircularAperture, aperture_photometry

from copy import deepcopy

from cStringIO import StringIO

from scipy.ndimage.interpolation import zoom,rotate
from scipy.signal import fftconvolve
from scipy.fftpack import fftn,ifftn

#Local Modules
from ..utilities import OffsetPosition, overlapadd2, read_table
from ..errors import GetCrProbs, GetCrTemplate, MakeCosmicRay
from ..galaxy_module import Sersic

class ImageData(object):
    def __init__(self, fname, shape, mode='r+'):
        self.fp = np.memmap(fname, dtype='float32', mode=mode, shape=shape)
    
    def __enter__(self):
        return self.fp
    
    def __exit__(self, exc_type, exc_value, traceback):
        del self.fp


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
        
        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            stream_handler = logging.StreamHandler(sys.stderr)
            stream_handler.setLevel(logging.DEBUG)
            stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
            self.logger = logging.getLogger()
            self.logger.addHandler(stream_handler)
    
        #Set unique ID and figure out where the numpy memmap will be stored
        self.out_path = kwargs.get('out_path', os.getcwd())
        self.name = kwargs.get('detname', "")
        
        self.oversample = kwargs.get('oversample', 1)
        psf_shape = kwargs.get('psf_shape', (0, 0))
        data = kwargs.get('data', None)
        if data is not None:
            base_shape = data.shape
        else:
            base_shape = kwargs.get('shape', (1, 1))
        self._init_dat(base_shape, psf_shape, data)
        
        #Get WCS values if present, or set up a default
        self.wcs = self._getWcs(**kwargs)
        self._prepRaDec()
        
        #Header
        self.header = kwargs.get('header', {})
        if self.header == {}:
            self.header = kwargs.get('imh', {})
        
        self._prepHeader()
        
        #History
        self.history = kwargs.get('history', [])
        
        #Zero Point. Necessary for output catalogues.
        self.zeropoint = kwargs.get('zeropoint', 0.)
        self.photflam = kwargs.get('photflam', 0.)
    
    def __del__(self):
        if os.path.exists(self.fname):
            os.remove(self.fname)

    def copy(self):
        other = AstroImage(out_path=self.out_path, detname=self.name, wcs=self.wcs, header=self.header, history=self.history,
                           xsize=self.xsize, ysize=self.ysize, zeropoint=self.zeropoint, photflam=self.photflam, 
                           logger=self.logger)
        with tempfile.NamedTemporaryFile(delete=False) as f:
            shutil.copy(self.fname, f.name)
            if os.path.exists(other.fname):
                os.remove(other.fname)
            other.fname = f.name
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
                with pyfits.open(file) as fits:
                    ext = kwargs.get('ext', 0)
                    dat = fits[ext].data
                    img._init_dat(base_shape=dat.shape, data=dat)
                    my_wcs = wcs.WCS(fits[ext].header)
                    for k,v in fits[ext].header.items():
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
                with pyfits.open(file) as fits:
                    ext = kwargs.get('ext', 0)
                    dat = fits[ext].data
                    img._init_dat(base_shape=dat.shape, data=dat)
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
        self.shape = (self.ysize, size)
    
    @property
    def ysize(self):
        return self.shape[0]
    
    @ysize.setter
    def ysize(self, size, offset=0):
        """Change the vertical size. The offset will be applied to the new image before adding"""
        self.crop(self.xsize, size, 0, offset)
        self.shape = (size, self.xsize)
    
    @property
    def xscale(self):
        return abs(self.wcs.wcs.cdelt[0])*3600.
    
    @property
    def yscale(self):
        return abs(self.wcs.wcs.cdelt[1])*3600.
    
    @property
    def scale(self):
        return [self.xscale, self.yscale]
    
    @property
    def rascale(self):
        return abs(self.wcs.wcs.cdelt[self.ranum])*3600.

    @property
    def decscale(self):
        return abs(self.wcs.wcs.cdelt[self.decnum])
    
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
    def hdu(self):
        """Output AstroImage as a FITS Primary HDU"""
        self._log("info","Creating Primary HDU from AstroImage %s" % (self.name))
        with ImageData(self.fname, self.shape, mode='r+') as dat:
            hdu = pyfits.PrimaryHDU(dat, header=self.wcs.to_header(relax=True))
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
        with ImageData(self.fname, self.shape, mode='r+') as dat:
            hdu = pyfits.ImageHDU(dat, header=self.wcs.to_header(relax=True), name=self.name)
        for k, v in self.header.iteritems():
            hdu.header[k] = v
        for item in self.history:
            hdu.header.add_history(item)
        self._log("info","Created Extension HDU from AstroImage %s" % (self.name))
        return hdu
    
    def toFits(self, outFile):
        """Create a FITS file from the current state of the AstroImage data."""
        self._log("info","Writing AstroImage %s to FITS" % (self.name))
        hdulist = pyfits.HDUList([self.hdu])
        hdulist.writeto(outFile, clobber=True)
    
    def updateHeader(self,k,v):
        """
        Updates a single keyword in the header dictionary, replacing the current value if there is
        one, otherwise adding a new value
        """
        self.header[k] = v
    
    def addHistory(self,v):
        """Adds an entry to the header history list."""
        self.history.append(v)
    
    def addTable(self, t, dist=False):
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
        to_keep = np.where((xs > 0) & (xs <= self.xsize) & (ys > 0) & (ys <= self.ysize))
        self._log("info", "Keeping {} items".format(len(xs[to_keep])))
        ot = None
        if len(xs[to_keep]) > 0:
            xs = xs[to_keep]
            ys = ys[to_keep]
            xfs, yfs = self.remap(xs, ys)
            fluxes = t['flux'][to_keep]
            types = t['type'][to_keep]
            ns = t['n'][to_keep]
            res = t['re'][to_keep]
            phis = t['phi'][to_keep]
            ratios = t['ratio'][to_keep]
            ids = t['id'][to_keep]
            old_notes = t['notes'][to_keep]
            notes = np.empty_like(xs, dtype="S6")
            notes = np.where(old_notes != "", old_notes,"")
            vegamags = -2.512 * np.log10(fluxes) - self.zeropoint
            stmags = -2.5 * np.log10(fluxes * self.photflam) - 21.10
            stars_idx = np.where(types == 'point')
            if len(xs[stars_idx]) > 0:
                self._log("info", "Writing {} stars".format(len(xs[stars_idx])))
                self.addPoints(xs[stars_idx], ys[stars_idx], fluxes[stars_idx])
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
                counter = 1
                total = len(gxs)
                for (x, y, flux, n, re, phi, ratio) in zip(gxs, gys, gfluxes, gns, gres, gphis, gratios):
                    self.addSersicProfile(x, y, flux, n, re, phi, ratio)
                    self._log("info", "Finished Galaxy {} of {}".format(counter, total))
                    counter += 1
            ot = Table()
            ot['x'] = Column(data=xfs, unit='pixels')
            ot['y'] = Column(data=yfs, unit='pixels')
            ot['type'] = Column(data=types)
            ot['vegamag'] = Column(data=vegamags)
            ot['stmag'] = Column(data=stmags)
            ot['countrate'] = Column(data=fluxes, unit='counts/s')
            ot['id'] = Column(data=ids)
            ot['notes'] = Column(data=notes)
        return ot
    
    def addCatalogue(self, cat, dist=False):
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
        self._log("info","Adding catalogue %s to AstroImage %s" % (catname, self.name))
        obsname = os.path.join(self.out_path, os.path.splitext(catname)[0]+"_observed_%s.txt" % (self.name))
        self.addHistory("Adding items from catalogue %s" % (cat))
        data = None
        with open(obsname, 'w') as outf:
            for i, t in enumerate(read_table(cat)):
                ot = self.addTable(t, dist)
                if ot is not None:
                    data = StringIO()
                    ot.write(data, format='ascii.ipac')
                    data.seek(0)
                    if i != 0: # get rid of the header lines
                        data.readline()
                        data.readline()
                        data.readline()
                        data.readline()
                    outf.write(data.read())
            if data is None:
                outf.write("No sources from catalogue were visible.\n")
        self._log("info","Added catalogue %s to AstroImage %s" % (catname, self.name))
        return obsname            
    
    def addPoints(self, xs, ys, rates):
        """Adds a set of point sources to the image given their co-ordinates and count rates."""
        self.addHistory("Adding %d point sources" % (len(xs)))
        self._log("info","Adding %d point sources to AstroImage %s" % (len(xs),self.name))
        xs = np.floor(xs).astype(int)
        ys = np.floor(ys).astype(int)
        with ImageData(self.fname, self.shape) as dat:
            dat[ys, xs] += rates
    
    def addSersicProfile(self,posX,posY,flux,n,re,phi,axialRatio):
        """
        Adds a single sersic profile to the image given its co-ordinates, count rate, and source 
        type.
        
        (posX,poxY) are the co-ordinates of the centre of the profile (pixels).
        flux is the total number of counts to add to the AstroImage.
        n is the Sersic profile index.
        re is the radius enclosing half of the total light (pixels).
        phi is the angle of the major axis (degrees east of north).
        axialRatio is the ratio of major axis to minor axis.
        """
        self.addHistory("Adding Sersic profile at (%f,%f) with flux %f, index %f, Re %f, Phi %f, and axial ratio %f" % (posX,posY,flux,n,re,phi,axialRatio))
        self._log("info","Adding sersic profile to AstroImage %s" % (self.name))
        from astropy.modeling.models import Sersic2D
        x, y, = np.meshgrid(np.arange(self.xsize), np.arange(self.ysize))
        mod = Sersic2D(amplitude=flux, r_eff=re, n=n, x_0=posX, y_0=posY, ellip=axialRatio, theta=(np.radians(phi) + 0.5*np.pi))
        with ImageData(self.fname, self.shape) as dat:
            img = mod(x, y)
            aperture = CircularAperture((posX, posY), re)
            flux_table = aperture_photometry(img, aperture)
            central_flux = flux_table['aperture_sum'][0]
            if central_flux != 0.:
                factor = flux / central_flux
                img *= factor
                new_flux_table = aperture_photometry(img, aperture)
                new_central_flux = new_flux_table['aperture_sum'][0]
                self._log("info", "Flux within half-light radius is {} ({} input)".format(new_central_flux, flux))
            else:
                self._log("info", "Flux within half-light radius is {} ({} input)".format(central_flux, flux))
            dat += img
        
#         s = Sersic(px,py,n,xs=self.xsize,ys=self.ysize,flux=flux,q=axialRatio,phi=p,re=re,lf=self._log)
#         with ImageData(self.fname, self.shape) as dat:
#             dat += s.image

    def convolve(self,other):
        """Convolves the AstroImage with another (provided) AstroImage, e.g. for PSF convolution."""
        self.addHistory("Convolving with file %s" % (other.name))
        self._log("info","Convolving AstroImage %s with %s" % (self.name,other.name))
        with tempfile.NamedTemporaryFile() as f, tempfile.NamedTemporaryFile(delete=False) as g:
            with ImageData(self.fname, self.shape, mode='r') as dat, ImageData(other.fname, other.shape, mode='r') as psf:
                fp_result = np.memmap(f.name, dtype='float32', mode='w+', shape=(self.shape[0]+psf.shape[0]-1, self.shape[1]+psf.shape[1]-1))
                sub_shape = (min(4095 - psf.shape[0], self.shape[0] + psf.shape[0] - 1), min(4095 - psf.shape[1], self.shape[1] + psf.shape[1] - 1))
                self._log('info', "PSF Shape: {}; Current Shape: {}".format(psf.shape, self.shape))
                self._log('info', "Choosing between 4095 - {} = {} and {} + {} - 1 = {}".format(psf.shape, 4095 - psf.shape[0],
                                                                                                psf.shape, self.shape,
                                                                                                psf.shape[0] + self.shape[0] - 1))
                self._log('info', "Using overlapping arrays of size {}".format(sub_shape))
                overlapadd2(dat, psf, sub_shape, y=fp_result)
                self._log('info', "Cropping convolved image down to detector size")
                centre = (fp_result.shape[0]//2, fp_result.shape[1]//2)
                half = (self.base_shape[0]//2, self.base_shape[1]//2)
                self._log('info', "Image Centre: {}; Image Half-size: {}".format(centre, half))
                self._log('info', "Taking [{}:{}, {}:{}]".format(centre[0]-half[0], centre[0]+half[0], centre[1]-half[1], centre[1]+half[1]))
                fp_crop = np.memmap(g.name, dtype='float32', mode='w+', shape=self.base_shape)
                fp_crop[:,:] = fp_result[centre[0]-half[0]:centre[0]+half[0], centre[1]-half[1]:centre[1]+half[1]]
                crpix = [half[0], half[1]]
                if self.wcs.sip is not None:
                    sip = wcs.Sip(self.wcs.sip.a, self.wcs.sip.b, None, None, crpix)
                else:
                    sip = None
                self.wcs = self._wcs(self.ra, self.dec, self.pa, self.scale, crpix=crpix, sip=sip)
                del fp_result
                del fp_crop
            if os.path.exists(self.fname):
                os.remove(self.fname)
            self.fname = g.name
            self.shape = self.base_shape

    def rotate(self,angle,reshape=False):
        """
        Rotate the image a number of radians as specified
        
        ..warning:: This function is not necessarily flux-conserving
        """
        self.addHistory("Rotating by %f degrees" % (angle))
        self._log("info","Rotating AstroImage %s by %f degrees" % (self.name,angle))
        self.pa = (self.pa + angle)%360.%360.
        with tempfile.NamedTemporaryFile(delete=False) as f:
            fp_result = np.memmap(f.name, dtype='float32', mode='w+', shape=self.shape)
            with ImageData(self.fname, self.shape, mode='r') as dat:
                rotate(dat, angle, order=5, reshape=reshape, output=fp_result)
            del fp_result
            if os.path.exists(self.fname):
                os.remove(self.fname)
            self.fname = f.name

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
    
    def _addWithOffset(self, other, offset_x, offset_y):
        """
        Add a np array to the current image with a given pixel offset. This allows for a smaller
        memory use for addWithAlignment.
        """
        (ys,xs) = other.shape
        yc = int(np.floor(ys/2))
        xc = int(np.floor(xs/2))
        
        ycenter = int(np.floor(self.ysize/2))
        xcenter = int(np.floor(self.xsize/2))
        
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
        
        #If any of these are false, the images are disjoint.
        if low_x < self.xsize and low_y < self.ysize and high_x > 0 and high_y > 0:
            with ImageData(self.fname, self.shape, mode='r+') as dat, ImageData(other.fname, other.shape, mode='r') as other_data:
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
        with tempfile.NamedTemporaryFile(delete=False) as f:
            shape_x = int(round(self.shape[1] * self.scale[0] / scale[0]))
            shape_y = int(round(self.shape[0] * self.scale[1] / scale[1]))
            new_shape = (shape_y, shape_x)
            self._log("info","New shape will be {}".format(new_shape))
            fp_result = np.memmap(f.name, dtype='float32', mode='w+', shape=new_shape)
            with ImageData(self.fname, self.shape, mode='r+') as dat:
                flux = dat.sum()
                zoom(dat, (np.array(self.scale)/np.array(scale)), fp_result)
            factor = flux / fp_result.sum()
            fp_result *= factor
            del fp_result
            if os.path.exists(self.fname):
                os.remove(self.fname)
            self.fname = f.name
            self.shape = new_shape
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
        self._log('info', "Pre-bin: (RA, DEC, PA) = ({}, {}, {})".format(self.ra, self.dec, self.pa))
        self.rescale([self.scale[0]*binx,self.scale[1]*biny])
        self.oversample = 1
        self._log('info', "Post-bin: (RA, DEC, PA) = ({}, {}, {})".format(self.ra, self.dec, self.pa))

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

        with tempfile.NamedTemporaryFile() as a, tempfile.NamedTemporaryFile() as n, ImageData(self.fname, self.shape, mode='r+') as dat:
            abs_data = np.memmap(a.name, dtype='float32', mode='w+', shape=self.shape)
            np.absolute(dat, abs_data)

            noise_data = np.memmap(n.name, dtype='float32', mode='w+', shape=self.shape)
            noise_data[:,:] = np.random.normal(size=self.shape) * np.sqrt(abs_data)
            del abs_data
            if absVal:
                noise_data[:,:] = np.abs(noise_data)
            mean, std = noise_data.mean(), noise_data.std()
            dat += noise_data
            del noise_data
        self.addHistory("Adding Poisson Noise with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding Poisson Noise with mean {} and standard deviation {}".format(mean, std))            

    def introduceReadnoise(self,readnoise):
        """
        Simulate kTC or read noise for detector.

        Parameters
        ----------
        readnoise: constant representing average read noise per pixel.
        """
        with tempfile.NamedTemporaryFile() as n, ImageData(self.fname, self.shape, mode='r+') as dat:
            noise_data = np.memmap(n.name, dtype='float32', mode='w+', shape=self.shape)
            noise_data[:,:] = readnoise * np.random.randn(self.ysize,self.xsize)            
            mean, std = noise_data.mean(), noise_data.std()
            dat += noise_data
            del noise_data
        self.addHistory("Adding Read noise with mean %f and standard deviation %f" % (mean, std))
        self._log("info","Adding readnoise with mean %f and STDEV %f" % (mean, std))

    def introduceFlatfieldResidual(self,flat):
        """
        Simulate flatfield correction error.
        
        flat: AstroImage containing flatfield error values.
        
        returns: mean, std. Mean and standard deviation of used portion of error array.
        """
        with ImageData(self.fname, self.shape, mode='r+') as dat, ImageData(flat.fname, flat.shape, mode='r') as flat_data:
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
        with ImageData(self.fname, self.shape, mode='r+') as dat, ImageData(dark.fname, dark.shape, mode='r') as dark_data:
            err = dark_data[:self.ysize,:self.xsize]
            mean, std = err.mean(), err.std()
            dat += err
        self.addHistory("Adding Dark residual with mean %f and standard deviation %f" % (mean, std))
        self._log("info","Adding Dark residual with mean %f and standard deviation %f" % (mean, std))
    
    def introduceCosmicRayResidual(self,pixel_size,exptime):
        """
        Simulate CR correction error.
        
        pixel_size: pixel size in microns (on detector)
        
        exptime: exposure time for this image.
        
        returns: mean,std: float. Mean and standard deviation of cosmic ray image that was added.
        """
        energies = (600.0,5000.0) # e- (gal, SS)
        rates = (5.0,5.0) # hits/cm^2/s (gal, SS)
        pixarea = pixel_size**2 / 1E8 # cm^2
        probs = GetCrProbs(rates, pixarea, exptime) # hits
        cr_size, cr_psf = GetCrTemplate()

        with tempfile.NamedTemporaryFile() as n, ImageData(self.fname, self.shape, mode='r+') as dat:
            noise_data = np.memmap(n.name, dtype='float32', mode='w+', shape=self.shape)
            noise_data.fill(0.)
            for i in range(len(energies)):
                noise_data += MakeCosmicRay(self.shape[1], self.shape[0], probs[i], energies[i], cr_size, cr_psf, verbose=False)
            noise_data *= 0.01
            mean, std = noise_data.mean(), noise_data.std()
            dat += noise_data
            del noise_data
        self.updateHeader('exptime', exptime)
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
            scale = kwargs.get('scale', [1., 1.])
            ra = kwargs.get('ra', 0.)
            offset_ra = kwargs.get('offset_ra', 0.)
            dec = kwargs.get('dec', 0.)
            offset_dec = kwargs.get('offset_dec', 0.)
            ra,dec = OffsetPosition(ra,dec,offset_ra,offset_dec)
            pa = kwargs.get('pa', 0.)
            wcs = self._wcs(ra,dec,pa,scale,sip=sip)
        return wcs
       
    def _sip(self,da,db,dap,dbp):
        """Create a SIP distortion model from the distortion arrays"""
        crpix = [int(np.floor(self.xsize/2.)),int(np.floor(self.ysize/2.))]
        sip = wcs.Sip(da,db,dap,dbp,crpix)
        return sip
    
    def _wcs(self,ra,dec,pa,scale,crpix=None,sip=None,ranum=0,decnum=1):
        """Create a WCS object given the scene information."""
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["",""]
        w.wcs.ctype[ranum] = "RA---TAN"
        w.wcs.ctype[decnum] = "DEC--TAN"
        if crpix is None:
            w.wcs.crpix = [int(np.floor(self.xsize/2.)),int(np.floor(self.ysize/2.))]
        else:
            w.wcs.crpix = crpix
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        self._log('info', "From {}: Detector {}: Detector Size: {}, CRPIX: {}".format(calframe[1][3], self.name, self.shape, w.wcs.crpix))
        w.wcs.crval = [0.,0.]
        w.wcs.crval[ranum] = ra
        w.wcs.crval[decnum] = dec
        scale_list = [-abs(scale[0])/3600., abs(scale[1]/3600.)]
        cdelt_array = np.array(scale_list)
        w.wcs.cdelt = cdelt_array
        cpa = np.cos(np.radians(pa))
        spa = np.sin(np.radians(pa))
        w.wcs.pc = np.array([[cpa, -spa], [-spa, cpa]])
        if sip is not None:
            w.wcs.ctype[ranum] = "RA---TAN-SIP"
            w.wcs.ctype[decnum] = "DEC--TAN-SIP"
            w.sip = sip
        self._log("info", "{}: (RA, DEC, PA) set to ({}, {}, {}), found to be ({}, {}, {})".format(self.name, ra, dec, pa,
                                                                                                   w.wcs.crval[ranum],
                                                                                                   w.wcs.crval[decnum],
                                                                                                   self._getPA(w, scale)))
        return w
    
    def _normalizeWCS(self,w):
        if w.wcs.lngtyp == 'RA' and w.wcs.lattyp == 'DEC':
            ranum = w.wcs.lng
            decnum = w.wcs.lat
        elif w.wcs.lattyp == 'RA' and w.wcs.lngtyp == 'DEC':
            ranum = w.wcs.lat
            decnum = w.wcs.lng
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
            with ImageData(self.fname, self.shape) as dat:
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
        with ImageData(self.fname, self.shape) as dat:
            dat += other
        return self
    
    def __mul__(self, other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        result = self.copy()
        with ImageData(result.fname, result.shape) as dat:
            dat *= other
        result.addHistory("Multiplied by %f" % (other))
        return result

    def __imul__(self,other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        self.addHistory("Multiplied by %f" % (other))
        with ImageData(self.fname, self.shape) as dat:
            dat *= other
        return self

    def __rmul__(self,other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        result = self.copy()
        with ImageData(result.fname, result.shape) as dat:
            dat *= other
        result.addHistory("Multiplied by %f" % (other))
        return result

    def _log(self,mtype,message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self,'logger'):
            getattr(self.logger,mtype)(message)
        else:
            sys.stderr.write("%s: %s\n" % (mtype,message))
    
    def _init_dat(self, base_shape, psf_shape=(0,0), data=None):
        if hasattr(self, 'fname') and self.fname is not None and os.path.exists(self.fname):
            os.remove(self.fname)
        with tempfile.NamedTemporaryFile(delete=False) as f:
            self.fname = f.name
            self.base_shape = base_shape
            self.shape = tuple(np.array(base_shape) + np.array(psf_shape))
            fp = np.memmap(self.fname, dtype='float32', mode='w+', shape=self.shape)
            fp.fill(0.)
            if data is not None:
                centre = tuple(np.array(self.shape)//2)
                half = tuple(np.array(base_shape)//2)
                fp[centre[0]-half[0]:centre[0]+self.base_shape[0]-half[0],centre[1]-half[1]:centre[1]+self.base_shape[1]-half[1]] = data
            del fp
    
    def remap(self, xs, ys):
        # Step 1 -- compensate for PSF adjustments
        adj_x = (self.shape[1] - self.base_shape[1]) // 2
        adj_y = (self.shape[0] - self.base_shape[0]) // 2
        x_outs = xs - adj_x
        y_outs = ys - adj_y
        # Step 2 -- handle oversample
        x_outs /= self.oversample
        y_outs /= self.oversample
        return x_outs, y_outs

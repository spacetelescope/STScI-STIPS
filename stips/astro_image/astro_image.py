__filetype__ = "base"

# External Modules
from astropy import units as u
from astropy import wcs
from astropy.convolution import convolve_fft
from astropy.io import fits
from astropy.table import Table, Column
from copy import deepcopy
import logging
import numpy as np
import os
from pandeia.engine import coords, profile, source
from scipy.ndimage.interpolation import zoom, rotate
import sys
import time

# Local Modules
from ..errors import GetCrProbs, GetCrTemplate, MakeCosmicRay
from ..utilities import StipsEnvironment
from ..utilities import OffsetPosition
from ..utilities import StipsDataTable
from ..utilities import SelectParameter
from ..utilities import sersic_lum
from ..utilities import rind
# PSF making functions
from ..utilities.makePSF import interpolate_epsf, make_epsf, place_source
# PSF constants
from ..utilities.makePSF import PSF_BOXSIZE, PSF_BRIGHT_BOXSIZE, PSF_EXTRA_BRIGHT_BOXSIZE, PSF_GRID_SIZE, PSF_UPSCALE

stips_version = StipsEnvironment.__stips__version__


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
        if kwargs.get('psf', True):
            self.psf_shape = (PSF_BOXSIZE, PSF_BOXSIZE)
            self.has_psf = True
            self.psf_grid_size = PSF_GRID_SIZE
        else:
            self.has_psf = False
            self.psf_shape = (0, 0)
            self.psf_grid_size = 0

        if 'parent' in kwargs:
            self.parent = kwargs['parent']
            self.logger = self.parent.logger
            self.out_path = self.parent.out_path
            self.prefix = self.parent.prefix
            self.seed = self.parent.seed
            self.fast_galaxy = self.parent.fast_galaxy
            self.convolve_galaxy = self.parent.convolve_galaxy
            self.telescope = self.parent.TELESCOPE.lower()
            self.instrument = self.parent.PSF_INSTRUMENT
            self.filter = self.parent.filter
            self.bright_limit = self.parent.bright_limit
            self.xbright_limit = self.parent.xbright_limit
            self.shape = self.parent.DETECTOR_SIZE
            self._scale = self.parent.SCALE
            self.zeropoint = self.parent.zeropoint
            self.photflam = self.parent.photflam
            self.photplam = self.parent.PHOTPLAM[self.filter]
            background = self.parent.background
            self.cat_type = self.parent.cat_type
        else:
            self.parent = None
            if 'logger' in kwargs:
                self.logger = kwargs['logger']
            else:
                self.logger = logging.getLogger('__stips__')
                log_level = SelectParameter("log_level")
                self.logger.setLevel(getattr(logging, log_level))
                if not len(self.logger.handlers):
                    stream_handler = logging.StreamHandler(sys.stderr)
                    format = '%(asctime)s %(levelname)s: %(message)s'
                    stream_handler.setFormatter(logging.Formatter(format))
                    self.logger.addHandler(stream_handler)
            self.out_path = SelectParameter('out_path', kwargs)
            self.bright_limit = SelectParameter('psf_bright_limit', kwargs)
            self.xbright_limit = SelectParameter('psf_xbright_limit', kwargs)
            self.fast_galaxy = SelectParameter('fast_galaxy', kwargs)
            self.convolve_galaxy = SelectParameter('convolve_galaxy', kwargs)
            self.shape = kwargs.get('shape', default['shape'])
            self._scale = kwargs.get('scale', np.array(default['scale']))
            self.prefix = kwargs.get('prefix', '')
            self.cat_type = SelectParameter('cat_type', kwargs)
            self.seed = SelectParameter('seed', kwargs)
            self.zeropoint = kwargs.get('zeropoint', default['zeropoint'])
            self.photflam = kwargs.get('photflam', default['photflam'])
            self.photplam = kwargs.get('photplam', default['photplam'])
            background = SelectParameter('background', kwargs)
            self.telescope = kwargs.get('telescope', default['telescope'])
            self.instrument = kwargs.get('instrument', default['instrument'])
            self.filter = kwargs.get('filter', default['filter'])

        # WCS origin is 0 for numpy/C and 1 for FITS/Fortran, do not change this one.
        self.wcs_origin = 0
        # Will only be used for the output catalog file, you can change this one.
        self.out_origin = 1

        self.name = kwargs.get('detname', default['detector'][self.instrument])
        self.detector = self.name

        data = kwargs.get('data', None)
        if data is not None:
            base_shape = np.array(data.shape)
        else:
            base_shape = self.shape
        self._init_dat(base_shape, self.psf_shape, data)

        # Get WCS values if present, or set up a default
        self.wcs = self._getWcs(**kwargs)
        self._prepRaDec()

        # Header
        self.header = kwargs.get('header', kwargs.get('imh', {}))
        self._prepHeader()
        if 'exptime' in self.header:
            self.exptime = self.header['exptime']
        else:
            self.exptime = kwargs.get('exptime', 1.)
        self.updateHeader('exptime', self.exptime)

        # History
        self.history = kwargs.get('history', [])

        # Special values for Sersic profile generation
        self.profile_multiplier = kwargs.get('profile_multiplier', 100.)
        self.noise_floor = max(background, kwargs.get('noise_floor', 1.))

    def copy(self):
        other = AstroImage(out_path=self.out_path, detname=self.name, wcs=self.wcs, header=self.header, history=self.history,
                           xsize=self.xsize, ysize=self.ysize, zeropoint=self.zeropoint, photflam=self.photflam,
                           logger=self.logger, data=self.data)
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
                    for k, v in inf[ext].header.items():
                        if k != '' and k not in my_wcs.wcs.to_header():
                            img.header[k] = v
                img.wcs = img._normalizeWCS(my_wcs)
                img._prepRaDec()
                img._prepHeader()
                img.updateHeader("ASTROIMAGEVALID", True)
                img.addHistory("Created from FITS file {}".format(os.path.split(file)[1]))
                img._log("info", "Created AstroImage {} from FITS file {}".format(img.name, os.path.split(file)[1]))
            except IOError as e:
                img.updateHeader("ASTROIMAGEVALID", False)
                img.addHistory("Attempted to create from invalid FITS file {}".format(os.path.split(file)[1]))
                img._log("warning", "Attempted to create AstroImage {} from invalid FITS file {}".format(img.name, os.path.split(file)[1]))
                img._log("warning", "Error is {}".format(e))
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
                img.addHistory("Data imported from FITS file {}".format(os.path.split(file)[1]))
                img._log("info", "Created AstroImage {} and imported data from FITS file {}".format(img.name, os.path.split(file)[1]))
            except IOError as e:
                img.updateHeader("ASTROIMAGEVALID", False)
                img.addHistory("Attempted to create from invalid FITS file {}".format(os.path.split(file)[1]))
                img._log("warning", "Attempted to create AstroImage {} from invalid FITS file {}".format(img.name, os.path.split(file)[1]))
                img._log("warning", "Error was {}".format(e))
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
        img._log("info", "Creating AstroImage {} from points".format(img.name))
        if img.xsize < np.ceil(np.max(xs)) + 1:
            img.xsize = np.ceil(np.max(xs)) + 1
        if img.ysize < np.ceil(np.max(ys)) + 1:
            img.ysize = np.ceil(np.max(ys)) + 1
        img.addPoints(xs, ys, rates)
        return img

    @classmethod
    def initFromProfile(cls, posX, posY, flux, n, re, phi, axialRatio, **kwargs):
        """Convenience to initialize a blank image and add a sersic profile to it"""
        img = cls(**kwargs)
        img._log("info", "Creating AstroImage {} from Sersic Profile".format(img.name))
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
            raise ValueError("WCS has longitude {} and latitude {}. Can't get RA".format(self.wcs.wcs.lngtyp, self.wcs.wcs.lattyp))

    @ra.setter
    def ra(self, ra):
        if self.wcs.wcs.lngtyp == 'RA':
            self.wcs.wcs.crval[self.wcs.wcs.lng] = ra % 360.
            self.addHistory("Set RA to {}".format(ra))
        elif self.wcs.wcs.lattyp == 'RA':
            self.wcs.wcs.crval[self.wcs.wcs.lat] = ra % 360.
            self.addHistory("Set RA to {}".format(ra))
        else:
            raise ValueError("WCS has longitude {} and latitude {}. Can't set RA".format(self.wcs.wcs.lngtyp, self.wcs.wcs.lattyp))

    @property
    def dec(self):
        if self.wcs.wcs.lngtyp == 'DEC':
            return self.wcs.wcs.crval[self.wcs.wcs.lng]
        elif self.wcs.wcs.lattyp == 'DEC':
            return self.wcs.wcs.crval[self.wcs.wcs.lat]
        else:
            raise ValueError("WCS has longitude {} and latitude {}. Can't get DEC".format(self.wcs.wcs.lngtyp, self.wcs.wcs.lattyp))

    @dec.setter
    def dec(self, dec):
        if self.wcs.wcs.lngtyp == 'DEC':
            self.wcs.wcs.crval[self.wcs.wcs.lng] = dec
            self.addHistory("Set DEC to {}".format(dec))
        elif self.wcs.wcs.lattyp == 'DEC':
            self.wcs.wcs.crval[self.wcs.wcs.lat] = dec
            self.addHistory("Set DEC to {}".format(dec))
        else:
            raise ValueError("WCS has longitude {} and latitude {}. Can't set DEC".format(self.wcs.wcs.lngtyp, self.wcs.wcs.lattyp))

    @property
    def pa(self):
        """WCS is normalized, so PA from angle will work"""
        return self._getPA(self.wcs, self.scale, self.decnum)

    @pa.setter
    def pa(self, pa):
        cpa = np.cos(np.radians(pa % 360.))
        spa = np.sin(np.radians(pa % 360.))
        self.wcs.wcs.pc = np.array([[cpa, -spa], [spa, cpa]])
        self.addHistory("Set PA to {}".format(pa))

    @property
    def hdu(self):
        """Output AstroImage as a FITS Primary HDU"""
        hdu = fits.PrimaryHDU(self.data, header=self.wcs.to_header(relax=True))
        if 'CDELT1' not in hdu.header:
            hdu.header['CDELT1'] = self.scale[0]/3600.
            hdu.header['CDELT2'] = self.scale[0]/3600.
        # Apparently astropy refuses to add the identity matrix to a header
        if ('PC1_1' not in hdu.header) and ('CD1_1' not in hdu.header):
            hdu.header['PC1_1'] = 1.
            hdu.header['PC1_2'] = 0.
            hdu.header['PC2_1'] = 0.
            hdu.header['PC2_2'] = 1.
        for k, v in self.header.items():
            if k != "ASTROIMAGEVALID":
                hdu.header[k] = v
        for item in self.history:
            hdu.header.add_history(item)
        self._log("info", "Created Primary HDU from AstroImage {}".format(self.name))
        return hdu

    @property
    def imageHdu(self):
        """Output AstroImage as a FITS Extension HDU"""
        self._log("info", "Creating Extension HDU from AstroImage {}".format(self.name))
        hdu = fits.ImageHDU(self.data, header=self.wcs.to_header(relax=True), name=self.name)
        if 'CDELT1' not in hdu.header:
            hdu.header['CDELT1'] = self.scale[0]/3600.
            hdu.header['CDELT2'] = self.scale[0]/3600.
        # Apparently astropy refuses to add the identity matrix to a header
        if ('PA1_1' not in hdu.header) and ('CD1_1' not in hdu.header):
            hdu.header['CD1_1'] = self.scale[0]/3600.
            hdu.header['CD1_2'] = 0.
            hdu.header['CD2_1'] = 0.
            hdu.header['CD2_2'] = self.scale[0]/3600.
        for k, v in self.header.items():
            hdu.header[k] = v
        for item in self.history:
            hdu.header.add_history(item)
        self._log("info", "Created Extension HDU from AstroImage {}".format(self.name))
        return hdu

    @property
    def psf_constructor(self):
        import webbpsf
        if not hasattr(webbpsf, self.telescope.lower()) and self.telescope.lower() == 'roman':
            return getattr(getattr(webbpsf, 'roman'), self.instrument)()
        if hasattr(webbpsf, self.instrument):
            return getattr(webbpsf, self.instrument)()
        return getattr(getattr(webbpsf, self.telescope), self.instrument)()

    def toFits(self, outFile):
        """Create a FITS file from the current state of the AstroImage data."""
        self._log("info", "Writing AstroImage {} to FITS".format(self.name))
        hdulist = fits.HDUList([self.hdu])
        hdulist.writeto(outFile, overwrite=True)

    def updateHeader(self, k, v):
        """
        Updates a single keyword in the header dictionary, replacing the current value if there is
        one, otherwise adding a new value
        """
        self.header[k] = v

    def addHistory(self, v):
        """Adds an entry to the header history list."""
        self.history.append(v)

    def addTable(self, t, dist=False, fast_galaxy=False, convolve_galaxy=True, *args, **kwargs):
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
        self._log("info", "Determining pixel co-ordinates")
        if dist and self.distorted:
            xs, ys = self.wcs.all_world2pix(t['ra'], t['dec'], self.wcs_origin,
                                            quiet=True, adaptive=True,
                                            detect_divergence=True)
        else:
            xs, ys = self.wcs.wcs_world2pix(t['ra'], t['dec'], self.wcs_origin)
        to_keep = np.where((xs >= 0) & (xs <= self.xsize) & (ys >= 0) & (ys <= self.ysize))
        self._log("info", "Keeping {} items".format(len(xs[to_keep])))
        ot = None
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
            vegamags = -2.5 * np.log10(fluxes) - self.zeropoint
            stmags = -2.5 * np.log10(fluxes * self.photflam) - 21.10
            stars_idx = np.where(types == 'point')
            if len(xs[stars_idx]) > 0:
                self._log("info", "Writing {} stars".format(len(xs[stars_idx])))
                self.addPSFs(xs[stars_idx], ys[stars_idx], fluxes[stars_idx], stmags[stars_idx], *args, **kwargs)
                fluxes_observed[stars_idx] = fluxes[stars_idx]
            gals_idx = np.where(types == 'sersic')
            if len(xs[gals_idx]) > 0:
                self._log("info", "Writing {} galaxies".format(len(xs[gals_idx])))
                gxs = xs[gals_idx]
                gys = ys[gals_idx]
                gfluxes = fluxes[gals_idx]
                gns = ns[gals_idx]
                gres = res[gals_idx]
                gphis = phis[gals_idx]
                gratios = ratios[gals_idx]
                gids = ids[gals_idx]
                counter = 1
                total = len(gxs)
                self._log('info', 'Starting Sersic Profiles at {}'.format(time.ctime()))
                for (x, y, flux, n, re, phi, ratio, id) in zip(gxs, gys, gfluxes, gns, gres, gphis, gratios, gids):
                    item_index = np.where(ids == id)[0][0]
                    self._log("info", "Index is {}".format(item_index))
                    if fast_galaxy:
                        central_flux = self.oldSersicProfile(x, y, flux, n, re, phi, ratio, convolve_galaxy, *args, **kwargs)
                    else:
                        central_flux = self.addSersicProfile(x, y, flux, n, re, phi, ratio, *args, **kwargs)
                    fluxes_observed[item_index] = central_flux
                    notes[item_index] = "{}: surface brightness {:.3f} yielded flux {:.3f}".format(notes[item_index], flux, central_flux)
                    self._log("info", "Finished Galaxy {} of {}".format(counter, total))
                    counter += 1
                self._log('info', 'Finishing Sersic Profiles at {}'.format(time.ctime()))
            ot = Table()
            ot['x'] = Column(data=xfs+self.out_origin, unit='pixel')
            ot['y'] = Column(data=yfs+self.out_origin, unit='pixel')
            ot['type'] = Column(data=types)
            ot['vegamag'] = Column(data=vegamags)
            ot['stmag'] = Column(data=stmags)
            ot['countrate'] = Column(data=fluxes_observed, unit=(u.photon/u.second))
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
        self._log("info", "Adding catalogue {} to AstroImage {}".format(catname, self.name))
        obs_file_name = "{}_observed_{}.{}".format(catbase, self.name, self.cat_type)
        obsname = os.path.join(self.out_path, obs_file_name)
        self.addHistory("Adding items from catalogue {}".format(catname))
        in_data = StipsDataTable.dataTableFromFile(cat)
        out_data = StipsDataTable.dataTableFromFile(obsname)
        out_data.meta = {'name': 'Observed Catalogue', 'type': 'observed', 'detector': self.name, 'input_catalogue': catname}
        counter = 0
        current_chunk = in_data.read_chunk()
        while current_chunk is not None:
            table_length = len(current_chunk['id'])
            out_chunk = self.addTable(current_chunk, dist, fast_galaxy = self.fast_galaxy, *args, **kwargs)
            if out_chunk is not None:
                out_data.write_chunk(out_chunk)
            counter += table_length
            current_chunk = in_data.read_chunk()
        self._log("info", "Added catalogue {} to AstroImage {}".format(catname, self.name))
        return obsname

    def make_epsf_array(self, psf_type='normal'):
        """
        Import the 9 PSFs per detector, one in each corner and middle
        of the detector. These 9 PSFs will be used to interpolate the
        best PSF at the location of a source. An ePSF will be calculated
        for each input PSF.

        Parameters
        ----------
        detector : str
            Name of WFI detector, e.g. 'SCA01'
        band : str
            Name of filter, e.g. 'F087'
        psf_type : str
            either 'normal', 'bright', or 'xbright' for
            dim, bright, and extra bright sources.

        Returns
        -------
        psf_array : list
            3 x 3 list with the corresponding input ePSFs.
        """

        have_psf = False

        # Assign prefix based on PSF size
        prefix = ''
        if psf_type != 'normal':
            prefix = psf_type+'_'

        psf_name = "{}psf_{}_{}_{}_{}.fits".format(prefix,
                                                   self.instrument,
                                                   stips_version,
                                                   self.filter,
                                                   self.detector.lower())

        psf_cache_dir = SelectParameter('psf_cache_location')
        psf_cache_name = SelectParameter('psf_cache_directory')
        if psf_cache_name not in psf_cache_dir:
            psf_cache_dir = os.path.join(psf_cache_dir, psf_cache_name)
        if SelectParameter('psf_cache_enable'):
            if not os.path.exists(psf_cache_dir):
                os.makedirs(psf_cache_dir)
            psf_file = os.path.join(psf_cache_dir, psf_name)
            self._log("info", "PSF File {} to be put at {}".format(psf_name, psf_cache_dir))
            self._log("info", "PSF File is {}".format(psf_file))
            if os.path.exists(psf_file):
                try:
                    from webbpsf.utils import to_griddedpsfmodel
                    self.psf = to_griddedpsfmodel(psf_file)
                    have_psf = True
                except Exception as e:
                    self._log("error", "Creating psf from file {}  failed with {}".format(psf_file, e))

        if not have_psf:
            ins = self.psf_constructor
            ins.filter = self.filter
            ins.detector = self.detector
            # Supersample the pixel scale to get WebbPSF to output
            # PSF models with even supersampling centered at the center of a pixel
            ins.pixelscale = self.scale[0] / PSF_UPSCALE
            # Figure out PSF size:
            #   size depends on type and is multiplied by the upscale because we made the
            #   pixels smaller above.
            psf_sizes = {'normal': PSF_BOXSIZE,
                         'bright': PSF_BRIGHT_BOXSIZE,
                         'xbright': PSF_EXTRA_BRIGHT_BOXSIZE}
            fov_pix = psf_sizes[psf_type] * PSF_UPSCALE
            if fov_pix % 2 == 0:
                fov_pix += 1

            num_psfs = self.psf_grid_size*self.psf_grid_size
            if SelectParameter('psf_cache_enable'):
                save = True
                overwrite = True
                psf_file = "{}psf_{}_{}_{}".format(prefix,
                                                   self.instrument,
                                                   stips_version,
                                                   self.filter)
            else:
                save = False
                overwrite = False
                psf_file = None
            msg = "{0}: Starting {1}x{1} PSF Grid creation at {2}"
            self._log("info", msg.format(self.name, self.psf_grid_size, time.ctime()))
            self.psf = ins.psf_grid(all_detectors=False, num_psfs=num_psfs,
                                    fov_pixels=fov_pix, normalize='last',
                                    oversample=1, save=save,
                                    outdir=psf_cache_dir, outfile=psf_file,
                                    overwrite=overwrite)
            msg = "{}: Finished PSF Grid creation at {}"
            self._log("info", msg.format(self.name, time.ctime()))

        psf_data = self.psf.data

        # Create ePSF
        epsf_0_0 = make_epsf(psf_data[0])
        epsf_0_1 = make_epsf(psf_data[1])
        epsf_0_2 = make_epsf(psf_data[2])
        epsf_1_0 = make_epsf(psf_data[3])
        epsf_1_1 = make_epsf(psf_data[4])
        epsf_1_2 = make_epsf(psf_data[5])
        epsf_2_0 = make_epsf(psf_data[6])
        epsf_2_1 = make_epsf(psf_data[7])
        epsf_2_2 = make_epsf(psf_data[8])

        # Reshape Array of ePSFs
        psf_array = [[epsf_0_0, epsf_0_1, epsf_0_2],
                     [epsf_1_0, epsf_1_1, epsf_1_2],
                     [epsf_2_0, epsf_2_1, epsf_2_2]]

        psf_middle = rind((psf_data[0].shape[0]-1) / 2)

        return psf_array, psf_middle

    def addPSFs(self, xs, ys, fluxes, mags, *args, **kwargs):
        """Adds a set of point sources to the image given their co-ordinates and count rates."""
        self.addHistory("Adding {} point sources".format(len(xs)))
        self._log("info", "Adding {} point sources to AstroImage {}".format(len(xs), self.name))

        # Add each source to the image using the ePSF routines
        image_size = self.data.shape[0]

        # Read input PSF files
        psf_array, psf_middle = self.make_epsf_array()
        boxsize = np.floor(psf_middle)/PSF_UPSCALE

        # Are there bright stars?
        are_bright = np.any(mags < self.bright_limit)
        are_xbright = np.any(mags < self.xbright_limit)
        # If so, generate extra ePSF arrays
        if are_bright:
            bright_psf_array, bright_psf_middle = self.make_epsf_array('bright')
            bright_boxsize = np.floor(bright_psf_middle)/PSF_UPSCALE
        if are_xbright:
            xbright_psf_array, xbright_psf_middle = self.make_epsf_array('xbright')
            xbright_boxsize = np.floor(xbright_psf_middle)/PSF_UPSCALE

        for k, (xpix, ypix, flux, mag) in enumerate(zip(xs, ys, fluxes, mags)):
            self.addHistory("Adding point source {} at {},{}".format(k+1, xpix, ypix))
            self._log("info", "Adding point source {} to AstroImage {},{}".format(k+1, xpix, ypix))

            if not self.has_psf:  # just add point source
                self._log("warning", "No PSF found, adding as point source")
                self.data[ypix, xpix] += flux
                continue  # break the loop

            # Create interpolated ePSF from input PSF files
            if mag > self.bright_limit:
                epsf = interpolate_epsf(xpix, ypix, psf_array, image_size)
                self.data = place_source(xpix, ypix, flux, self.data, epsf, boxsize=boxsize, psf_center=psf_middle)
            elif (self.xbright_limit < mag < self.bright_limit):
                self.addHistory("Placing Bright Source with mag = {}".format(mag))
                self._log("info", "Placing Bright Source with mag = {}".format(mag))
                epsf = interpolate_epsf(xpix, ypix, bright_psf_array, image_size)
                self.data = place_source(xpix, ypix, flux, self.data, epsf, boxsize=bright_boxsize, psf_center=bright_psf_middle)
            elif mag < self.xbright_limit:
                self.addHistory("Placing Extra Bright Source with mag = {}".format(mag))
                self._log("info", "Placing Extra Bright Source with mag = {}".format(mag))
                epsf = interpolate_epsf(xpix, ypix, xbright_psf_array, image_size)
                self.data = place_source(xpix, ypix, flux, self.data, epsf, boxsize=xbright_boxsize, psf_center=xbright_psf_middle)

    def cropToBaseSize(self):
        """
        Crops out the PSF overlap on both sides of the detector
        """
        if not self.has_psf:  # already done
            return
        self_y, self_x = self.shape
        base_y, base_x = self.base_shape

        current_centre = (rind(self.shape[0]/2), rind(self.shape[1]/2))
        base_centre = (rind(self.base_shape[0]/2), rind(self.base_shape[1]/2))

        self._log('info', "Cropping convolved image down to detector size")
        ly, hy = current_centre[0]-base_centre[0], current_centre[0]+base_centre[0]
        lx, hx = current_centre[1]-base_centre[1], current_centre[1]+base_centre[1]
        if hx-lx < self.base_shape[1]:
            hx += 1
        elif hx-lx > self.base_shape[1]:
            hx -= 1
        if hy-ly < self.base_shape[0]:
            hy += 1
        elif hy-ly > self.base_shape[0]:
            hy -= 1
        self._log('info', "Taking [{}:{}, {}:{}]".format(ly, hy, lx, hx))
        fp_crop = np.zeros(tuple(self.base_shape), dtype='float32')
        fp_crop[:, :] = self.data[ly:hy, lx:hx]
        crpix = [base_centre[0], base_centre[1]]
        if self.wcs.sip is not None:
            sip = wcs.Sip(self.wcs.sip.a, self.wcs.sip.b, None, None, crpix)
        else:
            sip = None
        self.wcs = self._wcs(self.ra, self.dec, self.pa, self.scale, crpix=crpix, sip=sip)
        self.data = fp_crop
        self.shape = self.base_shape
        self.has_psf = False

    def addSersicProfile(self, posX, posY, flux, n, re, phi, axialRatio,
                         *args, **kwargs):
        """
        Create a simulated Sersic profile, including PSF convolution, using
        Pandeia's SersicDistribution.

        Parameters
        ----------
        posX: float
            The X position of the profile center on the detector.
        posX: float
            The Y position of the profile center on the detector.
        flux: float
            The surface brightness (in counts/s/pixel) at the half-light radius.
        n: float
            The Sersic index. This must be between 0.3 and 6.2.
        re: float
            The half-light radius, in pixels.
        phi: float
            The position angle of the major axis, in degrees east of north.
        axialRatio: float
            The ratio between the major and minor axes.

        Returns
        -------
        central_flux: float
            The flux at the central pixel of the model
        """
        ix, iy = rind(posX), rind(posY)
        if flux == 0.:
            return 0.

        flux = sersic_lum(flux, re, n)

        source_dict = {'id': id}
        source_dict['shape'] = {'geometry': 'sersic'}
        source_dict['shape']['major'] = re
        source_dict['shape']['minor'] = re*axialRatio
        source_dict['shape']['sersic_index'] = n
        source_dict['spectrum'] = {'redshift': 0., 'sed': {'sed_type': 'flat'}}
        source_dict['normalization'] = {'type': 'none'}
        source_dict['position'] = {'orientation': 90-phi}
        source_dict['position']['x_offset'] = posX - self.xsize / 2
        source_dict['position']['y_offset'] = self.ysize / 2 - posY
        g = coords.Grid(1., 1., self.xsize, self.ysize)
        xc, yc = ix, iy

        # Roman is the only telescope supported
        src = source.Source('roman', config=source_dict)
        src.grid = g
        sersic = profile.SersicDistribution(src, 1)
        img = deepcopy(sersic.prof)*flux/np.sum(sersic.prof)
        central_flux = img[yc, xc]

        # Read input PSF files
        psf_array, psf_middle = self.make_epsf_array()
        boxsize = rind(np.floor(psf_middle)/PSF_UPSCALE)
        epsf = interpolate_epsf(posX, posY, psf_array, self.shape[0])
        psf_data = np.zeros((boxsize, boxsize), dtype=np.float32)
        psf_x, psf_y = rind(boxsize/2), rind(boxsize/2)
        psf_data = place_source(psf_x, psf_y, 1., psf_data, epsf, boxsize=boxsize, psf_center=psf_middle)
        result = convolve_fft(img, psf_data)
        self.data += result
        return central_flux


    def oldSersicProfile(self, posX, posY, flux, n, re, phi, axialRatio,
                         convolve_galaxy, *args, **kwargs):
        """
        Create a simulated Sersic profile, including PSF convolution, using
        the old (and faster) Astropy Sersic2D model.

        Parameters
        ----------
        posX: float
            The X position of the profile center on the detector.
        posX: float
            The Y position of the profile center on the detector.
        flux: float
            The surface brightness (in counts/s/pixel) at the half-light radius.
        n: float
            The Sersic index. This must be between 0.3 and 6.2.
        re: float
            The half-light radius, in pixels.
        phi: float
            The position angle of the major axis, in degrees east of north.
        axialRatio: float
            The ratio between the major and minor axes.
        convolve_galaxy : bool
            Whether to convolve the final image with the Roman PSF.
            [default: True]

        Returns
        -------
        central_flux: float
            The flux at the central pixel of the model
        """
        if flux == 0.:
            return 0.

        # Generate 2D Sersic profile
        from astropy.modeling.models import Sersic2D
        x, y = np.meshgrid(np.arange(self.xsize), np.arange(self.ysize))
        mod = Sersic2D(amplitude=flux, r_eff=re, n=n, x_0=posX, y_0=posY, ellip=(1.-axialRatio), theta=(np.radians(phi) + 0.5*np.pi))
        img = mod(x, y)

        # Convolve with Roman PSF
        if convolve_galaxy:
            psf_array, psf_middle = self.make_epsf_array()
            boxsize = rind(np.floor(psf_middle)/PSF_UPSCALE)
            epsf = interpolate_epsf(posX, posY, psf_array, self.shape[0])
            psf_data = np.zeros((boxsize, boxsize), dtype=np.float32)
            psf_x, psf_y = rind(boxsize/2), rind(boxsize/2)
            psf_data = place_source(psf_x, psf_y, 1., psf_data, epsf, boxsize=boxsize, psf_center=psf_middle)
            result = convolve_fft(img, psf_data)
            self.data += result
        else:
            self.data += img

        # Calculate central flux
        central_flux = mod(posX, posY)
        return central_flux

    def rotate(self, angle, reshape=False):
        """
        Rotate the image a number of degrees as specified

        ..warning:: This function is not necessarily flux-conserving
        """
        self.addHistory("Rotating by {} degrees".format(angle))
        self._log("info", "Rotating AstroImage {} by {} degrees".format(self.name, angle))
        self.pa = (self.pa + angle) % 360. % 360.
        self.data = rotate(self.data, angle, order=5, reshape=reshape)

    def addWithOffset(self, other, offset_x, offset_y):
        """
        Adds other to self, but with an offset.

        offset_n have units of pixels-of-self
        """
        self.addHistory("Adding image with offset of ({},{}) pixels".format(offset_x, offset_y))
        self._log("info", "Adding image {} with offset ({},{})".format(other.name, offset_x, offset_y))
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

        low_x = x_overlay - xc  # 0 on other array
        low_y = y_overlay - yc  # 0 on other array
        high_x = low_x + xs  # xs on other array
        high_y = low_y + ys  # ys on other array

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

        self._log("info", "Adding Other[{}:{},{}:{}] to AstroImage {}[{}:{},{}:{}]".format(ly, hy, lx, hx, self.name, low_y, high_y, low_x, high_x))
        return int(lx), int(ly), int(hx), int(hy), int(low_x), int(low_y), int(high_x), int(high_y)

    def _filled(self, offset_x, offset_y, size_x, size_y):
        """
        If an image with size (size_x, size_y) is added to the current image with the pixel offset
            (offset_x, offset_y), will that image fill the frame of the current image?
        """
        xs, ys = size_x, size_y
        lx, ly, hx, hy, low_x, low_y, high_x, high_y = self._findCentre(xs, ys, offset_x, offset_y)

        if low_y == 0 and low_x == 0 and high_y == self.shape[0] and high_x == self.shape[1]:
            return True
        return False

    def _addArrayWithOffset(self, other, offset_x, offset_y):
        """
        Add an AstroImage to the current image with a given pixel offset. This allows for a smaller
        memory use for addWithAlignment.
        """
        ys, xs = other.shape
        lx, ly, hx, hy, low_x, low_y, high_x, high_y = self._findCentre(xs, ys, offset_x, offset_y)

        # If any of these are false, the images are disjoint.
        if low_x < self.xsize and low_y < self.ysize and high_x > 0 and high_y > 0:
            self.data[low_y:high_y, low_x:high_x] += other[ly:hy, lx:hx]
        else:
            self.addHistory("Added image is disjoint")
            self._log("warning", "{}: Image is disjoint".format(self.name))

    def _addWithOffset(self, other, offset_x, offset_y):
        """
        Add an AstroImage to the current image with a given pixel offset. This allows for a smaller
        memory use for addWithAlignment.
        """
        ys, xs = other.shape
        lx, ly, hx, hy, low_x, low_y, high_x, high_y = self._findCentre(xs, ys, offset_x, offset_y)

        # If any of these are false, the images are disjoint.
        if low_x < self.xsize and low_y < self.ysize and high_x > 0 and high_y > 0:
            self.data[low_y:high_y, low_x:high_x] += other.data[ly:hy, lx:hx]
        else:
            self.addHistory("Added image is disjoint")
            self._log("warning", "{}: Image is disjoint".format(self.name))

    def addWithAlignment(self, other):
        """Adds other to self after aligning co-ordinates"""
        self.addHistory("Adding image aligned by WCS")
        self._log("info", "Adding image {} (RA,DEC,PA)=({},{},{}) aligned by WCS".format(other.name, other.ra, other.dec, other.pa))
        other_copy = other.copy()
        if abs(other.pa - self.pa) > 1.e-5:
            self._log("info", "Rotating other image by {} degrees".format((self.pa-other.pa)))
            other_copy.rotate((self.pa-other.pa), reshape=True)
        if abs(other.scale[0] - self.scale[0]) > 1.e-5 or abs(other.scale[1]-self.scale[1]) > 1.e-5:
            self._log("info", "Rescaling other from ({},{}) to ({},{})".format(other.scale[0], other.scale[1], self.scale[0], self.scale[1]))
            other_copy.rescale(self.scale)
            self._log("info", "Finished rescaling")
        ra = other.ra
        dec = other.dec
        pix = self.wcs.wcs_world2pix([ra], [dec], self.wcs_origin, ra_dec_order=True)
        px = pix[0][0]
        py = pix[1][0]
        offset_x = int(np.round(px - self.wcs.wcs.crpix[0]))
        offset_y = int(np.round(py - self.wcs.wcs.crpix[1]))
        self._log("info", "Other has centre co-ordinates ({},{}) for offset ({},{})".format(px, py, offset_x, offset_y))
        self._addWithOffset(other_copy, offset_x, offset_y)

    def rescale(self, scale):
        """Rescales the image to the provided plate scale."""
        self.addHistory("Rescaling to ({},{}) arcsec/pixel".format(scale[0], scale[1]))
        self._log("info", "Rescaling to ({},{}) arcsec/pixel".format(scale[0], scale[1]))
        shape_x = rind(self.shape[1] * self.scale[0] / scale[0])
        shape_y = rind(self.shape[0] * self.scale[1] / scale[1])
        new_shape = (shape_y, shape_x)
        self._log("info", "New shape will be {}".format(new_shape))
        fp_result = np.zeros(new_shape, dtype='float32')
        flux = self.data.sum()
        self._log("info", "Max flux is {}, sum is {}".format(np.max(self.data), flux))
        zoom(self.data, (np.array(self.scale)/np.array(scale)), fp_result)
        factor = flux / fp_result.sum()
        fp_result *= factor
        self._log("info", "Max flux is {}, sum is {}".format(np.max(fp_result), np.sum(fp_result)))
        self.data = fp_result
        self.shape = new_shape
        self._scale = scale
        self.wcs = self._wcs(self.ra, self.dec, self.pa, scale)
        self._prepHeader()

    def crop(self, xs, ys, offset_x, offset_y):
        """Crops the image. If the offset < 0, or offset+size > current size, pads with blank pixels"""
        self._addWithOffset(self.copy(), offset_x, offset_y)
        offset_value = np.array([offset_x, offset_y])
        pix_coords = np.array(offset_value+self.wcs.wcs.crpix)
        if self.distorted:
            world_coords = self.wcs.all_pix2world(pix_coords[0], pix_coords[1],
                                                  self.wcs_origin, ra_dec_order=True)
        else:
            world_coords = self.wcs.wcs_pix2world(pix_coords[0], pix_coords[1],
                                                  self.wcs_origin, ra_dec_order=True)
        self.wcs = self._wcs(world_coords[0], world_coords[1], self.pa, self.scale, self.wcs.sip)
        self._prepHeader()

    def setExptime(self, exptime):
        """
        Set the exposure time. Because the data is kept in units of counts/s, this does
        not affect the data values, but it does affect errors residuals (e.g. dark
        current and cosmic ray count)
        """
        factor = exptime / self.exptime
        self.data *= factor
        self.exptime = exptime
        self.updateHeader('exptime', self.exptime)

        #self.exptime = exptime
        #self.updateHeader('exptime', self.exptime)

    def setUnits(self):
        """
        Divide the image by the exposure time to the output is in e/s as opposed to e.
        """
        self.data /= self.exptime

    def addBackground(self, background):
        """
        Add a constant value per-pixel. Automatically scale that value based on oversampling.
        """
        self.addHistory("Added background of {} counts/s/pixel".format(background))
        self._log("info", "Added background of {} counts/s/pixel".format(background))
        self.data += background

    def introducePoissonNoise(self, absVal=False):
        """
        Generate Poisson noise and add it to the internal image.

        Based on `CALC_SHOT_NOISE` in `instrument__init.pro`
        from JWST IDL Simulator package.

        Parameters
        ----------
        absVal: bool, optional
            Take absolute value of `noiseData`.
        """
        abs_data = np.absolute(self.data)
        noise_data = np.random.RandomState(seed=self.seed).normal(size=self.shape) * np.sqrt(abs_data)
        if absVal:
            noise_data = np.abs(noise_data)
        mean, std = noise_data.mean(), noise_data.std()
        self.data += noise_data
        self.addHistory("Adding Poisson Noise with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding Poisson Noise with mean {} and standard deviation {}".format(mean, std))

    def introduceReadnoise(self, readnoise):
        """
        Simulate kTC or read noise for detector.

        Parameters
        ----------
        readnoise: constant representing average read noise per pixel.
        """
        noise_data = np.zeros(tuple(self.shape), dtype='float32')
        noise_data[:, :] = readnoise * np.random.RandomState(seed=self.seed).randn(self.ysize, self.xsize)
        mean, std = noise_data.mean(), noise_data.std()
        self.data += noise_data
        self.updateHeader("RDNOISE", readnoise)
        self.addHistory("Adding Read noise with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding readnoise with mean {} and STDEV {}".format(mean, std))

    def introduceFlatfieldResidual(self, flat):
        """
        Simulate flatfield correction error.

        flat: AstroImage containing flatfield error values.

        returns: mean, std. Mean and standard deviation of used portion of error array.
        """
        err = flat.data[:self.ysize, :self.xsize]
        mean, std = err.mean(), err.std()
        self.data *= err
        self.addHistory("Adding Flatfield residual with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding Flatfield residual with mean {} and standard deviation {}".format(mean, std))

    def introduceDarkResidual(self, dark):
        """
        Simulate dark correction error.

        dark: AstroImage containing dark residual.

        returns: mean,std: mean and standard deviation of dark error array.
        """
        err = dark.data[:self.ysize, :self.xsize]
        mean, std = err.mean(), err.std()
        self.data += err
        self.addHistory("Adding Dark residual with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding Dark residual with mean {} and standard deviation {}".format(mean, std))

    def introduceCosmicRayResidual(self, pixel_size):
        """
        Simulate CR correction error.

        pixel_size: pixel size in microns (on detector)

        returns: mean,std: float. Mean and standard deviation of cosmic ray image that was added.
        """
        energies = (600.0, 5000.0)  # e- (gal, SS)
        rates = (5.0, 5.0)  # hits/cm^2/s (gal, SS)
        pixarea = pixel_size**2 / 1E8  # cm^2
        probs = GetCrProbs(rates, pixarea, self.exptime)  # hits
        cr_size, cr_psf = GetCrTemplate()

        noise_data = np.zeros(self.shape, dtype='float32')
        for i in range(len(energies)):
            noise_data += MakeCosmicRay(self.shape[1], self.shape[0], probs[i], energies[i], cr_size, cr_psf, self.seed, verbose=False)
        noise_data *= 0.01
        mean, std = noise_data.mean(), noise_data.std()
        self.data += noise_data

        self.addHistory("Adding Cosmic Ray residual with mean {} and standard deviation {}".format(mean, std))
        self._log("info", "Adding Cosmic Ray residual with mean {} and standard deviation {}".format(mean, std))

    def _prepHeader(self):
        """
        Prepare the header WCS based on instrument and filter. For now, set RA and DEC to be
        identically zero.
        """
        self.header['DATE-OBS'] = time.strftime("%Y-%M-%d")
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

    def _getSip(self, **kwargs):
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
                sip = self._sip(dist_a, dist_b, dist_ap, dist_bp)
                self._log('info', "SIP A Array: {}".format(sip.a))
                self._log('info', "SIP B Array: {}".format(sip.b))
        return sip

    def _getWcs(self, **kwargs):
        """Makes a WCS from the input arguments"""
        sip = self._getSip(**kwargs)
        if 'wcs' in kwargs:
            wcs = kwargs['wcs']
            if sip is not None and wcs.sip is None:
                wcs.sip = sip
            wcs = self._normalizeWCS(wcs)
        else:
            # get pixel scale (if available)
            scale = self._scale
            ra = kwargs.get('ra', 0.)
            offset_ra = kwargs.get('offset_ra', 0.)
            dec = kwargs.get('dec', 0.)
            offset_dec = kwargs.get('offset_dec', 0.)
            ra, dec = OffsetPosition(ra, dec, offset_ra, offset_dec)
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
        w.wcs.ctype = ["", ""]
        w.wcs.ctype[ranum] = "RA---TAN"
        w.wcs.ctype[decnum] = "DEC--TAN"
        if crpix is None:
            w.wcs.crpix = [rind(self.xsize/2), rind(self.ysize/2)]
        else:
            w.wcs.crpix = crpix
        w.wcs.crval = [0., 0.]
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

    def _normalizeWCS(self, w):
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
            raise ValueError("Lattype = {} and Lngtype = {}: can't get RA, DEC".format(w.wcs.lngtyp, w.wcs.lattyp))
        ra = w.wcs.crval[ranum]
        dec = w.wcs.crval[decnum]
        if w.wcs.has_cd():
            scale = wcs.utils.proj_plane_pixel_scales(w)*3600.
        else:
            scale = w.wcs.cdelt*3600.
        pa = self._getPA(w, scale, decnum)
        return self._wcs(ra, dec, pa, scale, w.wcs.crpix, w.sip, ranum, decnum)

    def _getPA(self, wcs, scale, decnum=1):
        """With non-normalized WCS, find the angle of the pixel offset of increasing DEC by 1 arcsec"""
        offset_value = np.array([0., 0.])
        offset_value[decnum] += 1./3600.
        world_coords = np.array((wcs.wcs.crval, wcs.wcs.crval+offset_value))
        pix_coords = wcs.wcs_world2pix(world_coords, self.wcs_origin)
        offset = (pix_coords[1] - pix_coords[0])
#         offset = offset / scale # Divide by scale to correct for non-square pixels
        pa_north = np.degrees(np.arctan2(offset[0], offset[1])) % 360. % 360.
        return pa_north

    def __add__(self, other):
        """Adds one AstroImage to another."""
        result = self.copy()
        result.addWithAlignment(other)
        for k, v in other.header.iteritems():
            if k not in result.header:
                result.header[k] = v
        for item in other.history:
            result.history.append("{}: {}".format(other.name, item))
        result.history.append("Took {} and added {}".format(self.name, other.name))
        return result

    def __iadd__(self, other):
        """Adds an AstroImage to the current one. (i.e. the '+=' operator)"""
        if isinstance(other, int) or isinstance(other, float):
            self.addHistory("Added constant {}/pixel".format(other))
            self.data += other
            return self
        else:
            self.addWithAlignment(other)
            for k, v in other.header.iteritems():
                if k not in self.header:
                    self.header[k] = v
            for item in other.history:
                self.history.append("{}: {}".format(other.name, item))
            self.history.append("Added {}".format(other.name))
            return self

    def __radd__(self, other):
        """
        Adds an integer or floating-point constant to the AstroImage

        .. warning:: Assumes constant value per-pixel
        """
        self.addHistory("Added {}/pixel".format(other))
        self.data += other
        return self

    def __mul__(self, other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        result = self.copy()
        self.data *= other
        result.addHistory("Multiplied by {}".format(other))
        return result

    def __imul__(self, other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        self.addHistory("Multiplied by {}".format(other))
        self.data *= other
        return self

    def __rmul__(self, other):
        """Multiples an integer or floating-point constant to the AstroImage"""
        result = self.copy()
        self.data *= other
        result.addHistory("Multiplied by {}".format(other))
        return result

    @property
    def sum(self):
        """Returns the sum of the flux in the current image"""
        return np.sum(self.data)

    def _log(self, mtype, message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self, 'logger'):
            getattr(self.logger, mtype)(message)
        else:
            sys.stderr.write("{}: {}\n".format(mtype, message))

    def _init_dat(self, base_shape, psf_shape=(0, 0), data=None):
        self.base_shape = base_shape
        self.shape = np.array(base_shape) + np.array(psf_shape)
        t_shape = tuple(self.shape)
        self.data = np.zeros(t_shape, dtype='float32')
        if data is not None:
            centre = rind(self.shape/2)
            cy, cx = centre
            half = rind(base_shape/2)
            hy, hx = half
            self.data[cy-hy:cy+base_shape[0]-hy, cx-hx:cx+base_shape[1]-hx] = data

    def _remap(self, xs, ys):
        # Step 1 -- compensate for PSF adjustments
        adj = ((np.array(self.shape) - np.array(self.base_shape))/2).astype(np.int32)
        out_y = ys - adj[0]
        out_x = xs - adj[1]
        return out_x, out_y

    INSTRUMENT_DEFAULT = {
                            'telescope': 'roman',
                            'instrument': 'WFI',
                            'filter': 'F062',
                            'detector': {
                                            'WFI': 'SCA01',
                                            'NIRCam': "NRCA1",
                                            'MIRI': 'MIRI'
                                        },
                            'shape': (4088, 4088),
                            'scale': [0.11, 0.11],
                            'zeropoint': 21.0,
                            'photflam': 0.,
                            'photplam': 0.6700,
                         }

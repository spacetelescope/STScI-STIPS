"""
General CGI form functions.

:Author: Pey Lian Lim

:Organization: Space Telescope Science Institute

"""
from __future__ import absolute_import,division

# External modules
import importlib, inspect, os, shutil, socket, struct, sys, urllib, uuid
import numpy as np
import astropy.io.fits as pyfits
from numpy.fft import fft2, ifft2
from astropy.io import ascii
from astropy.table import Table
from jwst_backgrounds.jbt import background
from photutils.psf.models import GriddedPSFModel

from .. import __version__ as __stips__version__

#-----------
class classproperty(object):
    def __init__(self, f):
        self.f = classmethod(f)
    def __get__(self, *a):
        return self.f.__get__(*a)()
#-----------
class __grid__(object):
    @classproperty
    def __pandeia__version__(self):
        with open(GetStipsData(os.path.join("grid", "VERSION.txt")), "r") as inf:
            line = inf.readline()
            items = line.strip().split()
            return items[1]
    
    @classproperty
    def __stips__version__(self):
        with open(GetStipsData(os.path.join("grid", "VERSION.txt")), "r") as inf:
            line = inf.readline()
            line = inf.readline()
            items = line.strip().split()
            return items[1]

#-----------
class ImageData(object):
    def __init__(self, fname, shape, mode='r+', memmap=True):
        self.shape = tuple(shape)
        if memmap:
            self.fp = np.memmap(fname, dtype='float32', mode=mode, 
                                shape=self.shape)
        else:
            if isinstance(fname, np.ndarray):
                self.fp = fname
            else:
                self.fp = np.ndarray(self.shape, dtype='float32')
    
    def __enter__(self):
        return self.fp
    
    def __exit__(self, exc_type, exc_value, traceback):
        del self.fp
    
    @property
    def data(self):
        return self.fp

#-----------
class CachedJbtBackground(background, object):
    '''
    Version of the JBT Background class intended to work with a local cached file. Repeats the 
    entire __init__ function because there doesn't seem to be any way to give jbt a local file
    path instead of a URL.
    '''
    def __init__(self, ra, dec, wavelength, thresh=1.1):
        # global attributes
        self.cache_path = GetStipsData("background")
        if sys.version_info[0] >= 3:
            self.cache_url = "file://" + urllib.request.pathname2url(os.path.join(self.cache_path, "remote_cache/"))
        else:
            self.cache_url = urllib.pathname2url(os.path.join(self.cache_path, "remote_cache/"))
        self.local_path = os.path.join(self.cache_path, 'jbt_refdata')
        self.wave_file = 'std_spectrum_wavelengths.txt' # The wavelength grid of the background cache
        self.thermal_file = 'thermal_curve_jwst_jrigby_cchen_1.1a.csv' # The constant (not time variable) thermal self-emission curve
        self.nside = 128  # Healpy parameter, from generate_backgroundmodel_cache.c .  
        self.wave_array,self.thermal_bg = self.read_static_data()
        self.sl_nwave = self.wave_array.size  # Size of wavelength array
        
        # input parameters
        self.ra = ra
        self.dec = dec
        self.wavelength = wavelength
        self.thresh = thresh

        # Load variable content
        self.cache_file = self.myfile_from_healpix(ra, dec)
        self.bkg_data = self.read_bkg_data(self.cache_file)
                
        # Interpolate bathtub curve and package it    
        self.make_bathtub(wavelength)

    def read_bkg_data(self, cache_file, verbose=False):
        """
        Method for reading one JWST background file, and parsing it.    
        
        Schema of each binary file in the cache:
        ----------------------------------------
        JRR verified the schema against the source code, generate_stray_light_with_threads.c. 
        The cache uses a Healpix RING tesselation, with NSIDE=128.  Every point on the sky 
        (tesselated tile) corresponds to one binary file, whose name includes its healpix 
        pixel number, in a directory corresponding to the first 4 digits of the healpix number.  

        - double RA
        - double DEC
        - double pos[3]
        - double nonzodi_bg[SL_NWAVE]
        - int[366] date_map  : This maps dates to indices.  NOTE: There are 366 days, not 365!
        - for each day in FOR:
            - double zodi_bg[SL_NWAVE]
            - double stray_light_bg[SL_NWAVE]
        
        parameters
        ----------
        cache_file: string
        
        attributes
        ----------
        
        """

        # Read the background file via http 
        if sys.version_info[0] >= 3:
            # Python 3
            # Read the background cache version
            version_file = urllib.request.urlopen(self.cache_url + 'VERSION')
            sbet_file = urllib.request.urlopen(self.cache_url + cache_file)
        else:
            # Python 2
            # Read the background cache version
            version_file = urllib.urlopen(self.cache_url + 'VERSION')
            sbet_file = urllib.urlopen(self.cache_url + cache_file)

        self.cache_version = version_file.readlines()[0].decode('utf-8')[:-1]
        sbet_data = sbet_file.read()
        
        # Unpack the constant first part
        if verbose: 
            print("File has", len(sbet_data), "bytes, which is", len(sbet_data)/8., "doubles")
            
        size_calendar = struct.calcsize("366i") # bytes, not doubles
        partA = struct.unpack(str(5 + self.sl_nwave)+'d', sbet_data[0: (5 + self.sl_nwave)*8])
        ra = partA[0]
        dec = partA[1]
        pos = partA[2:5]
        nonzodi_bg = np.array(partA[5:5+self.sl_nwave])

        # Unpack the calendar dates - the dates go from 0 to 365 days.
        date_map = np.array(struct.unpack('366i', sbet_data[(5 + self.sl_nwave)*8  : (5 + self.sl_nwave)*8 + size_calendar]))
        if verbose: 
            print("Out of", len(date_map), "days, these many are legal:", np.sum(date_map >=0))

        calendar = np.where(date_map >=0)[0]

        Ndays = len(calendar) 
        if verbose: 
            print(len(date_map), Ndays)

        # Unpack part B, the time-variable part
        zodi_bg        = np.zeros((Ndays,self.sl_nwave))
        stray_light_bg = np.zeros((Ndays,self.sl_nwave))
        perday = self.sl_nwave*2
        partB= struct.unpack(str((len(calendar))*self.sl_nwave*2)+'d', sbet_data[perday*Ndays*-8 : ])

        # The index dd in zodi_bg[dd, : ] corresponds to the calendar day lookup[dd]
        for dd in range(0, int(Ndays)):
            br1 = dd*perday
            br2 = br1 + self.sl_nwave
            br3 = br2 + self.sl_nwave
            zodi_bg[dd, ] = partB[br1 : br2]
            stray_light_bg[dd, ] = partB[br2 : br3]

        # Expand static background components to the same shape as zodi_bg
        total_bg = np.tile(nonzodi_bg + self.thermal_bg,(Ndays,1)) + stray_light_bg + zodi_bg

        # pack everything up as a dict
        return {'calendar':calendar, 'ra':ra, 'dec':dec, 'pos':pos, 'wave_array':self.wave_array, 'nonzodi_bg':nonzodi_bg, 
                'thermal_bg':self.thermal_bg, 'zodi_bg':zodi_bg, 'stray_light_bg':stray_light_bg, 'total_bg':total_bg} 

#-----------
def GetStipsData(to_retrieve):
    """
    Retrieve a file from the stips_data directory. Will also print out a warning if the directory
    can't be found.
    """
    stips_version = __stips__version__.strip().replace(".", "")
    if "dev" in stips_version:
        stips_version = stips_version[:stips_version.find("dev")]
    if "stips_data" not in os.environ:
        msg = "ERROR: stips_data environment variable not found. STIPS "
        msg += "requires the data directory to function. Please download the "
        msg += "STIPS data from <https://stsci.box.com/v/stips-data-{}> "
        msg += "and set the stips_data environment variable to point to it.\n"
        sys.stderr.write(msg.format(stips_version))
        raise EnvironmentError("stips_data environment variable not found.")
    stips_data_base = os.environ["stips_data"]
    if not os.path.exists(stips_data_base):
        msg = "ERROR: stips_data directory at {} not ".format(stips_data_base)
        msg += "found. STIPS requires the data directory to function. "
        msg += "Please make sure that the STIPS data directory exists.\n"
        sys.stderr.write(msg)
        raise FileNotFoundError("${stips_data} does not exist.")
    retrieval_file = os.path.join(stips_data_base, to_retrieve)
    if not os.path.exists(retrieval_file):
        msg = "ERROR: STIPS data file {} not found. ".format(retrieval_file)
        msg += "STIPS requires the data directory to function. Please download "
        msg += "a new copy of the STIPS data from "
        msg += "<https://stsci.box.com/v/stips-data-{}>.\n"
        sys.stderr.write(msg.format(stips_version))
        raise FileNotFoundError("File {} does not exist.".format(retrieval_file))
    return retrieval_file

#-----------
def internet(host="8.8.8.8", port=53, timeout=3):
    """
    Host: 8.8.8.8 (google-public-dns-a.google.com)
    OpenPort: 53/tcp
    Service: domain (DNS/TCP)
    """
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except Exception as ex:
        return False

#-----------
def OffsetPosition(in_ra,in_dec,delta_ra,delta_dec):
    """
    Offset a position given in decimal degrees.

    Parameters
    ----------
    in_ra: float
        Initial RA (decimal degrees).

    in_dec: float
        Initial DEC (demical degrees).

    delta_ra: float
        Offset in RA (decimal degrees).
    
    delta_dec: float
        Offset in DEC (decimal degrees).

    Returns
    -------
    ra: float
        Offset RA.
    
    dec: float
        Offset DEC.
    """
    ra = in_ra
    dec = in_dec + delta_dec
    if dec > 90.:
        dec = 180 - dec
        ra = 180 + ra
    if dec < -90.:
        dec = -180 - dec
        ra = 180 + ra
    ra = ra + delta_ra
    if ra > 360.:
        ra = ra - 360.
    if ra < 0.:
        ra = ra + 360.
    return ra,dec

#-----------
def InstrumentList(excludes=[]):
    """
    Looks through the available instrument files, and imports all instruments that do not match
    one of the names in 'excludes', or whose telescope does not match one of the names in 'excludes'
    
    Parameters
    ----------
    
    excludes: list of strings
        Class and Telescope names to exclude when building the array.

    Returns
    -------
    instruments: dictionary
        dictionary of instrument classes (by instrument name)
    """
    instruments = {}
    from ..instruments import __all__ as all_instruments
    for f in all_instruments:
        if f == "__init__":
            continue
        module = importlib.import_module(".instruments."+f, "stips")
        for name, obj in inspect.getmembers(module):
            if inspect.isclass(obj) and hasattr(obj, "TELESCOPE") and hasattr(obj, "INSTRUMENT") and hasattr(obj, "DETECTOR"):
                telescope = getattr(obj, "TELESCOPE").lower()
                instrument = getattr(obj, "INSTRUMENT").lower()
                detector = getattr(obj, "DETECTOR").lower()
                add = True
                for exclude in excludes:
                    excl = exclude.lower()
                    if excl == telescope or excl == instrument or excl == detector or excl == name.lower():
                        add = False
                if add:
                    instruments[name] = obj
    return instruments

def read_metadata(filename, n_lines=100000, format='ascii.ipac'):
    lines = []
    for i, line in enumerate(open(filename, 'r')):
        lines.append(line.rstrip(os.linesep))
        if i == n_lines:
            break
    t = Table.read(lines, format=format)
    return t

def read_table(filename, n_chunk=100000, format="ipac"):
    """
    Chunk reader to (hopefully) not take up so much memory when reading very large tables.
    """
    names = None
    lines = []
    for i, line in enumerate(open(filename, 'r')):
        lines.append(line.rstrip(os.linesep))
        if i % n_chunk == n_chunk - 1:
            if i < n_chunk:  # first
                chunk = ascii.read(lines, format=format, guess=False)
                names = chunk.colnames
                yield chunk
            else:
                chunk = ascii.read(lines, format='no_header', names=names, guess=False)
                yield chunk
            lines = []
    if lines:
        if names is not None:
            yield ascii.read(lines, format='no_header', names=names, guess=False)
        else:
            yield ascii.read(lines, format=format, guess=False)


class Percenter(object):
    def __init__(self, total):
        self.total = total
        self.current = 0
    @property
    def percent(self):
        return "{:.2f}%".format(100.*self.current/self.total)
    def incr(self):
        self.current += 1



def computation(arr, Hf, pos, Nfft, y, ys, shape, lock, path, memmap):
    start_y, end_y, start_x, end_x, thisend_y, thisend_x = pos
    if isinstance(Hf, GriddedPSFModel):
        # Generate PSF
        y, x = np.mgrid[start_y:end_y, start_x:end_x]
        y_0, x_0 = (start_y+end_y)//2, (start_x+end_x)//2
        # Make a grid that's the size of the PSF fov_pix around its centre.
        ys2, xs2 = shape[0]//2, shape[1]//2
        y, x = np.mgrid[y_0-ys2:y_0+ys2, x_0-xs2:x_0+xs2]
        psf = Hf.evaluate(x=x, y=y, flux=1, x_0=x_0, y_0=y_0)
        fftHf = fft2(psf, Nfft)
        conv = np.real(ifft2(fftHf * fft2(arr, Nfft)))
    else:
        conv = np.real(ifft2(Hf * fft2(arr, Nfft)))
    lock.acquire()
    with ImageData(y, ys, memmap=memmap) as dat:
        dat[start_y:thisend_y, start_x:thisend_x] += (conv[:(thisend_y-start_y), :(thisend_x-start_x)])
    lock.release()
    return "[{}, {}]".format(start_y, start_x)


def overlapaddparallel(Amat, amat_shape,
                       Hmat, hmat_shape,
                       L=None, Nfft=None, y=None, verbose=False, logger=None, 
                       state_setter=None, base_state="", path=None, cores=None,
                       memmap=True):
    """
    Fast two-dimensional linear convolution via the overlap-add method.
    The overlap-add method is well-suited to convolving a very large array,
    `Amat`, with a much smaller filter array, `Hmat` by breaking the large
    convolution into many smaller `L`-sized sub-convolutions, and evaluating
    these using the FFT. The computational savings over the straightforward
    two-dimensional convolution via, say, scipy.signal.convolve2d, can be
    substantial for large Amat and/or Hmat.
    Parameters
    ----------
    Amat, Hmat : array_like
        Two-dimensional input arrays to be convolved. For computational
        purposes, Amat should be larger.
    L : sequence of two ints, optional
        Length of sub-convolution to use. This should ideally be as large as
        possible to obtain maximum speedup: that is, if you have the memory to
        compute the linear convolution using FFT2, i.e., replacing spatial
        convolution with frequency-domain multiplication, then let `L =
        np.array(Amat.shape) + np.array(Hmat.shape) - 1`. Usually, though, you
        are considering overlap-add because you can't afford a batch transform
        on your arrays, so make `L` as big as you can.
    Nfft : sequence of two ints, optional
        Size of the two-dimensional FFT to use for each `L`-sized
        sub-convolution. If omitted, defaults to the next power-of-two, that
        is, the next power-of-two on or after `L + np.array(Hmat.shape) - 1`.
        If you choose this default, try to avoid `L` such that this minimum, `L
        + np.array(Hmat.shape) - 1`, is just a bit greater than a power-of-two,
        since you'll be paying for a 4x larger FFT than you need.
    y : array_like, optional
        Storage for the output. Useful if using a memory-mapped file, e.g.
    verbose : boolean, optional
        If True, prints a message for each `L`-sized subconvolution.
    Returns
    -------
    y : same as passed in, or ndarray if no `y` passed in
        The `np.array(Amat.shape) + np.array(Hmat.shape) - 1`-sized
        two-dimensional array containing the linear convolution. Should be
        within machine precision of, e.g., `scipy.signal.convolve2d(Amat,
        Hmat, 'full')`.
    Raises
    ------
    ValueError if `L` and `Nfft` aren't two-element, and too small: both
    elements of `L` must be greater than zero, and `Nfft`'s must be greater
    than `L + np.array(Hmat.shape) - 1`. Also if `Amat` or `Hmat` aren't
    two-dimensional arrays, or if `y` doesn't have the correct size to store
    the output of the linear convolution.
    References
    ----------
    Wikipedia is only semi-unhelpful on this topic: see "Overlap-add method".
    """
    
    import multiprocessing.dummy as multiprocessing
    
    if amat_shape is None:
        amat_shape = Amat.shape
    if hmat_shape is None:
        hmat_shape = Hmat.shape
    
    M = np.array(hmat_shape)
    Na = np.array(amat_shape)
    
    ys = (amat_shape[0] + hmat_shape[0] - 1, amat_shape[1] + hmat_shape[1] - 1)
    
    if path is None:
        path = os.getcwd()
    
    if L is None:
        L = M * 100
    else:
        L = np.array(L)

    if Nfft is None:
        Nfft = 2 ** np.ceil(np.log2(L + M - 1)).astype(int)
    else:
        Nfft = np.array(Nfft, dtype=int)

    if not (np.all(L > 0) and L.size == 2):
        raise ValueError('L must have two positive elements')
    if not (np.all(Nfft >= L + M - 1) and Nfft.size == 2):
        msg = 'Nfft must have two elements >= L + M - 1 where M = Hmat.shape'
        raise ValueError(msg)
    if not (Amat.ndim <= 2):
        raise ValueError('Amat must be a 2D array')
    if hasattr(Hmat, 'ndim') and not (Hmat.ndim <= 2):
        raise ValueError('Hmat must be a 2D array')

    if isinstance(Hmat, GriddedPSFModel):
        Hf = Hmat
    else:
        Hf = fft2(Hmat, Nfft)
    
    pool = multiprocessing.Pool(processes=cores)
    m = multiprocessing.Manager()
    lock = m.Lock()
    print_lock = m.Lock()
    results = []
    logger.info("Starting job server with {} workers".format(pool._processes))
        
    (XDIM, YDIM) = (1, 0)
    adjust = lambda x: x                           # no adjuster
    start = [0, 0]
    endd = [0, 0]

    total_boxes = (Na[XDIM] // L[XDIM] + 1) * (Na[YDIM] // L[YDIM] + 1)
    percent_done = Percenter(total_boxes)
    def closing_log(pos):
        print_lock.acquire()
        if verbose and logger is not None:
            logger.info("Finishing box {}".format(pos))
        if verbose and state_setter is not None:
            percent_done.incr()
            state_setter(base_state + " {} done".format(percent_done.percent))
        print_lock.release()

    while start[XDIM] <= Na[XDIM]:
        endd[XDIM] = min(start[XDIM] + L[XDIM], Na[XDIM])
        start[YDIM] = 0
        while start[YDIM] <= Na[YDIM]:
            if verbose and logger is not None:
                logger.info("Starting box {}".format(start))
            endd[YDIM] = min(start[YDIM] + L[YDIM], Na[YDIM])
            thisend = np.minimum(Na + M - 1, start + Nfft)
            pos = (start[YDIM], endd[YDIM], start[XDIM], endd[XDIM], thisend[YDIM], thisend[XDIM])
            sub_arr = np.empty_like(Amat[start[YDIM]:endd[YDIM], start[XDIM]:endd[XDIM]])
            sub_arr[:,:] = Amat[start[YDIM]:endd[YDIM], start[XDIM]:endd[XDIM]]
            res = pool.apply_async(computation, args=(sub_arr, Hf, pos, Nfft, y, ys, adjust, lock, path, memmap), callback=closing_log)
            results.append(res)
            start[YDIM] += L[YDIM]
        start[XDIM] += L[XDIM]
    pool.close()
    pool.join()
#     for result in results:
#         logger.info("Result: {}".format(result.get()))
#         logger.info("Success: {}".format(result.successful()))


def overlapadd2(Amat, amat_shape, 
                Hmat, hmat_shape,
                L=None, Nfft=None, y=None, verbose=False, logger=None, 
                state_setter=None, base_state=""):
    """
    Fast two-dimensional linear convolution via the overlap-add method.
    The overlap-add method is well-suited to convolving a very large array,
    `Amat`, with a much smaller filter array, `Hmat` by breaking the large
    convolution into many smaller `L`-sized sub-convolutions, and evaluating
    these using the FFT. The computational savings over the straightforward
    two-dimensional convolution via, say, scipy.signal.convolve2d, can be
    substantial for large Amat and/or Hmat.
    Parameters
    ----------
    Amat, Hmat : array_like
        Two-dimensional input arrays to be convolved. For computational
        purposes, Amat should be larger.
    L : sequence of two ints, optional
        Length of sub-convolution to use. This should ideally be as large as
        possible to obtain maximum speedup: that is, if you have the memory to
        compute the linear convolution using FFT2, i.e., replacing spatial
        convolution with frequency-domain multiplication, then let `L =
        np.array(Amat.shape) + np.array(Hmat.shape) - 1`. Usually, though, you
        are considering overlap-add because you can't afford a batch transform
        on your arrays, so make `L` as big as you can.
    Nfft : sequence of two ints, optional
        Size of the two-dimensional FFT to use for each `L`-sized
        sub-convolution. If omitted, defaults to the next power-of-two, that
        is, the next power-of-two on or after `L + np.array(Hmat.shape) - 1`.
        If you choose this default, try to avoid `L` such that this minimum, `L
        + np.array(Hmat.shape) - 1`, is just a bit greater than a power-of-two,
        since you'll be paying for a 4x larger FFT than you need.
    y : array_like, optional
        Storage for the output. Useful if using a memory-mapped file, e.g.
    verbose : boolean, optional
        If True, prints a message for each `L`-sized subconvolution.
    Returns
    -------
    y : same as passed in, or ndarray if no `y` passed in
        The `np.array(Amat.shape) + np.array(Hmat.shape) - 1`-sized
        two-dimensional array containing the linear convolution. Should be
        within machine precision of, e.g., `scipy.signal.convolve2d(Amat,
        Hmat, 'full')`.
    Raises
    ------
    ValueError if `L` and `Nfft` aren't two-element, and too small: both
    elements of `L` must be greater than zero, and `Nfft`'s must be greater
    than `L + np.array(Hmat.shape) - 1`. Also if `Amat` or `Hmat` aren't
    two-dimensional arrays, or if `y` doesn't have the correct size to store
    the output of the linear convolution.
    References
    ----------
    Wikipedia is only semi-unhelpful on this topic: see "Overlap-add method".
    """
    if amat_shape is None:
        amat_shape = Amat.shape
    if hmat_shape is None:
        hmat_shape = Hmat.shape

    M = np.array(hmat_shape)
    Na = np.array(amat_shape)

    if y is None:
        y = np.zeros(M + Na - 1, dtype=Amat.dtype)
    elif y.shape != tuple(M + Na - 1):
        raise ValueError('y given has incorrect dimensions', M + Na - 1)

    if L is None:
        L = M * 100
    else:
        L = np.array(L)

    if Nfft is None:
        Nfft = 2 ** np.ceil(np.log2(L + M - 1)).astype(int)
    else:
        Nfft = np.array(Nfft, dtype=int)

    if not (np.all(L > 0) and L.size == 2):
        raise ValueError('L must have two positive elements')
    if not (np.all(Nfft >= L + M - 1) and Nfft.size == 2):
        raise ValueError('Nfft must have two elements >= L + M - 1 where M = Hmat.shape')
    if not (Amat.ndim <= 2):
        raise ValueError('Amat must be a 2D array')
    if hasattr(Hmat, 'ndim') and not (Hmat.ndim <= 2):
        raise ValueError('Hmat must be a 2D array')

    if isinstance(Hmat, GriddedPSFModel):
        Hf = Hmat
    else:
        Hf = fft2(Hmat, Nfft)

    (XDIM, YDIM) = (1, 0)
    start = [0, 0]
    endd = [0, 0]
    total_boxes = (Na[XDIM] // L[XDIM] + 1) * (Na[YDIM] // L[YDIM] + 1)
    current_box = 0
    while start[XDIM] <= Na[XDIM]:
        endd[XDIM] = min(start[XDIM] + L[XDIM], Na[XDIM])
        start[YDIM] = 0
        while start[YDIM] <= Na[YDIM]:
            if verbose and logger is not None:
                logger.info("Starting box {}".format(start))
            if verbose and state_setter is not None:
                new_state = base_state + " {:.2f}% done"
                state_setter(new_state.format((current_box/total_boxes)*100.))
            endd[YDIM] = min(start[YDIM] + L[YDIM], Na[YDIM])
            thisend = np.minimum(Na + M - 1, start + Nfft)
            Asub = Amat[start[YDIM]:endd[YDIM], start[XDIM]:endd[XDIM]]
            Af = fft2(Asub, Nfft)
            
            if isinstance(Hf, GriddedPSFModel):
                # Generate PSF
                yg, xg = np.mgrid[start[YDIM]:endd[YDIM],start[XDIM]:endd[XDIM]]
                y_0 = (start[YDIM]+endd[YDIM])//2
                x_0 = (start[XDIM]+endd[XDIM])//2
                # Make a grid that's the size of the PSF fov_pix around its centre.
                ys2, xs2 = hmat_shape[0]//2, hmat_shape[1]//2
                yg, xg = np.mgrid[y_0-ys2:y_0+ys2, x_0-xs2:x_0+xs2]
                psf = Hf.evaluate(x=xg, y=yg, flux=1, x_0=x_0, y_0=y_0)
                yt = np.real(ifft2(fft2(psf, Nfft) * Af))
            else:
                yt = np.real(ifft2(Hf * Af))
            ys = yt[:(thisend[YDIM]-start[YDIM]), :(thisend[XDIM]-start[XDIM])]
            y[start[YDIM]:thisend[YDIM], start[XDIM]:thisend[XDIM]] += ys[:,:]
            start[YDIM] += L[YDIM]
            current_box += 1
        start[XDIM] += L[XDIM]
    return y


def test():
    from scipy.signal import convolve2d
    A = np.random.randn(33, 55)
    H = np.random.randn(4, 5)
    gold = convolve2d(A, H)
    assert(np.allclose(gold, overlapadd2(A, H, L=[12, 12])))
    assert(np.allclose(gold, overlapadd2(A, H, L=[12, 120])))
    assert(np.allclose(gold, overlapadd2(A, H, L=[90, 120])))

    assert(np.allclose(gold, overlapadd2(H, A, L=[190, 220])))
    assert(np.allclose(gold, overlapadd2(H, A, L=[1, 1])))

    assert(np.allclose(gold, overlapadd2(H, A)))
    assert(np.allclose(gold, overlapadd2(A, H)))

    assert(np.allclose(convolve2d(A.T, H.T),
                       overlapadd2(A.T, H.T, L=[190, 220])))

    A = np.random.randn(33, 55) + 1j * np.random.randn(33, 55)
    H = np.random.randn(4, 5) + 1j * np.random.randn(4, 5)
    gold = convolve2d(A, H)
    assert(np.allclose(gold, overlapadd2(H, A)))
    assert(np.allclose(gold, overlapadd2(A, H)))

    print('Test passed!')


if __name__ == '__main__':
    test()

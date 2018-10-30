"""
General CGI form functions.

:Author: Pey Lian Lim

:Organization: Space Telescope Science Institute

"""
from __future__ import absolute_import,division

# External modules
import importlib, inspect, os, shutil, socket, sys, urllib, uuid
import numpy as np
import astropy.io.fits as pyfits
import multiprocessing.dummy as multiprocessing
from numpy.fft import fft2, ifft2
from astropy.io import ascii
from astropy.table import Table
from jwst_backgrounds.jbt import background


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
    def __init__(self, fname, shape, mode='r+'):
        self.fp = np.memmap(fname, dtype='float32', mode=mode, shape=shape)
    
    def __enter__(self):
        return self.fp
    
    def __exit__(self, exc_type, exc_value, traceback):
        del self.fp

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
            self.cache_url = urllib.request.pathname2url(os.path.join(self.cache_path, "remote_cache/"))
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

#-----------
def GetStipsData(to_retrieve):
    """
    Retrieve a file from the stips_data directory. Will also print out a warning if the directory
    can't be found.
    """
    local_data_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "data"))
    if "stips_data" not in os.environ:
        sys.stderr.write("WARNING: stips_data environment variable not found. Falling back on local STIPS data.\n")
        sys.stderr.write("WARNING: STIPS local data may be older than data available via the stips_data environment variable.\n")
    stips_data_base = os.environ.get("stips_data", local_data_dir)
    if not os.path.exists(stips_data_base):
        sys.stderr.write("ERROR: stips_data directory at {} not found. STIPS requires the stips_data directory to function correctly.\n".format(stips_data_base))
        sys.stderr.write("ERROR: Please make sure that the stips_data environment variable exists and points to the location of stips_data.\n")
    retrieval_file = os.path.join(stips_data_base, to_retrieve)
    if not os.path.exists(retrieval_file):
        sys.stderr.write("ERROR: STIPS data file {} not found. STIPS requires the stips_data directory to function correctly.\n".format(retrieval_file))
        sys.stderr.write("ERROR: Please make sure that the stips_data environment variable exists and points to the location of stips_data.\n")
        sys.stderr.write("ERROR: Please try downloading the stips_data directory again.\n")
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



def computation(arr, Hf, pos, Nfft, y, ys, adjust, lock, path):
    start_y, end_y, start_x, end_x, thisend_y, thisend_x = pos
    conv = adjust(ifft2(Hf * fft2(arr, Nfft)))        
    lock.acquire()
    with ImageData(y, ys) as dat:
        dat[start_y:thisend_y, start_x:thisend_x] += (conv[:(thisend_y-start_y), :(thisend_x-start_x)])
    lock.release()
    return "[{}, {}]".format(start_y, start_x)


def overlapaddparallel(Amat, Hmat, L=None, Nfft=None, y=None, verbose=False, logger=None, state_setter=None, base_state="", path=None):
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
    
    M = np.array(Hmat.shape)
    Na = np.array(Amat.shape)
    
    ys = (Amat.shape[0]+Hmat.shape[0]-1, Amat.shape[1]+Hmat.shape[1]-1)
    
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
        raise ValueError('Nfft must have two elements >= L + M - 1 where M = Hmat.shape')
    if not (Amat.ndim <= 2 and Hmat.ndim <= 2):
        raise ValueError('Amat and Hmat must be 2D arrays')

    Hf = fft2(Hmat, Nfft)
    
    pool = multiprocessing.Pool()
    m = multiprocessing.Manager()
    lock = m.Lock()
    print_lock = m.Lock()
    results = []
    logger.info("Starting job server with {} workers".format(pool._processes))
        
    (XDIM, YDIM) = (1, 0)
    adjust = lambda x: x                           # no adjuster
    if np.isrealobj(Amat) and np.isrealobj(Hmat):  # unless inputs are real
        adjust = np.real                           # then ensure real
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
            res = pool.apply_async(computation, args=(sub_arr, Hf, pos, Nfft, y, ys, adjust, lock, path), callback=closing_log)
            results.append(res)
            start[YDIM] += L[YDIM]
        start[XDIM] += L[XDIM]
    pool.close()
    pool.join()
#     for result in results:
#         logger.info("Result: {}".format(result.get()))
#         logger.info("Success: {}".format(result.successful()))


def overlapadd2(Amat, Hmat, L=None, Nfft=None, y=None, verbose=False, logger=None, state_setter=None, base_state=""):
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
    M = np.array(Hmat.shape)
    Na = np.array(Amat.shape)

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
    if not (Amat.ndim <= 2 and Hmat.ndim <= 2):
        raise ValueError('Amat and Hmat must be 2D arrays')

    Hf = fft2(Hmat, Nfft)

    (XDIM, YDIM) = (1, 0)
    adjust = lambda x: x                           # no adjuster
    if np.isrealobj(Amat) and np.isrealobj(Hmat):  # unless inputs are real
        adjust = np.real                           # then ensure real
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
                state_setter(base_state + " {:.2f}% done".format((current_box/total_boxes)*100.))
            endd[YDIM] = min(start[YDIM] + L[YDIM], Na[YDIM])
            thisend = np.minimum(Na + M - 1, start + Nfft)
            yt = adjust(ifft2(Hf * fft2(Amat[start[YDIM] : endd[YDIM], start[XDIM] : endd[XDIM]], Nfft)))
            y[start[YDIM] : thisend[YDIM], start[XDIM] : thisend[XDIM]] += (yt[:(thisend[YDIM] - start[YDIM]), :(thisend[XDIM] - start[XDIM])])
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

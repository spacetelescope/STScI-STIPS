"""
Unit conversion functions.

*THIS IS FOR CGI ONLY AND DOES NOT HAVE SYNPHOT*

.. note:: Jansky = 1E-23 erg/s/cm^2/Hz

:Author: Pey Lian Lim

:Organization: Space Telescope Science Institute

:History:
    * 2010/11/23 PLL created this module.
    * 2011/02/18 PLL added `Convert2Jy` and `Convert2ergs`.
    * 2011/06/02 PLL modified `GetUresp`.
    * 2011/08/17 PLL added `RadiusParsec2pix`, `DistanceModulus` and `RescaleArray`.
    * 2011/11/02 PLL added function to calculate expected galaxy counts.
    * 2011/12/01 PLL added function to get NIRCam mode based on filter.
    * 2012/02/14 PLL removed Synphot for CGI.
    * 2014/12/03 BAQ removed old un-needed functions as part of moving to instrument classes.

"""
import numpy as np

def RadiiUnknown2Arcsec(radii, rad_units, distances):
    """
    Given offsets from zero in some unit system (either parsecs or arcseconds),
    each with an attached distance, convert those offsets into arcseconds.

    Parameters
    ----------
    radii: array_like
        Radii of objects in unknown units.
    
    rad_units: string
        Units of the radius. Value will be one of:
            - pc
            - arcsec

    distances: array_like
        Distances of stars in parsec.

    Returns
    -------
    arcsecs: array_like
        Radii in arcseconds.
    """
    # Calculate angular size in parcsecs
    if rad_units == "arcsec": #no conversion needed
        arcsecs = radii
    elif rad_units == "pc":
        r = np.arctan(radii/distances) # rad
        arcsecs = r * 0.5 * 1296000.0 / np.pi # rad to arcsec
    
    return arcsecs

def RadiiUnknown2Parsec(radii, rad_units, distances):
    """
    Given offsets from zero in some unit system (either parsecs or arcseconds),
    each with an attached distance, convert those offsets into arcseconds.

    Parameters
    ----------
    radii: array_like
        Radii of objects in unknown units.
    
    rad_units: string
        Units of the radius. Value will be one of:
            - pc
            - arcsec

    distances: array_like
        Distances of stars in parsec.

    Returns
    -------
    parsecs: array_like
        Radii in parsecs.
    """
    # Calculate angular size in parcsecs
    if rad_units == "pc": #no conversion needed
        parsecs = radii
    elif rad_units == "arcsec":
        parsecs = np.tan(2*np.pi*radii/1296000.0)*distances
    
    return parsecs

#-----------
def RescaleArray(arr, newLo, newHi):
    """
    Rescale array to new limits.

    Parameters
    ----------
    arr: array_like
        Data to rescale.

    newLo, newHi: float
        New lower and upper limits.

    Returns
    -------
    newArr: array_like
        Rescale `arr`.
        
    """
    a_min, a_max = arr.min().astype('float'), arr.max().astype('float')
    frac = (arr - a_min) / (a_max - a_min)
    newArr = frac * (newHi - newLo) + newLo
    return newArr

#-----------
def DivideInterval(divisor,reverse=False):
    """
    Takes an interval expressed as a string in the form "low,high,interval" and
        then performs black magic to turn it into a linear array.
        
        In particular, depending on the value of interval:
        
            interval is an integer X:
            interval is 'n' followed by an integer X:
            
                Use np.linspace(low,high,interval) to produce X intervals,
                    including both low and high, and evenly divided between them.
            
            interval is 'i' followed by an integer X:
            
                Use np.linspace(low,high,((high-low)/interval)+1) to produce
                    intervals *of* X from low to high, and including both. Uses
                    'linspace' rather than 'arange' because 'linspace' does not
                    produce weird effects at start and end points, and will
                    actually make zero 0. rather than 1.58395e-15 or some such.
            
            interval is 'd':
            
                Use np.linspace() called multiple times to produce decades.
                    For example, '1,10,d' would produce 1,2,3,...,10, whilst
                    '1,100,d' would produce 1,2,3,...,9,10,20,30,...,90,100.
                    '1,150,d' would produce 1,2,3,...,9,10,20,30,...,90,100,110,...,140,150
                    By analogy up succeeding powers of 10. Note that results 
                    will *always* include high, and may be odd if the first
                    number is not a power of 10. For example,
                    '2,250,d' produces 2,4,6,...,18,20,40,60,...,180,200,220,250.
                    So be careful and, if you want better behaviour, program
                    it yourself.
            
            interval is 'd5':
            
                Use np.linspace() called multiple times to produce half-decades.
                    Exactly as before, except that it includes half-intervals, so
                    '1,100,d5' would give 1,1.5,2,2.5,3,...,9.5,10,15,20,25,...,95,100.
                    Same caveats apply as before. If you want d with an arbitrary number
                    following it, then, again, program it yourself.
    """
    items = divisor.split(",")
    low = float(items[0])
    high = float(items[1])
    interval = items[2]
    if interval == "d": #divide by decades
        interval = np.array(())
        start = low
        while start <= high:
            end = start * 10
            intervals = 10
            if end > high:
                intervals = 91
            arr = np.linspace(start,end,intervals)
            if end > high:
                arr = arr[np.where(arr<=high)]
            interval = np.append(interval,arr[:-1])
            start = end
        interval = np.append(interval,high)
    elif interval == "d5": #divide by half-decades
        interval = np.array(())
        start = low
        while start <= high:
            end = start * 10
            intervals = 19
            if end > high:
                intervals = 181
            arr = np.linspace(start,end,intervals)
            if end > high:
                arr = arr[np.where(arr<=high)]
            interval = np.append(interval,arr[:-1])
            start = end
        interval = np.append(interval,high)
    elif interval[0] == "n": #Number of intervals follows the n
        interval = np.linspace(low,high,int(interval[1:]))
    elif interval[0] == "i": #Interval value follows the n.
        #Note: this is a bit hacky, but linspace gives overall a better set.
        interval = np.linspace(low,high,int(round((high-low)/float(interval[1:])))+1)
    else: #assume number of intervals, but with no 'n' (backwards compatibility)
        interval = np.linspace(low,high,int(interval))
    if reverse:
        interval = np.flipud(interval)
    return interval

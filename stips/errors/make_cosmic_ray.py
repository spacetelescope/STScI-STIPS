"""
Functions to simulate cosmic rays.

:Authors: Pey Lian Lim (Python); Mike Regan (IDL)

:Organization: Space Telescope Science Institute

:History:
    * 2010/08/17 PLL converted from IDL to Python.
    * 2014/12/03 BAQ moved GaussPsf2D into this file, as this file is the only place it's used.
"""

# External modules
import numpy as np

#-----------
def GaussPsf2D(npix, fwhm, normalize=True):
    """
    Parameters
    ----------
    npix: int
        Number of pixels for each dimension.
        Just one number to make all sizes equal.

    fwhm: float
        FWHM (pixels) in each dimension.
        Single number to make all the same.

    normalize: bool, optional
        Normalized so total PSF is 1.

    Returns
    -------
    psf: array_like
        Gaussian point spread function.
    """

    # Initialize PSF params
    cntrd = (npix - 1.0) * 0.5
    st_dev = 0.5 * fwhm / np.sqrt( 2.0 * np.log(2) )

    # Rene Breton 2011-10-20
    # https://groups.google.com/group/astropy-dev/browse_thread/thread/5ee6cd662236e382
    x, y = np.indices([npix,npix]) - (npix-1)*0.5
    psf = np.exp( -0.5 * ((x**2 + y**2)/st_dev**2) ) 
        
    # Normalize
    if normalize: psf /= psf.sum()

    return psf

#-----------
def MakeCosmicRay(xSize, ySize, crProb, crElectrons, crSize, crPsf, seed, verbose=True):
    """
    Simulate cosmic rays.

    Parameters
    ----------
    xSize: int
        Output X size (pix).

    ySize: int
        Output Y size (pix).

    crProb: float
        Probability of CR.

    crElectrons: float
        CR energy (e-).

    crSize: int
        Size of CR (pix) from `GetCrShape`.

    crPsf: array_like
        PSF of CR from `GetCrShape`.

    verbose: bool, optional
        Print info to screen.

    Returns
    -------
    crArray: array_like
        CR to populate (e-).
    """
    crArray = np.zeros((ySize,xSize))

    # See if a CR hit
    crPoisson = np.random.RandomState(seed=seed).poisson(lam=crProb, size=xSize*ySize)
    cr_locs = np.where(crPoisson >= 1)
    num_crs = len(cr_locs[0])
    if num_crs == 0:
        if verbose: 
            print('No CR in poisson')
        return crArray

    # CR locations
    cr_x = cr_locs[0] % xSize
    cr_y = cr_locs[0] / ySize

    # Only want full CR on detector
    crXMin = cr_x - crSize
    crXMax = cr_x + crSize
    crYMin = cr_y - crSize
    crYMax = cr_y + crSize
    idxGood = np.where((crXMin > 0) & (crXMax < xSize) & (crYMin > 0) & (crYMax < ySize))
    numGood = len(idxGood[0])
    if numGood == 0:
        if verbose: 
            print('No CR fully inside frame')
        return crArray

    # CR energy in e-
    cr_electrons = crElectrons + np.sqrt(crElectrons) * np.random.RandomState(seed=seed).randn(numGood)

    # Populate CR
    for i in range(numGood):
        ii = idxGood[0][i]
        x1 = crXMin[ii]
        x2 = crXMax[ii] + 1
        y1 = crYMin[ii]
        y2 = crYMax[ii] + 1
        cr_sim = cr_electrons[i] * crPsf
        crArray[int(y1):int(y2),int(x1):int(x2)] += cr_sim
    # End of i loop

    return crArray

#-----------
def GetCrTemplate(fwhm=0.9):
    """
    Define cosmic rays size and PSF.

    Parameters
    ----------
    fwhm: float, optional
        FWHM of cosmic rays (pix).

    Returns
    -------
    cr_size: int
        Cosmic rays size (pix).

    cr_psf: array_like
        Cosmic rays gaussian PSF.
    """
    cr_size = int(2 * fwhm)
    cr_psf = GaussPsf2D(2*cr_size+1, fwhm)
    return cr_size, cr_psf

#-----------
def GetCrProbs(rates, pixarea, exptime):
    """
    Calculate probabilities of CR hits.

    Parameters
    ----------
    rates: tuple of float
        CR hit rates (hits/cm^2/s).

    pixarea: float
        Pixel area (cm^2).

    exptime: float
        Exposure time (s).

    Returns
    -------
    probs: tuple of float
        Probabilities for respective elements in `rates`.
    """
    m = pixarea * exptime # cm^2 s
    probs = []
    for r in rates: probs.append( r * m) # hits
    return tuple(probs)

"""
This file contains all the functions required to create the ePSF. Including creating the
ePSF grid, IPC correction, and functions to place stars in the output images.

ePSF functions translated here were originally written by Andrea Bellini in Fortran.

Current assumptions:
- The input PSF is square, the size of the array x == y.
- The ePSF will be oversampled by a fixed factor of 4.
- The output image is square.

Pending:
- Write logging functions in.
- Implement error catching.
- make_epsf() and real_psf() have not been optimized.
- Write tests for all functions.

:Author: Sebastian Gomez
:Organization: Space Telescope Science Institute
"""

# External modules
import numpy as np
from scipy import ndimage

# Convolution Constants
# Inter Pixel Capacitance
IPC = np.array([[0.21,  1.62,  0.20],
                [1.88, 91.59,  1.87],
                [0.21,  1.66,  0.22]]) / 100.0
# Convolution array for creating a 4x upscaled ePSF, resulting
# in the summ of the inner 3 pixels and a corresponding fraction
# of the edges and corners to create an effective 4x4 pixels.
EPSF4 = np.array([[0.25, 0.5, 0.5, 0.5, 0.25],
                  [0.50, 1.0, 1.0, 1.0, 0.50],
                  [0.50, 1.0, 1.0, 1.0, 0.50],
                  [0.50, 1.0, 1.0, 1.0, 0.50],
                  [0.25, 0.5, 0.5, 0.5, 0.25]])

# The ePSF is *always* 4x upscaled
PSF_UPSCALE = 4

# There are 3 PSF sizes.
#   - general PSF: 44 pixels
#   - bright source PSF:
#   - very bright source PSF:
PSF_BOXSIZE = 44
PSF_BRIGHT_BOXSIZE = 88
PSF_EXTRA_BRIGHT_BOXSIZE = 176

# The PSF grid is always 3x3
PSF_GRID_SIZE = 3


def make_epsf(psf_in):
    """
    Create an ePSF from an input array of a specific detector
    corner and filter. The PSF will be convolved with an effective
    filter of 4x4 pixels, and then the IPC calculation applied to
    each pixel.

    Parameters
    ----------
    psf_in : numpy.ndarray
        2D array with PSF image

    Returns
    -------
    psf_final : numpy.ndarray
        2D array with ePSF of same shape as psf_in

    """

    # Empty array to fill later
    psf_out = np.zeros_like(psf_in)

    # Assuming the PSF is square
    size = psf_in.shape[0]

    # Create ePSF
    psf_out = ndimage.convolve(psf_in, EPSF4, mode='constant', cval=0.0)

    # Apply correct scaling to edges of the image, currently done in
    # individual steps for computational reasons.
    # Multiply columns x16
    psf_out[0:size,      0] = psf_in[0:size,      0]*16.0
    psf_out[0:size,      1] = psf_in[0:size,      1]*16.0
    psf_out[0:size,      2] = psf_in[0:size,      2]*16.0
    psf_out[0:size, size-1] = psf_in[0:size, size-1]*16.0
    psf_out[0:size, size-2] = psf_in[0:size, size-2]*16.0
    psf_out[0:size, size-3] = psf_in[0:size, size-3]*16.0
    # Multiply rows x16, while excluding the corners that were
    # already multiplied in the previous step
    psf_out[0, 3:(size-PSF_UPSCALE)+1] = psf_in[0, 3:(size-PSF_UPSCALE)+1]*16.0
    psf_out[1, 3:(size-PSF_UPSCALE)+1] = psf_in[1, 3:(size-PSF_UPSCALE)+1]*16.0
    psf_out[2, 3:(size-PSF_UPSCALE)+1] = psf_in[2, 3:(size-PSF_UPSCALE)+1]*16.0
    psf_out[size-1, 3:(size-PSF_UPSCALE)+1] = psf_in[size-1, 3:(size-PSF_UPSCALE)+1]*16.0
    psf_out[size-2, 3:(size-PSF_UPSCALE)+1] = psf_in[size-2, 3:(size-PSF_UPSCALE)+1]*16.0
    psf_out[size-3, 3:(size-PSF_UPSCALE)+1] = psf_in[size-3, 3:(size-PSF_UPSCALE)+1]*16.0

    # Empty Array
    psf_final = np.zeros_like(psf_in)

    # Do calculation for every pixel in PSF image
    x_range = range(PSF_UPSCALE, size-PSF_UPSCALE-1)
    y_range = range(PSF_UPSCALE, size-PSF_UPSCALE-1)

    # Add IPC to output PSF
    for idx in x_range:
        for idy in y_range:
            pcen = psf_out[idy, idx]

            psf_final[idy-PSF_UPSCALE, idx-PSF_UPSCALE] += pcen*IPC[0, 0]
            psf_final[idy-PSF_UPSCALE,             idx] += pcen*IPC[0, 1]
            psf_final[idy-PSF_UPSCALE, idx+PSF_UPSCALE] += pcen*IPC[0, 2]

            psf_final[idy, idx-PSF_UPSCALE] += pcen*IPC[1, 0]
            psf_final[idy,             idx] += pcen*IPC[1, 1]
            psf_final[idy, idx+PSF_UPSCALE] += pcen*IPC[1, 2]

            psf_final[idy+PSF_UPSCALE, idx-PSF_UPSCALE] += pcen*IPC[2, 0]
            psf_final[idy+PSF_UPSCALE,             idx] += pcen*IPC[2, 1]
            psf_final[idy+PSF_UPSCALE, idx+PSF_UPSCALE] += pcen*IPC[2, 2]

    return psf_final


def bicubic(epsf, iy, ix, fx, fy):
    """
    Perform the bi-cubic interpolation of an image from a center
    pixel and a pixel phase. Functions originally obtained from
    Andrea Bellini's Fortran code.

    This function is significantly faster than the more complex
    scipy.interpolate.interp2d, which differs from this at the level
    of ~0.5%

    Parameters
    ----------
    epsf : numpy.ndarray
        2D array with ePSF image
    ix : int
        Reference x pixel
    iy : int
        Reference y pixel
    fx : int
        Pixel phase long x direction
    fy : int
        Pixel phase long y direction

    Returns
    -------
    rpsf_phot : float
        Output value of fractional light
    """

    # Lower Left Value
    A1 = (epsf[iy,     ix])
    B1 = (epsf[iy,   ix+1]-epsf[iy,   ix-1])/2
    C1 = (epsf[iy+1,   ix]-epsf[iy-1,   ix])/2
    D1 = (epsf[iy,   ix+1]+epsf[iy,   ix-1]-2*A1)/2
    E1 = (epsf[iy+1, ix+1]-A1)
    F1 = (epsf[iy+1,   ix]+epsf[iy-1,   ix]-2*A1)/2
    V1 = (A1
          + B1*(fx)
          + C1*(fy)
          + D1*(fx)**2
          + E1*(fx)*(fy)
          + F1*(fy)**2)

    # Lower Right Value
    A2 = (epsf[iy,   ix+1])
    B2 = (epsf[iy,   ix+2]-epsf[iy,     ix])/2
    C2 = (epsf[iy+1, ix+1]-epsf[iy-1, ix+1])/2
    D2 = (epsf[iy,   ix+2]+epsf[iy,     ix]-2*A2)/2
    E2 = -(epsf[iy+1,   ix]-A2)
    F2 = (epsf[iy+1, ix+1]+epsf[iy-1, ix+1]-2*A2)/2
    V2 = (A2
          + B2*(fx-1)
          + C2*(fy)
          + D2*(fx-1)**2
          + E2*(fx-1)*(fy)
          + F2*(fy)**2)

    # Upper Left Value
    A3 = (epsf[iy+1,   ix])
    B3 = (epsf[iy+1, ix+1]-epsf[iy+1, ix-1])/2
    C3 = (epsf[iy+2,   ix]-epsf[iy,     ix])/2
    D3 = (epsf[iy+1, ix+1]+epsf[iy+1, ix-1]-2*A3)/2
    E3 = -(epsf[iy,   ix+1]-A3)
    F3 = (epsf[iy+2,   ix]+epsf[iy,     ix]-2*A3)/2
    V3 = (A3
          + B3*(fx)
          + C3*(fy-1)
          + D3*(fx)**2
          + E3*(fx)*(fy-1)
          + F3*(fy-1)**2)

    # Upper Right Value
    A4 = (epsf[iy+1, ix+1])
    B4 = (epsf[iy+1, ix+2]-epsf[iy+1,   ix])/2
    C4 = (epsf[iy+2, ix+1]-epsf[iy,   ix+1])/2
    D4 = (epsf[iy+1, ix+2]+epsf[iy+1,   ix]-2*A4)/2
    E4 = (epsf[iy,     ix]-A4)
    F4 = (epsf[iy+2, ix+1]+epsf[iy,   ix+1]-2*A4)/2
    V4 = (A4
          + B4*(fx-1)
          + C4*(fy-1)
          + D4*(fx-1)**2
          + E4*(fx-1)*(fy-1)
          + F4*(fy-1)**2)

    rpsf_phot = ((1-fx)*(1-fy)*V1
                 + (fx)*(1-fy)*V2
                 + (1-fx)*(fy)*V3
                 + (fx)*(fy)*V4)

    return rpsf_phot


def real_psf(dx, dy, epsf, psf_center=177, boxsize=PSF_BOXSIZE):
    """
    Calculate the fraction of light from a epsf that should
    fall on a given pixel (x, y). This function assumes the
    input ePSF was oversampled by a factor of 4.

    Parameters
    ----------
    dx : float
        relative location of source along x within boxsize
        (0) is the center of boxsize.
    dy : float
        relative location of source along y within boxsize
        (0) is the center of boxsize.
    epsf : numpy.ndarray
        2D array with ePSF image
    psf_center : int
        center of the input psf model
    boxsize : int
        size of PSF box.

    Returns
    -------
    rpsf_phot : float
        Output value of fractional light
    """

    # Relative location of pixel within ePSF (for a factor of 4 oversample)
    rx = psf_center + dx*4
    ry = psf_center + dy*4
    ix = int(rx)
    iy = int(ry)
    fx = rx-ix  # Pixel Phase
    fy = ry-iy  # Pixel Phase
    dd = np.sqrt(dx**2+dy**2)
    rpsf_phot = 0.
    # If the pixel location is outside the boxsize return 0
    if (np.abs(dx) > boxsize) or (np.abs(dy) > boxsize):
        return rpsf_phot
    # Else, this will be the case for most of the pixels
    if (dd > 4.0):
        # Bi-linear interpolation
        # Weigh by the area substended by the 4 adjacent
        # pixels to the reference pixel in question.
        rpsf_phot = ((1-fx)*(1-fy)*epsf[iy,    ix]
                     + (fx)*(1-fy)*epsf[iy,  ix+1]
                     + (1-fx)*(fy)*epsf[iy+1,  ix]
                     + (fx)*(fy)*epsf[iy+1, ix+1])
        return rpsf_phot

    # Do Bi-cubic interpolation for the innermost pixels
    rpsf_phot = bicubic(epsf, iy, ix, fx, fy)

    return rpsf_phot


def place_source(xpix, ypix, flux, image, epsf, boxsize=PSF_BOXSIZE, psf_center=177):
    """
    Place a source into image.

    Parameters
    ----------
    xpix : float
        x location of source
    ypix : float
        y location of source
    flux : float
        Flux of the source
    image : numpy.ndarray
        Empty (or not) 2D array of image where to place sources
    epsf : numpy.ndarray
        2D array with ePSF

    Returns
    -------
    image : numpy.ndarray
        2D array of image where sources were placed
    """

    # Get image size, assuming square
    image_size = image.shape[0]

    # Place if source is within image and boxise
    if (xpix > 0) & (xpix < image_size) & (ypix > 0) & (ypix < image_size):
        # Generate a box around the location of the source with size boxsize
        max_y = min(round(ypix+boxsize), image_size)
        min_y = max(0, round(ypix-boxsize))
        max_x = min(round(xpix+boxsize), image_size)
        min_x = max(0, round(xpix-boxsize))
        # Apply real_psf for every pixel in the box, dx/dy are the pixel
        # positions within the box.
        for j in range(min_y, max_y):
            dy = j-ypix
            for i in range(min_x, max_x):
                dx = i-xpix
                ff = real_psf(dx, dy, epsf, psf_center=psf_center, boxsize=boxsize)
                ffa = (flux*ff)
                image[j, i] += ffa

    return image


def interpolate_epsf(xpix, ypix, psf_array, image_size):
    """
    Interpolate the input ePSFs at the location of a specified
    source.

    Parameters
    ----------
    xpix : float
        x location of source
    ypix : float
        y location of source
    psf_array : list
        3 x 3 list with input ePSFs
    image_size : int
        Image size in pixels

    Returns
    -------
    epsf : np.array
        Interpolated ePSF at the location xpix, ypix
    """

    # Pixel location of half image
    half_image = round(image_size/2)

    # If star is in Lower Left Quadrant
    if (xpix <= half_image and ypix <= half_image):
        xf = ((xpix+4)/(half_image+4))
        yf = ((ypix+4)/(half_image+4))
        epsf = (xf*yf*psf_array[1][1] +
                (1-xf)*(1-yf)*psf_array[0][0] +
                (xf)*(1-yf)*psf_array[1][0] +
                (1-xf)*(yf)*psf_array[0][1])

    # If star is in Lower Right Quadrant
    elif (xpix > half_image and ypix <= half_image):
        xf = ((xpix+4-(half_image+4))/(half_image+4))
        yf = ((ypix+4)/(half_image+4))
        epsf = (xf*yf*psf_array[2][1] +
                (1-xf)*(1-yf)*psf_array[1][0] +
                (xf)*(1-yf)*psf_array[2][0] +
                (1-xf)*(yf)*psf_array[1][1])

    # If star is in Upper Left Quadrant
    elif (xpix <= 1.0*half_image and ypix > half_image):
        xf = ((xpix+4)/(half_image+4))
        yf = ((ypix+4-(half_image+4))/(half_image+4))
        epsf = (xf*yf*psf_array[1][2] +
                (1-xf)*(1-yf)*psf_array[0][1] +
                (xf)*(1-yf)*psf_array[1][1] +
                (1-xf)*(yf)*psf_array[0][2])

    # If star is in Upper Right Quadrant
    elif (xpix > half_image and ypix > half_image):
        xf = ((xpix+4-(half_image+4))/(half_image+4))
        yf = ((ypix+4-(half_image+4))/(half_image+4))
        epsf = (xf*yf*psf_array[2][2] +
                (1-xf)*(1-yf)*psf_array[1][1] +
                (xf)*(1-yf)*psf_array[2][1] +
                (1-xf)*(yf)*psf_array[1][2])

    return epsf

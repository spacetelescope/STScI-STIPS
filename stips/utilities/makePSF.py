"""
This file contains all the functions required to create the ePSF, intended to replace the
original GriddedPSFModel in STIPS. Including creating the ePSF grid, IPC correction, and
functions to place stars in the output images. No oversampling is happening here, so the 
oversample keyword is meaningless.

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
- Allow sources in the edges of the image (those break the code).
- Generate Webb PSFs as opposed to having this fixed one.

:Author: Sebastian Gomez
:Organization: Space Telescope Science Institute
"""

from __future__ import absolute_import,division

# External modules
from astropy.io import fits
import numpy as np
from scipy import ndimage
#import glob, importlib, inspect, os, requests, shutil, socket, struct, sys, tarfile 
#import urllib, uuid, yaml
#import numpy as np
#import astropy.io.fits as pyfits
#from astropy.io import ascii
#from astropy.table import Table
#from jwst_backgrounds.jbt import background
#from numpy.fft import fft2, ifft2
#from pathlib import Path
#from photutils.psf.models import GriddedPSFModel

#from .. import __version__ as __stips__version__

rind = lambda x : np.round(x).astype(int)

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

    # Create ePSF by summing inner 3 pixels and corresponding fraction of
    # edges and corners to create an effective 4x4 pixels.
    epsf_array = np.array([[0.25, 0.5 , 0.5 , 0.5 , 0.25],
                           [0.5 , 1.  , 1.  , 1.  , 0.5 ],
                           [0.5 , 1.  , 1.  , 1.  , 0.5 ],
                           [0.5 , 1.  , 1.  , 1.  , 0.5 ],
                           [0.25, 0.5 , 0.5 , 0.5 , 0.25]])

    psf_out = ndimage.convolve(psf_in, epsf_array, mode='constant', cval=0.0)

    # Apply correct scaling to edges of the image
    #for idy in range(0:size-1):
    psf_out[0:size,     0]=psf_in[0:size,0     ]*16.0
    psf_out[0:size,     1]=psf_in[0:size,1     ]*16.0
    psf_out[0:size,     2]=psf_in[0:size,2     ]*16.0
    psf_out[0:size,size-1]=psf_in[0:size,size-1]*16.0
    psf_out[0:size,size-2]=psf_in[0:size,size-2]*16.0
    psf_out[0:size,size-3]=psf_in[0:size,size-3]*16.0

    #for idx in range(3, 357):
    psf_out[0     ,3:(size-4)+1]=psf_in[0     ,3:(size-4)+1]*16.0
    psf_out[1     ,3:(size-4)+1]=psf_in[1     ,3:(size-4)+1]*16.0
    psf_out[2     ,3:(size-4)+1]=psf_in[2     ,3:(size-4)+1]*16.0
    psf_out[size-1,3:(size-4)+1]=psf_in[size-1,3:(size-4)+1]*16.0
    psf_out[size-2,3:(size-4)+1]=psf_in[size-2,3:(size-4)+1]*16.0
    psf_out[size-3,3:(size-4)+1]=psf_in[size-3,3:(size-4)+1]*16.0

    # Add IPC to output PSF
    ipc_array = np.array([[ 0.21,  1.62,  0.2 ],
                          [ 1.88, 91.59,  1.87],
                          [ 0.21,  1.66,  0.22]]) / 100.0

    # Empty Array
    psf_final = np.zeros_like(psf_in)

    # Do calculation for every pixel in PSF image
    x_range = range(4, size-5)
    y_range = range(4, size-5)

    for idx in x_range:
        for idy in y_range:
            pcen = psf_out[idy, idx]

            psf_final[idy-4,idx-4] += pcen*ipc_array[0,0]
            psf_final[idy-4,idx  ] += pcen*ipc_array[0,1]
            psf_final[idy-4,idx+4] += pcen*ipc_array[0,2]

            psf_final[idy  ,idx-4] += pcen*ipc_array[1,0]
            psf_final[idy  ,idx  ] += pcen*ipc_array[1,1]
            psf_final[idy  ,idx+4] += pcen*ipc_array[1,2]

            psf_final[idy+4,idx-4] += pcen*ipc_array[2,0]
            psf_final[idy+4,idx  ] += pcen*ipc_array[2,1]
            psf_final[idy+4,idx+4] += pcen*ipc_array[2,2]

    return psf_final

def real_psf(x,y,psf,psf_center=180):
    """
    Calculate the fraction of light from a psf that should 
    fall on a given pixel (x, y). This function assumes the 
    input ePSF was oversampled by a factor of 4.

    Parameters
    ----------
    x : float
        x location of source
    y : float
        y location of source
    psf : numpy.ndarray
        2D array with PSF image
    psf_center : int
        center of the input psf model

    Returns
    -------
    rpsf_phot : numpy.ndarray
        2D array with source placed
    """

    rx = psf_center + x*4
    ry = psf_center + y*4
    ix = int(rx)
    iy = int(ry)
    fx = rx-ix
    fy = ry-iy
    dd = np.sqrt(x**2+y**2)
    rpsf_phot = 0.
    if (np.abs(x) > 44) or (np.abs(y) > 44):
        return rpsf_phot
    if (dd > 4.0):
        rpsf_phot = ( (1-fx)*(1-fy)*psf[iy  ,ix  ]
                    + ( fx )*(1-fy)*psf[iy  ,ix+1]
                    + (1-fx)*( fy )*psf[iy+1,ix  ]
                    + ( fx )*( fy )*psf[iy+1,ix+1])
        return rpsf_phot

    # Bi-cubic interpolation
    A1 =  psf[iy  ,ix  ]
    B1 = (psf[iy  ,ix+1]-psf[iy  ,ix-1])/2
    C1 = (psf[iy+1,ix  ]-psf[iy-1,ix  ])/2
    D1 = (psf[iy  ,ix+1]+psf[iy  ,ix-1]-2*A1)/2
    F1 = (psf[iy+1,ix  ]+psf[iy-1,ix  ]-2*A1)/2
    E1 = (psf[iy+1,ix+1]-A1)
    A2 =  psf[iy  ,ix+1]
    B2 = (psf[iy  ,ix+2]-psf[iy  ,ix  ])/2
    C2 = (psf[iy+1,ix+1]-psf[iy-1,ix+1])/2
    D2 = (psf[iy  ,ix+2]+psf[iy  ,ix  ]-2*A2)/2
    F2 = (psf[iy+1,ix+1]+psf[iy-1,ix+1]-2*A2)/2
    E2 =-(psf[iy+1,ix  ]-A2)
    A3 =  psf[iy+1,ix  ]
    B3 = (psf[iy+1,ix+1]-psf[iy+1,ix-1])/2
    C3 = (psf[iy+2,ix  ]-psf[iy  ,ix  ])/2
    D3 = (psf[iy+1,ix+1]+psf[iy+1,ix-1]-2*A3)/2
    F3 = (psf[iy+2,ix  ]+psf[iy  ,ix  ]-2*A3)/2
    E3 =-(psf[iy  ,ix+1]-A3)
    A4 =  psf[iy+1,ix+1]
    B4 = (psf[iy+1,ix+2]-psf[iy+1,ix  ])/2
    C4 = (psf[iy+2,ix+1]-psf[iy  ,ix+1])/2
    D4 = (psf[iy+1,ix+2]+psf[iy+1,ix  ]-2*A4)/2
    F4 = (psf[iy+2,ix+1]+psf[iy  ,ix+1]-2*A4)/2
    E4 = (psf[iy  ,ix  ]-A4)
    V1 = (A1
        + B1*( fx )
        + C1*( fy )
        + D1*( fx )**2
        + E1*( fx )*( fy )
        + F1*( fy )**2)
    V2 = (A2
        + B2*(fx-1)
        + C2*( fy )
        + D2*(fx-1)**2
        + E2*(fx-1)*( fy )
        + F2*( fy )**2)
    V3 = (A3
        + B3*( fx )
        + C3*(fy-1)
        + D3*( fx )**2
        + E3*( fx )*(fy-1)
        + F3*(fy-1)**2)
    V4 = (A4
        + B4*(fx-1)
        + C4*(fy-1)
        + D4*(fx-1)**2
        + E4*(fx-1)*(fy-1)
        + F4*(fy-1)**2)
    rpsf_phot = ((1-fx)*(1-fy)*V1
               + ( fx )*(1-fy)*V2
               + (1-fx)*( fy )*V3
               + ( fx )*( fy )*V4)
    return rpsf_phot

def place_source(xpix, ypix, flux, image, psf, boxsize = 43):
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
    psf : numpy.ndarray
        2D array with ePSF

    Returns
    -------
    image : numpy.ndarray
        2D array of image where sources were placed
    """

    # Get image size, assuming square
    image_size = image.shape[0]

    # Apply if source is within image and boxise
    if (xpix > boxsize) & (xpix < image_size - boxsize) & (ypix > boxsize) & (ypix < image_size - boxsize):
        for j in range(round(ypix-boxsize),round(ypix+boxsize)):
            dy = j-ypix
            for i in range(round(xpix-boxsize),round(xpix+boxsize)):
                dx = i-xpix
                if (i < 1) or (i > image_size) or (j < 1) or (j > image_size):
                    pass
                else:
                    ff = real_psf(dx,dy,psf)
                    ffa=(flux*ff)
                    image[j,i]+=ffa

    return image

def get_psf(psf_file = 'working_psf.fits'):
    """
    Read in a PSF file
    """
    psf_data = fits.open(psf_file)
    return psf_data[0].data


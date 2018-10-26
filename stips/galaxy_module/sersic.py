# Copyright 2013 Justin Bird:
# http://physics.drexel.edu/~jbird/galflex/
#
# This file is part of GalFlex: The flexion simulation module
#
# GalFlex is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GalFlex is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GalFlex.  If not, see <http://www.gnu.org/licenses/>
#


"""
Sersic Profile creator.

This code has been somewhat modified from Justin Bird's GalFlex code by Brian York. The original
copyright and license notice are provided above.
"""

import numpy as np
import os
import sys

class baseImage(object):
    """Base class for defining the interface with which all GalFlex image classes access their
    shared methods and attributes.

        Images are composed of a 2-d, row-contiguous numpy array of NxN shape and a pixel scale.
    Representing sky coordinates, this scale is defined by a parameter `xmax`, and the image
    x,y axes run from -xmax to xmax, with the origin centered at (0,0).

    An Image can be constructed from either an existing NxN numpy array, or created from a
    number of predefined profiles (such as a Sersic):

        xmax = 20.0                 #Coordinates range from -20.0 to 20.0 arcsec on both axes
        array = numpy.ones((100,100))                   #Sample NxN array
        gal1 = baseImage(array, xmax)                   #Create instance

        xmax = 15.0
        n = 1.3                                             #Sersic index
        gal2 = Sersic(xmax,n,N=500,flux=1.e5,re=2.3)        #Create Sersic intensity profile


    The actual numpy array can be accessed by using the instance `.image` attribute.

    The individual elements in the `.image` attribute are indexed as im.image[y,x], matching the
    standard numpy convention of [row,column].

    @param image            An NxN numpy array
    @param xmax         Largest coordinate value along an array axis which runs from:
                        (-xmax) to (xmax).
    """

    instanceNum = 0
    instances = []
    
    def __init__(self,image,xmax):
        self.image = image
        self.N = len(image) #Image size
        self.xmax = xmax
        #Initialize method to run when script ends
        #Keep track of number of instances (for plot flagging)...
        baseImage.instanceNum+=1


    def __add__(self, instance):
        """Overload add operator for adding instance images together.

        @param instance     Instance of class baseImage.
        """
        try:
            m = self.image + instance.image
            return baseImage(m,self.xmax)
        except Exception:
            raise ValueError("Unable to add images. Make sure term is a baseImage instance")

    def __mul__(self, scale):
        """Overload mulitplication operator for scaling images by a scalar

        @param scale        Scalar that multiplies the entire image
        """
        try:
            self.image * float(scale)
        except Exception:
            raise ValueError("Unable to scale image")

    def copy(self):
        """Return a copy of instance
        """
        import copy
        baseImage.instanceNum+=1
        return copy.deepcopy(self)
        

    def getFlux(self):
        """Return sum of pixel values of image. 

        This may not be equal to the 'true' flux if the intensity profile is poorly sampled or
        is rapidly varying (i.e. high Sersic index)
        """
        return np.sum(self.image)
        
    def setFlux(self,new_flux):
        """Set flux of image to the specified value

        @param new_flux     Simple sum of photon counts in image.
        """
        self.image = self.image * float(new_flux)/np.sum(self.image)

    def addGaussianNoise(self,sigma=1.0,mean=0.0):
        """Add Gaussian noise to the image, centered at `mean` with standard deviation `sigma`

        Performs on a per-pixel basis.

        @param sigma        Standard deviation of distribution  [default = 1.0]
        @param mean         Mean of Gaussian distribution       [default = 0.0]
        """
        noise = np.random.RandomState(seed=self.seed).normal(loc=float(mean),scale=float(sigma),size=(self.N,self.N))
        self.image += noise
        return self.image

    def addPoissonNoise(self,sky_level=0.0):
        """Add Poisson noise to the image, based on the pixel fluxes.

        A sky_level can be specified per pixel, where the Poisson noise is derived from the
        sum of the image and the sky, and then the sky_level is subtracted from the image again.

        @param sky_level        Value per pixel of an (assumed) subtracted sky image [default = 0.0]
        """

        totalCount = self.image+float(sky_level)
        noise =  np.random.RandomState(seed=self.seed).poisson(totalCount)
        self.image = noise - float(sky_level)
        return self.image

    def writeFITS(self, filename, dir=None, overwrite=True):
        """Writes image to a FITS file using pyfits.

        This can be called directly as `galflex.writeFITS(image,...)` or as a baseImage method:
        `image.writeFITS(...)`. If called directly, accepts either a single instance or a list of
        instances as the first argument.

        Automatically creates a directory 'output' in the run directory where the FITS file is
        saved. If 'dir' is sprecified, that directory will replace 'output'

        @param filename     The name of the file to write to
        @param dir          Optionally a directory name can be provided if the filename does not 
                            already include it.     [Default = None]
        @param overwrite    Setting `overwrite=True` when `filename` is given will silently overwrite 
                            existing files.         [Default = True]
        """

        import pyfits                       #Keep pyfits optional
        import os

        hdu = pyfits.PrimaryHDU(self.image)
        hdulist = pyfits.HDUList([hdu])

        #If directory is specified, join filename to that directory.
        #   Otherwise create an 'output' directory                          
        if dir:
            filename = os.path.join(dir,filename)
        else:
            if not os.path.exists("output"):        # Make output directory
                os.makedirs("output")
            filename = 'output/'+str(filename)
    
        printname = str(filename)+'.fits'

        #Overwrite existing filename if overwrite=True
        if os.path.isfile(printname):
            if overwrite:
                os.remove(printname)
            else:
                raise IOError('File {} already exists'.format(printname))
        hdulist.writeto(printname, overwrite=True) # Save images as FITS file
        print("FITS file destination: {}".format(printname))


    def add(self,*args):
        """Another way for adding images together. Accepts baseImage class instances.

        Takes arguments in comma format, i.e. gal.add(img1,img2,img3...)

        @param *args    Either a baseImage class or a list of baseImage classes.
        """
        m = self.image
        for i in range(len(args)):
            try:
                m += args[i].image
            except Exception:
                print("Addition {} failed".format(i+1))
                pass
        return baseImage(m,self.xmax)   
            
class Sersic(baseImage):
    """Initialization class for a Sersic profile baseImage. Inherits all baseImage methods.

    Supports an elliptical Sersic intensity profile of form:

    I(r) ~ exp(-b(n)*(r/re)^(1/n)-1)

    or

    I(r) ~ exp(-(r/r0)^(1/n))

    depending on the input arguments. 

    User must specify either: 
        - The scale radius 'r0', where the intensity drops to 1/e the central value
        - The half-light-radius 're', which is the radius of the isophote containing
          half the total flux

    Uses an approximate expression for b(n):
        b(n) = 1.992*n - 0.3271
    which is a valid approximation in the index range 0.5<n<8.0.

    Some exponent values correspond to well-known models:
        n = 0.5: Gaussian profile
        n = 1.0: Exponential profile
        n = 4.0: deVaucouleurs profile

    @param xc   x centre of the image, where xc = N/2. would be a centred image.
    @param yc   y centre of the image, where yc = N/2. would be a centred image
    @param n    The exponent of intensity profile.
    @param xs   X Size of the image to produce (pixels)                             [default=500]
    @param ys   Y Size of the image to produce (pixels)                             [default=500]
    @param flux Total photon flux in image.                                         [default=1.0]
    @param q    Axis ratio of galaxy image (0.0 to 1.0), defined by x:(y*q).        [default=1.0]   
    @param phi  Orientation angle of galaxy image, in radians.                      [default=0.0]
    
    MUST specify one of either:
    
    @param re   Effective radius of the isophote containing half the total flux (half-light-radius).
    @param r0   Radius where the intensity drops to 1/e the central peak value.
    """
    def __init__(self,xc,yc,n,xs=500,ys=500,flux=1.0,q=1.0,phi=0.0, **kwargs):
        self.xc = float(xc)
        self.yc = float(yc)
        self.n = float(n)
        self.xs = int(xs)
        self.ys = int(ys)
        self.flux = float(flux)
        self.q = float(q)
        self.phi = float(phi)
        baseImage.instanceNum+=1
        
        max_offset = abs(np.array([self.xc,self.yc])).max()
        full_size = (max_offset*2.)**2 #full symmetrical image around centre
        image_size = self.xs * self.ys #fraction of image that's included
        self.flux = self.flux * image_size / full_size
        
#         assert xc >= 0., "ValueError: X centre must be greater than zero..."
#         assert yc >= 0., "ValueError: Y centre must be greater than zero..."
#         assert xc < xs, "ValueError: X centre must be in the image..."
#         assert yc < ys, "ValueError: Y centre must be in the image..."
        assert n > 0, "ValueError: Sersic Index must be larger than zero..."
        assert (q >= 0.0 and q<=1.0), "ValueError: Axis ratio q must be between 0.0 and 1.0"
        if flux<0.0:
            print("Warning: Negative total flux value entered")
    
        bn = 1.992*n - 0.3271 #Approximation valid for 0.5<n<8.0

        #Calculate both half-light-radius and scale radius
        try:
            r0 = kwargs['r0']           
            re = kwargs['r0'] * bn**n
        except Exception:
            try:
                re = kwargs['re']   #Half-light-radius
                r0 = re / bn**n     #Scale radius
            except Exception:
                raise ValueError("Invalid scale parameter")
        self.re = re
        self.r0 = r0
        
        #List of coordinates along an array axis
        coords_x = np.arange(xs)
        coords_y = np.arange(ys)
        X,Y = np.meshgrid(coords_x-xc,coords_y-yc) #Create 2 matrices: x values, y values of the coordinates
        x2 = np.cos(phi)*X+np.sin(phi)*Y
        y2 = -np.sin(phi)*X+np.cos(phi)*Y
        r = np.sqrt((x2**2+(y2*q)**2))
        img = np.exp(-bn*((r/re)**(1.0/n)-1.0)) #Calculate based on half-light-radius
        img = img.astype('float32')
        #Set total flux to flux
        preflux = img.sum()
        img = img*self.flux/preflux
        self.image = img

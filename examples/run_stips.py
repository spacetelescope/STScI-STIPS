import numpy as np
from astropy.io import fits
from astropy import table
from stips.observation_module import ObservationModule

def AB2Jy(mag):
    # Convert AB magnitudes to Jansky
    jy = 10**(0.4*(8.9-mag))
    return jy

def create_catalog(N,filename,ra,dec,fov,mags,band,overwrite=False,seed=None,output=None):
    """
    Method to create a catalog of N point sources distributed 
    uniformily over the detector (and beyond) around coordinates
    RA and DEC in degrees of size fov. The function will save an
    output .fits file that STIPS will be able to read.

    All sources assumed to be point sources with no rotation.

    Parameters
    ----------
    N : int
        Number of sources to generate.
    filename : str
        Name of output catalog fits file.
    ra : float
        RA of the field center, in degrees.
    dec : float
        DEC of the field center, in degrees.
    fov : float
        Field of view of the detected, in degrees.
    mags : tuple
        Range of magnitudes for sources.
    band : str
        Name of filter for sources, will be used by STIPS.
    overwrite : bool
        Overwrite existing .fits file?
    seed : int, default None
        Seed of random number generator.
    output : str, default None
        Return either a 'fits' table or a 'astropy' table.

    Returns
    -------
    hdut : astropy.io.fits.hdu.table.BinTableHDU
        fits table containing all the randomly generated sources
    """

    # Make a catalog of N sources
    cat = np.zeros(N,dtype=[('RA', np.float32),('Dec', np.float32),('mag', np.float32)])

    # Set seed if not None
    if isinstance(seed,(int,np.integer)):
        np.random.seed(seed)

    # Create uniform distribution of coordinates
    cs=np.cos(np.deg2rad(dec))
    lo=ra-fov/(2.0*cs)
    hi=ra+fov/(2.0*cs)
    cat['RA']=np.random.uniform(low=lo,high=hi,size=N)
    cat['Dec']=np.random.uniform(low=dec-fov/2,high=dec+fov/2,size=N)

    # Create random magnitudes
    cat['mag']=np.random.uniform(low=mags[0],high=mags[1],size=N)

    # Create the catalog in a format suitable for STIPS
    cols=[]
    cols.append(fits.Column(name='id'   ,array=np.arange(1,N+1,dtype=int),format='K' ))
    cols.append(fits.Column(name='ra'   ,array=cat['RA']                 ,format='E' ))
    cols.append(fits.Column(name='dec'  ,array=cat['Dec']                ,format='D' ))   
    cols.append(fits.Column(name='flux' ,array=AB2Jy(cat['mag'])         ,format='D' ))
    cols.append(fits.Column(name='type' ,array=['point']*N               ,format='8A'))
    cols.append(fits.Column(name='n'    ,array=np.empty(N)               ,format='D' ))
    cols.append(fits.Column(name='re'   ,array=np.empty(N)               ,format='D' ))
    cols.append(fits.Column(name='phi'  ,array=np.empty(N)               ,format='D' ))
    cols.append(fits.Column(name='ratio',array=np.empty(N)               ,format='D' ))
    cols.append(fits.Column(name='notes',array=['']*N                    ,format='8A'))
    cols.append(fits.Column(name='units',array=['j']*N                   ,format='8A'))
   
    # Create output fits table
    hdut = fits.BinTableHDU.from_columns(cols)

    # Add obligatory header keywords
    hdut.header['TYPE']='mixed' # specifies the catalog type for STIPS
    hdut.header['FILTER']=band  # specifies which filter this flux is in

    # Format filename
    if '.fits' not in filename:
        savename = filename + '.fits'
    elif filename.strip()[-5:] == '.fits':
        savename = filename.strip()
    else:
        raise ValueError(f'Format for {filename} not recognized.')

    # Write to disk
    hdut.writeto(savename,overwrite=overwrite)

    # Return catalog if requested
    if output:
        if output == 'fits':
            return hdut
        elif output == 'astropy':
            return table.Table(hdut.data)
        else:
            raise ValueError(f"param 'output' must be 'fits' or 'astropy', {output} not recognized.")

# Scene Parameters
filename = 'STIPS_cat.fits'
RA, DEC    = 90, 30 # center of the field
FILTER     = 'F129' # the name of the filter to use
EXPTIME    = 1000.0 # exposure time in sec
BACKGROUND = 1.5    # background level in e/s

# Create random sources catalog
create_catalog(500,filename,RA,DEC,fov=0.15,mags=(15,22),band=FILTER,overwrite=True,seed=42)

# Build observation parameters
obs = {'instrument'       : 'WFI', 
       'filters'          : [FILTER], 
       'detectors'        : 1,
       'distortion'       : False, 
       'oversample'       : 1,
       'pupil_mask'       : '', 
       'background'       : 'custom', 
       'custom_background': BACKGROUND, 
       'observations_id'  : 1, 
       'exptime'          : EXPTIME,
       'offsets'          : [{'offset_id'    : 1    , 
                              'offset_centre': False,
                              'offset_ra'    : 0.0  , 
                              'offset_dec'   : 0.0  ,
                              'offset_pa'    : 0.0  }]}
# Create observation object
obm = ObservationModule(obs,
                        ra    = RA,
                        dec   = DEC,
                        pa    = 0,
                        seed  = 42,
                        cores = 6)

# Initialize the local instrument
obm.nextObservation()
# Add catalog with sources
cat_name = obm.addCatalogue(filename)
# Add error to image
psf_file = obm.addError()
# Call the final method
fits_file, mosaic_file, params = obm.finalize(mosaic=False)
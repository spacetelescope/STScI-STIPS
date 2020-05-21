from __future__ import absolute_import, division

# External modules
import logging, os, sys

import numpy as np

from astropy.table import Table, Column
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from cStringIO import StringIO

# Local modules
from ..stellar_module import StarGenerator
from ..utilities import GetStipsData, OffsetPosition, StipsDataTable
from .convert_units import (DivideInterval, RadiiUnknown2Arcsec, RadiiUnknown2Parsec, RescaleArray)

#-----------
class SceneModule(object):

    #-----------
    def __init__(self, **kwargs):
        """
        Noiseless scene generator module.

        :Author: Brian York

        :Organization: Space Telescope Science Institute

        :History:
            * 2010/10/19 PLL created this module.
            * 2011/06/14 PLL added single star simulation.
            * 2011/06/28 PLL reorganized functions.
            * 2011/10/28 PLL added galaxies simulation.
            * 2014/02/14 BY modified the code to be instrument-independent

        Examples
        --------
        >>> from stips import SceneModule


        Parameters
        ----------
        self: obj
            Class instance.

        **kwargs: dictionary
            Additional arguments needed to make the scene

        """
        self.out_path = kwargs.get('out_path', os.getcwd())
        self.prefix = kwargs.get('out_prefix', 'sim')
        self.cat_type = kwargs.get('cat_type', 'fits')
        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            self.logger = logging.getLogger('__stips__')
            self.logger.setLevel(logging.INFO)
            if not len(self.logger.handlers):
                stream_handler = logging.StreamHandler(sys.stderr)
                stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))# [in %(pathname)s:%(lineno)d]'))
                self.logger.addHandler(stream_handler)
        if 'scene_general' in kwargs:
            self.ra = kwargs['scene_general'].get('ra', 0.0)
            self.dec = kwargs['scene_general'].get('dec', 0.0)
            self.seed = kwargs['scene_general'].get('seed', 0)
        else:
            self.ra = kwargs.get('ra', 0.0)
            self.dec = kwargs.get('dec', 0.0)
            self.seed = kwargs.get('seed', 0)

        self.params =  [ 'Random seed: %d' % (self.seed) ]
        self.params += [ 'Centre (RA,DEC) = (%f,%f)' % (self.ra,self.dec) ]
        self.catalogues = {}
    
    #-----------
    def CreatePopulation(self, pop, id=0):
        """
        Generate a stellar population.

        Output list will have these columns:
            # ID
            # RA
            # DEC
            # Distance
            # Age
            # Metallicity
            # Mass
            # Count Rate (in the chosen instrument/filter), Absolute
            # Count Rate (in the chosen instrument/filter), Observed

        Parameters
        ----------
        self: obj
            Class instance.

        pop: dictionary
            Information about the population. Includes:
                n_stars: int
                    Number of stars
                age_low,age_high: floating point
                    Minimum and maximum ages.
                z_low,z_high: floating point
                    Minimum and maximum metallicities.
                imf: string
                    Initial Mass Function for the population
                alpha: float
                    Exponent (if imf = 'powerlaw')
                binary_fraction: float
                    Binary Fraction
                distribution: string
                    Stellar distribution in the sky (e.g. power law, inverse power law, uniform, etc.)
                clustered: bool
                    Cluster higher masses to centre?
                radius: float
                    Radius in (units)
                radius_units: string
                    Units of radius (above)
                distance_low,distance_high: floating point
                    Minimum and maximum distance (in kpc) of the 'cluster'
                offset_ra,offset_dec: floating point
                    Offset of the cluster from the scene centre in arcseconds
        
        Returns
        -------
        
        outList: string
            The catalogue file produced
        """
        star_chunk = 100000
        age_bins = DivideInterval("1.0e6,1.35e10,d5")
        met_bins = DivideInterval("-2.5,0.5,i0.1")
        out_file = "{}_stars_{:03d}.{}".format(self.prefix, id, self.cat_type)
        outList = os.path.join(self.out_path, out_file)
        if os.path.isfile(outList): 
            os.remove(outList) # No append
        data_table = StipsDataTable.dataTableFromFile(outList)

        self._log("info","Creating catalogue %s" % (outList))
        n_stars = int(pop['n_stars'])
        age_l = age_bins[(np.abs(age_bins-float(pop['age_low']))).argmin()]
        age_h = age_bins[(np.abs(age_bins-float(pop['age_high']))).argmin()]
        age_bins = age_bins[np.where((age_bins>=age_l) & (age_bins<=age_h))]
        met_l = met_bins[(np.abs(met_bins-float(pop['z_low']))).argmin()]
        met_h = met_bins[(np.abs(met_bins-float(pop['z_high']))).argmin()]
        met_bins = met_bins[np.where((met_bins>=met_l) & (met_bins<=met_h))]
        imf = pop['imf']
        alpha = abs(float(pop['alpha']))
        distribution = pop['distribution']
        clustered = pop['clustered']
        radius = float(pop['radius'])
        rad_units = pop['radius_units']
        dist_l = float(pop['distance_low']) * 1.e3 #convert kpc to pc
        dist_h = float(pop['distance_high']) * 1.e3 #convert kpc to pc
        binary_fraction = float(pop['binary_fraction'])
        offset_ra = float(pop['offset_ra'])/3600. #offset in RA arcseconds, convert to degrees.
        offset_dec = float(pop['offset_dec'])/3600. #offset in DEC arcseconds, convert to degrees.
        
        metadata = {'type': 'phoenix', 'id': id, 'n_stars': n_stars, 'age_l': age_l, 'age_h': age_h,
                    'met_l': met_l, 'met_h': met_h, 'imf': imf, 'alpha': alpha,
                    'distribution': distribution, 'clustered': clustered, 'radius': radius,
                    'radius_units': rad_units, 'dist_l': dist_l, 'dist_h': dist_h, 
                    'offset_ra': offset_ra, 'offset_dec': offset_dec, 
                    'name': 'Phoenix Stellar Population Table'}
        data_table.meta = metadata

        self._log("info","Creating age and metallicity numbers")
#         ages = np.random.RandomState(seed=self.seed).random_sample(size=len(age_bins))
        ages = np.random.random_sample(size=len(age_bins))
        ages /= ages.sum()
#         mets = np.random.RandomState(seed=self.seed).random_sample(size=len(met_bins))
        mets = np.random.random_sample(size=len(met_bins))
        mets /= mets.sum()
        self._log("info","Created age and metallicity numbers")
        
        self._log("info","Creating stars")
        #Generate star masses
        datasets = 0
        total = 0
        for i, age in enumerate(age_bins):
            self.logger.info("Age %g",age)
            n_age = int(round(n_stars * ages[i]))
            for j, met in enumerate(met_bins):
                self.logger.info("Metallicity %f",met)
                num_stars = int(round(n_age * mets[j]))
                if num_stars == 0:
                    continue
                self.logger.info("Creating %d stars",num_stars)
                stargen = StarGenerator(age, met, imf=imf, alpha=alpha, seed=self.seed, logger=self.logger)
                all_masses, all_rates, all_temps, all_gravs = stargen.make_cluster(num_stars)
                all_x, all_y, all_z = self._MakeCoords(num_stars, radius, func=distribution, scale=2.8, do_z=True)
#                 all_distances = np.random.RandomState(seed=self.seed).uniform(low=dist_l, high=dist_h, size=num_stars)
                all_distances = np.random.uniform(low=dist_l, high=dist_h, size=num_stars)
                if clustered:
                    all_x, all_y, all_z = self._CenterObjByMass(all_x, all_y, all_masses, z=all_z)
#                 all_binaries = np.random.RandomState(seed=self.seed).binomial(1,binary_fraction,len(all_masses))
                all_binaries = np.random.binomial(1,binary_fraction,len(all_masses))
                idx = np.where(all_binaries==1)[0]
                mb, rb, tb, gb = stargen.make_cluster(len(idx))
                xb, yb, zb = all_x[idx], all_y[idx], all_z[idx]
                db = all_distances[idx]
                all_masses = np.insert(all_masses, idx, mb)
                all_rates = np.insert(all_rates, idx, rb)
                all_temps = np.insert(all_temps, idx, tb)
                all_gravs = np.insert(all_gravs, idx, gb)
                all_x = np.insert(all_x, idx, xb)
                all_y = np.insert(all_y, idx, yb)
                all_z = np.insert(all_z, idx, zb)
                all_distances = np.insert(all_distances, idx, db)
                all_binaries = np.insert(all_binaries, idx+1, 0)
                num_stars += len(idx)
                cached_ra = 0.
                cached_dec = 0.
                cached_distance = 0.
                cached = False
                for k in range(num_stars // star_chunk + 1):
                    xl, xh = k * star_chunk, min(k * star_chunk + star_chunk, num_stars-1)
                    star_set = xh - xl
                    self._log("info", "Chunk {}: {} stars".format(k+1, star_set))
                    masses = all_masses[xl:xh]
                    rates = all_rates[xl:xh]
                    temps = all_temps[xl:xh]
                    gravs = all_gravs[xl:xh]
                    x, y, z = all_x[xl:xh], all_y[xl:xh], all_z[xl:xh]
                    distances = all_distances[xl:xh]
                    binaries = all_binaries[xl:xh]
                    ids = np.arange(total + xl, total + xh) + 1
                    x = RadiiUnknown2Arcsec(x, rad_units, distances)
                    y = RadiiUnknown2Arcsec(y ,rad_units, distances)
                    z = RadiiUnknown2Parsec(z, rad_units, distances)
                    distances += z
                    ras = x/3600. #decimal degrees
                    decs = y/3600. #decimal degrees
                    base_ra,base_dec = OffsetPosition(self.ra,self.dec,offset_ra,offset_dec)
                    decs += base_dec
                    idxg = np.where(decs>90.)
                    idxl = np.where(decs<-90.)
                    decs[idxg] = 180. - decs[idxg]
                    ras[idxg] = 180. + ras[idxg]
                    decs[idxl] = -180. - decs[idxl]
                    ras[idxl] = 180. + ras[idxl]
                    ras = (ras + base_ra)%360
                    apparent_rates = rates + (5.0 * np.log10(distances) - 5.0)

                    t = Table()
                    t['id'] = Column(data=ids, format="%8d")
                    t['ra'] = Column(data=ras, unit='degrees', format="%17.9e")
                    t['dec'] = Column(data=decs, unit='degrees', format="%17.9e")
                    t['distance'] = Column(data=distances, unit='pc', format="%17.9e")
                    t['age'] = Column(data=np.full_like(ids, age), unit='years', format="%17d")
                    t['metallicity'] = Column(data=np.full_like(ras, met), format="%4.1f")
                    t['mass'] = Column(data=masses,unit='solar masses', format="%17.9e")
                    t['teff'] = Column(data=temps, unit='K', format="%13.8f")
                    t['log_g'] = Column(data=gravs, format="%12.9f")
                    t['binary'] = Column(data=binaries, format="%3d")
                    t['dataset'] = Column(data=np.full_like(ids, datasets), format="%6d")
                    t['absolute'] = Column(data=rates,unit='johnson,i', format="%14.6e")
                    t['apparent'] = Column(data=apparent_rates,unit='johnson,i', format="%12.4e")
                    data_table.write_chunk(t)
                    del t
                datasets += 1
                total += num_stars

        self._log("info","Done creating catalogue")
        return outList

    #-----------
    def CreateGalaxies(self, gals, id=0):
        """
        Generate galaxies list.

        Output list will have these columns:
            # ID
            # RA
            # DEC
            # Redshift
            # Model
            # Age
            # Profile
            # Half-flux_radius
            # Axial_ratio
            # Position_angle
            # Johnson,V absolute
            # Johnson,V apparent

        Parameters
        ----------
        self: obj
            Class instance.

        gals: dictionary
            Information about the galaxies. Includes:
                n_gals: int
                    Number of galaxies
                z_low,z_high: float
                    Minimum and maximum redshifts (converted to distances?).
                rad_low,rad_high: float
                    Minimum and maximum galactic half-light radii (in arcseconds)
                sb_v_low, sb_v_high: float
                    Minimum and maximum V-band average surface brightness within rad
                distribution: string
                    Stellar distribution in the sky (e.g. power law, inverse power law, uniform, etc.)
                clustered: bool
                    Cluster higher masses to centre?
                radius: float
                    Radius in (units)
                radius_units: string
                    Units of radius (above)
                offset_ra,offset_dec: float
                    Offset of cluster from scene centre in mas

        Returns
        -------
        
        outList: string
            The catalogue file produced
        """
        bc95_models = np.array(('a','b','c','d','e','f'))
        bc95_ages = np.array(("10E5","25E5","50E5","76E5","10E6","25E6","50E6","10E7","50E7","10E8","50E8","10E9"))
        out_file = "{}_gals_{:03d}.{}".format(self.prefix, id, self.cat_type)
        outList = os.path.join(self.out_path, out_file)
        if os.path.isfile(outList): 
            os.remove(outList) # No append
        data_table = StipsDataTable.dataTableFromFile(outList)
        
        # Write star list (overwrite)
        self.logger.info("Creating catalogue %s",outList)
        # Generate galaxy list
        n_gals = int(gals['n_gals'])
        z_l = float(gals['z_low'])
        z_h = float(gals['z_high'])
        r_l = float(gals['rad_low'])
        r_h = float(gals['rad_high'])
        m_l = float(gals['sb_v_low'])
        m_h = float(gals['sb_v_high'])
        distribution = gals['distribution']
        clustered = gals['clustered']
        radius = float(gals['radius'])
        rad_units = gals['radius_units']
        offset_ra = float(gals['offset_ra'])/3600. #offset in RA arcseconds, convert to degrees.
        offset_dec = float(gals['offset_dec'])/3600. #offset in DEC arcseconds, convert to degrees.
        self._log("info","Wrote preamble")
        self._log("info","Parameters are: {}".format(gals))

        ids = np.arange(n_gals)

        # Roughly 50% spiral, 50% elliptical
        ellipRatio = 0.5
#         binoDist = np.random.RandomState(seed=self.seed).binomial(1, ellipRatio, n_gals)
        binoDist = np.random.binomial(1, ellipRatio, n_gals)
        idx_ellip = np.where(binoDist == 1)
        idx_spiral = np.where(binoDist != 1)
        types = np.array( ['expdisk'] * n_gals )
        types[idx_ellip] = 'devauc'
        n_ellip = len( idx_ellip[0] )
        n_spiral = n_gals - n_ellip

        # Axial ratio
        # Spiral     = 0.1 to 1
        # Elliptical = 0.5 to 1
        axialRatioSpiralMin, axialRatioSpiralMax = 0.1, 1.0
        axialRatioEllipMin,  axialRatioEllipMax  = 0.5, 1.0
        axials = np.zeros(n_gals)
#         axials[idx_spiral] = np.random.RandomState(seed=self.seed).uniform(axialRatioSpiralMin, axialRatioSpiralMax, n_spiral)
#         axials[idx_ellip] = np.random.RandomState(seed=self.seed).uniform(axialRatioEllipMin,  axialRatioEllipMax, n_ellip)
        axials[idx_spiral] = np.random.uniform(axialRatioSpiralMin, axialRatioSpiralMax, n_spiral)
        axials[idx_ellip] = np.random.uniform(axialRatioEllipMin,  axialRatioEllipMax, n_ellip)
        
        # Position angle
        posAngleAlgo = 'uniform'
#         angles = np.random.RandomState(seed=self.seed).uniform(0.0, 359.9, n_gals)
        angles = np.random.uniform(0.0, 359.9, n_gals)

        # Half-flux radius - uniform
#         rads = np.random.RandomState(seed=self.seed).uniform(r_l, r_h, n_gals)
        rads = np.random.uniform(r_l, r_h, n_gals)
        
        # Redshifts
        # If both z_low and z_high are zero, do local galaxies. Distance is 0.5 Mpc -- 50 Mpc.
        # In the future, offer an option for straight distance or redshift.
        if z_l == 0. and z_h == 0.:
            z_label = "distance"
#             distances = np.random.RandomState(seed=self.seed).uniform(5.e5, 5.e7, n_gals)
            distances = np.random.uniform(5.e5, 5.e7, n_gals)
            zs = distances / 1.e3
            convs = np.log10(distances)
        else:
            z_label = "redshift"
#             zs = np.random.RandomState(seed=self.seed).uniform(z_l, z_h, n_gals)
            zs = np.random.uniform(z_l, z_h, n_gals)
            distances = np.array(cosmo.comoving_distance(zs).to(u.pc))
            convs = np.log10(np.array(cosmo.luminosity_distance(zs).to(u.pc)))

        # Luminosity function - power law
        lumPow = -1.8
#         vmags = np.random.RandomState(seed=self.seed).power(np.abs(lumPow)+1.0, size=n_gals)
        vmags = np.random.power(np.abs(lumPow)+1.0, size=n_gals)
        if lumPow < 0: vmags = 1.0 - vmags
        vmags = RescaleArray(vmags, m_l, m_h)
        vmags_abs = vmags - 5*(convs-1.)
        
#         models = np.random.RandomState(seed=self.seed).choice(bc95_models,size=n_gals)
#         ages = np.random.RandomState(seed=self.seed).choice(bc95_ages,size=n_gals)
        models = np.random.choice(bc95_models,size=n_gals)
        ages = np.random.choice(bc95_ages,size=n_gals)

        self._log("info","Making Co-ordinates")
        x,y = self._MakeCoords(n_gals,radius,func=distribution,scale=2.8)
        x = RadiiUnknown2Arcsec(x,rad_units,distances)
        y = RadiiUnknown2Arcsec(y,rad_units,distances)
        if clustered:
            self._log("info","Clustering")
            x,y = self._CenterObjByMass(x,y,1/vmags)
        self._log("info","Converting Co-ordinates into RA,DEC")
        ras = x/3600. #decimal degrees
        decs = y/3600. #decimal degrees
        base_ra,base_dec = OffsetPosition(self.ra,self.dec,offset_ra,offset_dec)
        decs += base_dec
        idxg = np.where(decs>90.)
        idxl = np.where(decs<-90.)
        decs[idxg] = 180. - decs[idxg]
        ras[idxg] = 180. + ras[idxg]
        decs[idxl] = -180. - decs[idxl]
        ras[idxl] = 180. + ras[idxl]
        ras = (ras + base_ra)%360
        
        metadata = {'type': 'bc95', 'id': id, 'n_gals': n_gals, 'z_l': z_l, 'z_h': z_h, 
                    'radius_l': r_l, 'radius_h': r_h, 'sb_v_l': m_l, 'sb_v_h': m_h,
                    'distribution': distribution, 'clustered': clustered, 'radius': radius,
                    'radius_units': rad_units, 'offset_ra': offset_ra, 'offset_dec': offset_dec,
                    'name': 'Galaxy Population Table'}
        data_table.meta = metadata
        t = Table()
        t['id'] = Column(data=ids)
        t['ra'] = Column(data=ras,unit='degrees')
        t['dec'] = Column(data=decs,unit='degrees')
        t[z_label] = Column(data=zs)
        t['model'] = Column(data=models)
        t['age'] = Column(data=ages,unit='years')
        t['profile'] = Column(data=types)
        t['radius'] = Column(data=rads)
        t['axial_ratio'] = Column(data=axials,unit='degrees')
        t['pa'] = Column(data=angles,unit='degrees')
        t['absolute_surface_brightness'] = Column(data=vmags_abs,unit='johnson,v')
        t['apparent_surface_brightness'] = Column(data=vmags,unit='johnson,v')
        data_table.write_chunk(t)
        self._log("info","Done creating catalogue")
        return outList

    #-----------
    def _CenterObjByMass(self, x, y, mass, z=None):
        """
        Place slightly more massive stars near image center
        to simulate mass segragation.

        Parameters
        ----------
        self: obj
            Class instance.

        x, y: array_like
            Initial coordinates of object placement.

        mass: array_like
            Stellar masses.
        
        z: array_like, optional
            Initial z co-ordinates of object placement. If provided, return 3D co-ordinates.
        
        Returns
        -------
        new_x, new_y, [new_z]: array_like
            Re-ordered `x` and `y`.

        """
        x_cen = 0.
        y_cen = 0.
        z_cen = 0.
        n_stars = x.size

        # Central coordinates will have smallest values
        dx = x - x_cen
        dy = y - y_cen
        if z is not None:
            dz = z - z_cen
            d = np.sqrt(dy*dy + dx*dx + dz*dz)
        else:
            d = np.sqrt(dy*dy + dx*dx)
        i_sorted = np.argsort(d)

        # Segregate mass
        m_cut = mass.max() * 0.75 # Arbitrary cut of low/high masses
        j_lo = np.where( mass <= m_cut )[0]
        j_hi = np.where( mass >  m_cut )[0]

        # Place high masses with slightly more preference for center
        j_sorted = np.zeros(n_stars, dtype='int')
        k_lo, k_hi = 0, 0
        x_fac = n_stars // j_hi.size
        for i in range(n_stars):
            if i%x_fac == 0 and k_hi < j_hi.size:
                j_sorted[i] = j_hi[k_hi]
                k_hi += 1
            else:
                j_sorted[i] = j_lo[k_lo]
                k_lo += 1

        # Match sorted coordinates to masses
        new_x = np.zeros(n_stars)
        new_y = x.copy()
        if z is not None:
            new_z = x.copy()
        for i in range(n_stars):
            new_x[ j_sorted[i] ] = x[ i_sorted[i] ]
            new_y[ j_sorted[i] ] = y[ i_sorted[i] ]
            if z is not None:
                new_z[ j_sorted[i] ] = z[ i_sorted[i] ]
        
        if z is not None:
            return new_x, new_y, new_z
        else:
            return new_x, new_y

    #-----------
    def _MakeCoords(self, numObj, radMax, func='exp', scale=1.0, do_z=False):
        """
        Position objects at image center according to given
        density profile.

        Parameters
        ----------
        self: obj
            Class instance.

        numObj: int
            Number of objects to generate.

        radMax: float
            Radius of object collection.

        func: {'exp', 'invpow', 'regpow', 'uniform'}
            Exponential, inverse power law, regular power
            law, or uniform.

        scale: float
            Scale parameter for random generator.
            Definition depends on `func`. Not used for
            uniform distribution.
        
        do_z: bool
            If true, provide 3D co-ordinates. Else, provide 2D co-ordinates.

        Returns
        -------
        x, y[, z]: array_like
            Coordinates for `numObj` objects.
            
        """
        self.logger.info("Creating {} objects, max radius {}, function {}, scale {}".format(numObj, radMax, func, scale))
        x = np.array([])
        y = x.copy()
        if do_z:
            z = x.copy()

        isPolar = True

        # Random radii for polar coordinates, except uniform
        if func == 'exp':
#             r_arr = np.random.RandomState(seed=self.seed).exponential(scale=scale, size=numObj)
            r_arr = np.random.exponential(scale=scale, size=numObj)
        elif func == 'invpow':
#             r_arr = 1.0 - np.random.RandomState(seed=self.seed).power(scale, size=numObj)
            r_arr = 1.0 - np.random.power(scale, size=numObj)
        elif func == 'regpow':
#             r_arr = np.random.RandomState(seed=self.seed).power(scale, size=numObj)
            r_arr = np.random.power(scale, size=numObj)
        elif func == 'uniform':
            isPolar = False
#             x = np.random.RandomState(seed=self.seed).uniform(-radMax, radMax, numObj)
#             y = np.random.RandomState(seed=self.seed).uniform(-radMax, radMax, numObj)
            x = np.random.uniform(-radMax, radMax, numObj)
            y = np.random.uniform(-radMax, radMax, numObj)
            if do_z:
#                 z = np.random.RandomState(seed=self.seed).uniform(-radMax, radMax, numObj)
                z = np.random.uniform(-radMax, radMax, numObj)
        else:
            raise ValueError('Invalid _MakeCoords func')

        # For polar coordinates only
        if isPolar:
            r_arr = RescaleArray(r_arr, 0, radMax)

            if do_z:
                # Random angles for spherical co-ordinates
#                 t_arr = np.random.RandomState(seed=self.seed).uniform(0,np.pi,numObj)
#                 p_arr = np.random.RandomState(seed=self.seed).uniform(0,np.pi*2,numObj)
                t_arr = np.random.uniform(0,np.pi,numObj)
                p_arr = np.random.uniform(0,np.pi*2,numObj)
                
                x = r_arr * np.sin(t_arr) * np.cos(p_arr)
                y = r_arr * np.sin(t_arr) * np.sin(p_arr)
                z = r_arr * np.cos(t_arr)
            else:
                # Random angles for polar coordinates
#                 t_arr = np.random.RandomState(seed=self.seed).uniform(0, np.pi*2, numObj)
                t_arr = np.random.uniform(0, np.pi*2, numObj)
                
                # Convert to cartesian coordinates
                x = r_arr * np.cos(t_arr)
                y = r_arr * np.sin(t_arr)

        if do_z:
            return x,y,z
        else:
            return x, y
    
    def _log(self,mtype,message):
        """
        Checks if a logger exists. Else prints.
        """
        if hasattr(self,'logger'):
            getattr(self.logger,mtype)(message)
        else:
            sys.stderr.write("%s: %s\n" % (mtype,message))

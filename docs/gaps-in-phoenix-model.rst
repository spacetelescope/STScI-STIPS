
Phoenix model grid
===================

Grids for the pre-computed Phoenix models are available for each Roman filter. These grids will be used when the input table type is set to 'type=phoenix' for point sources.  The grids cover stellar parameters for a range of effective temperature (teff), surface gravity (logg), and metallicity (Z). 
    * Effective temperature spans from 2000K to 70000K 
        * at 100K increment from 2000K to 7000K, 
        * at 200K increment from 7000K to 12000K, 
    	* at 500K increment from 12000K to 20000K, and 
    	* at 1000K increment from 20000K to 70000K. 
    *	Surface gravity spans from 0.0 to 5.5 at 0.5 increment
    *	Metallicity spans from -4.0 to 0.5 dex 
    	* at 0.5 dex increment from -4.0 to 0.0 dex
    	* above 0.0 dex, two points are available at 0.3 dex and 0.5 dex
 
However, there are gaps in the model where the Phoenix models are not available, due to the set of grid point being non-physical. The gaps exist at
    *	At logg 0.0 
    	* With Z from -4.0 dex to -2.0 dex and teff above 3000K
    	* With Z from -1.5 dex to 0.5 dex and teff above 12500K
    *	At logg 0.5 
    	* With Z from -4.0 dex to -2.0 dex and teff above 3000K
    	* With Z from -1.5 dex to 0.5 dex and teff above 15500K
    *	At logg 1.0, teff above 21000K for all Z
    *	At logg 1.5, teff above 26000K for all Z
    *	At logg 2.0, teff above 31000K for all Z
    *	At logg 2.5, 
    	* Teff above 41000K for Z's between -4.0 dex to -0.5dex, 0.3 dex, and 0.5 dex
    	* Teff 41000K, 43000K, 45000K, 47000K, 49000K, and above 51000K for 0.0 dex 
    *	At logg 3.0, teff above 41000K for all Z
    *	At logg 3.5, teff above 51000K for all Z
    *	At logg 5.0, 
    	* Teff above 10200K for Z’s between -4.0 dex to -1.5 dex 
    	* Teff above 10000K for Z‘s between -1.0 dex and 0.3 dex
    	* Teff above 3000K for Z of 0.5 dex
    *	At logg 5.5
    	* Teff above 5100K for Z’s between -4.0dex to 0.3dex
    	* Teff above 3000K for Z of 0.5 dex


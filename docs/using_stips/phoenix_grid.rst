Phoenix Model Grid and Gaps
===========================

Grids for the pre-computed Phoenix models are available for each Roman filter. These grids will be used when the input table type is set to :ref:`'type=phoenix' <phoenix-tables>` for point sources. The grids cover stellar parameters for a range of effective temperature (teff), surface gravity (logg), and metallicity (Z).

* Effective temperature spans from 2,000 K to 70,000 K:

    * at 100 K increment from 2,000 K to 7,000 K

    * at 200 K increment from 7,000 K to 12,000 K

    * at 500 K increment from 12,000 K to 20,000 K

    * at 1000 K increment from 20,000 K to 70,000 K

* Surface gravity spans from 0.0 to 5.5 at 0.5 increment

* Metallicity spans from -4.0 to 0.5 dex

    * at 0.5 dex increment from -4.0 to 0.0 dex

    * above 0.0 dex, two points are available at 0.3 dex and 0.5 dex

Gaps exist where the Phoenix models are not available due to the combination of grid points being non-physical. These are the locations of the gaps:

*   At logg 0.0:

    * With teff above 3,000 K and Z from -4.0 dex to -2.0 dex

    * With teff above 12,500 K and Z from -1.5 dex to 0.5 dex and

*   At logg 0.5:

    * With teff above 3,000 K and Z from -4.0 dex to -2.0 dex

    * With teff above 15,500 K and Z from -1.5 dex to 0.5 dex

*   At logg 1.0:

    * With teff above 21,000K and all Z

*   At logg 1.5:

    * With teff above 26,000K and all Z

*   At logg 2.0:

    * With teff above 31,000K and all Z

*   At logg 2.5:

    * With teff above 41,000 K and Z between -4.0 dex to -0.5 dex, at 0.3 dex, or at 0.5 dex

    * With teff at 41,000 K, 43,000 K, 45,000 K, 47,000 K, 49,000 K, or above 51,000 K and Z at 0.0 dex

*   At logg 3.0:

    * With teff above 41,000 K and all Z

*   At logg 3.5:

    * With teff above 51,000 K and all Z

*   At logg 5.0:

    * With teff above 10,200 K and Z between -4.0 dex to -1.5 dex

    * With teff above 10,000 K and Z between -1.0 dex to 0.3 dex

    * With teff above 3,000 K and Z at 0.5 dex

*   At logg 5.5:

    * With Teff above 5,100 K and Z between -4.0 dex to 0.3 dex

    * With Teff above 3,000 K and Z at 0.5 dex

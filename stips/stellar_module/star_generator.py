"""
Various codes to work with the initial mass function

Written by Adam Ginsburg. Distributed under the MIT license (below):

--------------------------------------------------------------------------------
The MIT License (MIT)

Copyright (c) 2013 Adam Ginsburg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
--------------------------------------------------------------------------------

Modified by Brian York starting November 14, 2013.
"""
import os
import numpy as np
import sqlite3

from scipy.interpolate import RegularGridInterpolator

from esutil import sqlite_util

from ..utilities import GetStipsData

class StarGenerator(object):
    def __init__(self, age, metallicity, **kwargs):
        """
        Stellar population generator, given an IMF, age, and metallicity.
        
        :Author: Brian York. IMF code written by Adam Ginsburg, and released under the MIT license.

        :Organization: Space Telescope Science Institute

        :History:
            * 2013/11/14 BY created this module.

        Examples
        --------

        Parameters
        ----------
        self: obj
            Class instance.

        imf: string
            IMF function to use. Currently one of Salpeter, Kroupa, Chabrier, Schechter, Modified Schechter

        age: float
            Population age. Will be rounded to the nearest match in the database.
            Minimum: 1.0e6
            Maximum: 13.5e9
        
        metallicity: float
            Population metallicity. Will be rounded to the nearest match in the database.
            Minimum: -2.5
            Maximum: +0.5
        *****
        So, here's what it's going to be:
            - give it an age and a metallicity and an IMF
                - it'll figure out minimum and maximum masses based on the grid
            - give it a number of stars
                - it'll give back masses and absolute magnitude count rates
            - At least, that's the plan
        """
        massfunctions = {   
                            'powerlaw':self.powerlaw,
                            'kroupa':self.kroupa, 
                            'salpeter':self.salpeter, 
                            'chabrier':self.chabrier, 
                            'schechter':self.schechter,
                            'modified_schechter':self.modified_schechter
                        }
        self.metallicity = metallicity
        self.age = age
        self.dbname = GetStipsData(kwargs.get('dbname', 'IsochroneGrid.db'))
        self.gridpath = GetStipsData(kwargs.get('gridpath', 'grid'))
        self.imf = massfunctions[kwargs.get('imf', 'powerlaw').lower().replace(" ","_")]
        self.alpha = kwargs.get('alpha', -2.35)
        self.logger = kwargs.get('logger', None)
        self.seed = kwargs.get('seed', 1234)
        db = sqlite3.connect(self.dbname)
        c = db.cursor()
        c.execute("""SELECT age FROM star_info ORDER BY ABS(? - age) LIMIT 1""",(self.age,))
        res = c.fetchall()
        self.age = res[0][0]
        c.execute("""SELECT Z FROM star_info ORDER BY ABS(? - Z) LIMIT 1""",(self.metallicity,))
        res = c.fetchall()
        self.metallicity = res[0][0]
        c.execute("""SELECT MAX(m_ini) FROM star_info WHERE age = ? and Z = ?""",(self.age,self.metallicity))
        res = c.fetchall()
        self.mmax = res[0][0]
        c.execute("""SELECT MIN(m_ini) FROM star_info WHERE age = ? and Z = ?""",(self.age,self.metallicity))
        res = c.fetchall()
        self.mmin = res[0][0]
        db.close()
    
    # four codes for dn/dlog(m)
    def powerlaw(self,m,integral=False,alpha=-2.35,**kwargs):
        """
        Note that Salpeter (1955) is a power law curve with alpha=2.35
        """
        if integral: alpha -= 1
        return m**(-alpha)

    def salpeter(self,m,alpha=2.35, integral=False,**kwargs):
        """
        the Salpeter 1955 IMF: dn/dm ~ m^-2.35
        """
        return self.powerlaw(m,alpha=2.35,integral=integral)

    def kroupa(self,m, integral=False,**kwargs):
        """
        Kroupa 2001 IMF (http://arxiv.org/abs/astro-ph/0009005, http://adsabs.harvard.edu/abs/2001MNRAS.322..231K)
        """
        exp1 = 0.3
        exp2 = 1.3
        exp3 = 2.3
        if integral: 
            exp1 += 1
            exp2 -= 1
            exp3 -= 1
        zeta = (m**exp1 / 0.08**exp1 * 0.08**-exp2)*(m<0.08)
        zeta += (m**-exp2) * (m>=0.08) * (m<0.5)
        zeta += (m**-exp3 / 0.5**-exp3 * 0.5**-exp2) * (m>=0.5)
        return zeta

    def chabrier(self,m, integral=False,**kwargs):
        """
        Chabrier 2003 IMF
        http://adsabs.harvard.edu/abs/2003PASP..115..763C
        (only valid for m < 1 msun)

        not sure which of these to use...

        integral is NOT IMPLEMENTED
        """
        if integral: log("error", "Chabrier integral NOT IMPLEMENTED")
        # This system MF can be parameterized by the same type of lognormal form as
        # the single MF (eq. [17]), with the same normalization at 1 Msun, with the
        # coefficients (Chabrier 2003)
        return 0.86 * np.exp(-1*(np.log10(m)-np.log10(0.22))**2/(2*0.57**2))
        # This analytic form for the disk MF for single objects below 1 Msun, within these uncertainties, is given by the following lognormal form (Chabrier 2003):
        return 0.158 * np.exp(-1*(np.log10(m)-np.log10(0.08))**2/(2*0.69**2))

    def schechter(self,m,A=1,beta=2,m0=100, integral=False,**kwargs):
        """
        A Schechter function with arbitrary defaults
        (integral may not be correct - exponent hasn't been dealt with at all)

        $$ A m^{-\\beta} e^{-m/m_0} $$
        
        Parameters
        ----------
            m : np.ndarray
                List of masses for which to compute the Schechter function
            A : float
                Arbitrary amplitude of the Schechter function
            beta : float
                Power law exponent
            m0 : float
                Characteristic mass (mass at which exponential decay takes over)
        
        Returns
        -------
            p(m) - the (unnormalized) probability of an object of a given mass
            as a function of that object's mass
            (though you could interpret mass as anything, it's just a number)

        """
        if integral: beta -= 1
        return A*m**-beta * np.exp(-m/m0)

    def modified_schechter(self,m,m1=0.5, **kwargs):
        """
        A Schechter function with a low-level exponential cutoff
        "
        Parameters
        ----------
            m : np.ndarray
                List of masses for which to compute the Schechter function
            m1 : float
                Characteristic minimum mass (exponential decay below this mass)
            ** See schecter for other parameters ** 

        Returns
        -------
            p(m) - the (unnormalized) probability of an object of a given mass
            as a function of that object's mass
            (though you could interpret mass as anything, it's just a number)
        """
        return self.schechter(m, **kwargs) * np.exp(-m1/m)

    def inverse_imf(self,p,nbins=1000,**kwargs):
        """
        Inverse mass function
        """
     
        masses = np.logspace(np.log10(self.mmin),np.log10(self.mmax),nbins)
        if self.alpha is not None:
            mf = self.imf(masses,alpha=self.alpha,integral=True,**kwargs)
        else:
            mf = self.imf(masses,integral=True,**kwargs)
        mfcum = mf.cumsum()
        mfcum /= mfcum.max() # normalize to sum (cdf)

        return np.interp(p,mfcum,masses)

    def make_cluster_masses(self,num_stars,nbins=1000,**kwargs):
        """
        Sample from an IMF to make a cluster.  Returns the masses of all stars in the cluster
        
        kwargs are passed to `inverse_imf`
        """

        masses = self.inverse_imf(np.random.RandomState(seed=self.seed).random_sample(num_stars), nbins=nbins, **kwargs)

        return masses
    
    def get_star_info(self):
        stmt = "SELECT m_ini,Te,log_g,johnson_i_abs FROM star_info WHERE star_info.age = %f AND star_info.Z = %f" % (self.age,self.metallicity)
        sc = sqlite_util.SqliteConnection(self.dbname)
        arr = sc.execute(stmt,asarray=True)
        return arr['m_ini'],arr['te'],arr['log_g'],arr['johnson_i_abs']
    
    def make_cluster_rates(self,masses,instrument,filter,bandpass=None,refs=None):
        try:
            coords = np.load(os.path.join(self.gridpath, 'input.npy'), allow_pickle=True)
        except UnicodeError:
            coords = np.load(os.path.join(self.gridpath, 'input.npy'), allow_pickle=True, encoding='bytes')
        m, t, g, i = self.get_star_info()
        temps = np.interp(masses,m,t)
        gravs = np.interp(masses,m,g)
        mags = np.interp(masses,m,i)
        metals = np.full_like(mags, self.metallicity)
        if os.path.exists(os.path.join(self.gridpath, 'result_{}_{}.npy'.format(instrument.lower(), filter.lower()))):
            values = np.load(os.path.join(self.gridpath, 'result_{}_{}.npy'.format(instrument.lower(), filter.lower())), allow_pickle=True)
            interpolation_function = RegularGridInterpolator(tuple([x for x in coords]), values)
            try:
                countrates = interpolation_function(np.array((metals, gravs, temps, mags)).T)
            except ValueError as v:
                self.log('error', 'Exception caught when interpolating: {}'.format(v))
                min_mag = coords[-1][0]
                max_mag = coords[-1][-1]
                interpolation_function = RegularGridInterpolator(tuple([x for x in coords]), values, bounds_error=False, fill_value=0.)
                mags_min = np.full_like(mags, min_mag)
                mags_max = np.full_like(mags, max_mag)
                countrates = interpolation_function(np.array((metals, gravs, temps, mags)).T)
                countrates_min = interpolation_function(np.array((metals, gravs, temps, mags_min)).T)
                countrates_min = countrates_min * np.power(10, -(mags-min_mag)/2.512)
                countrates_max = interpolation_function(np.array((metals, gravs, temps, mags_max)).T)
                countrates_max = countrates_max * np.power(10, -(mags-max_mag)/2.512)
                countrates[np.where(mags < mags_min)] = countrates_min[np.where(mags < mags_min)]
                countrates[np.where(mags > mags_max)] = countrates_max[np.where(mags > mags_max)]
        else:
            self.log('warning', 'Could not find result file "result_{}_{}.npy" from {}'.format(instrument.lower(), filter.lower(), self.gridpath))
#             raise FileNotFoundError('Could not find result file "result_{}_{}.npy" from {}'.format(instrument.lower(), filter.lower(), self.gridpath))
            import pysynphot as ps
            countrates = np.array(())
            ps.setref(**refs)
            johnson_i = ps.ObsBandpass('johnson,i')
            for te, log_g, z, j_i in zip(temps, gravs, metals, mags):
                spectrum = ps.Icat('phoenix', te, z, log_g)
                spectrum = spectrum.renorm(j_i, 'vegamag', johnson_i)
                obs = ps.Observation(spectrum, bandpass, binset=spectrum.wave)
                countrates = np.append(countrates, obs.countrate())
                self.log('info', 'Manually created star {} of {}'.format(len(countrates), len(temps)))
        return countrates

    def make_cluster_mags(self,masses):
        db = sqlite3.connect(self.dbname)
        sc = sqlite_util.SqliteConnection(self.dbname)
        stmt = "SELECT m_ini,johnson_i_abs,Te,log_g FROM star_info WHERE star_info.age = %f AND star_info.Z = %f" % (self.age,self.metallicity)
        arr = sc.execute(stmt, asarray=True)
        star_masses = arr['m_ini']
        star_i = arr['johnson_i_abs']
        star_t = arr['te']
        star_g = arr['log_g']
        countrates = np.interp(masses, star_masses, star_i)
        temps = np.interp(masses, star_masses, star_t)
        gravs = np.interp(masses, star_masses, star_g)
        return (countrates, temps, gravs)

    def make_cluster(self,num_stars,nbins=1000,**kwargs):
        """
        Sample from an IMF to make a cluster. Returns the masses and count rates of all stars 
        in the cluster
        
        kwargs are passed to `inverse_imf`
        """

        masses = self.make_cluster_masses(num_stars,nbins=nbins,**kwargs)
        mags, temps, gravs = self.make_cluster_mags(masses)

        return masses, mags, temps, gravs
    
    def log(self, priority, message):
        if self.logger is not None:
            getattr(self.logger, priority)(message)
        else:
            print("{}: {}".format(priority, message))

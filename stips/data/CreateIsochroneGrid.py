#!/usr/bin/env python

import glob
import logging
import math
import mechanize
import numpy as np
import os
import re
import sys
import sqlite3
import warnings
from argparse import ArgumentParser

def getZeta(Z,base_zeta = 0.019):
    """
    Returns the fractional metal content of a star based on [M/H]. Truncates the range of possible
    values based on the limits of CMD 2.5.
    
    Inputs:
    -------
    Z: floating point.
    
        The [M/H] value to convert
    
    base_zeta: floating point, optional, default = 0.019
    
        The base of the [M/H] scale. Default is 0.019 for a solar base
    
    Output:
    -------
    
    zeta: floating point
    
        Metal content of the star in question
    """
    zeta_limits = [0.0001,0.06]
    zeta = 10**Z * base_zeta
    if zeta < zeta_limits[0]:
        zeta = zeta_limits[0]
    if zeta > zeta_limits[1]:
        zeta = zeta_limits[1]
    return zeta

def getZ(zeta,base_zeta = 0.019):
    """
    Returns [M/H] based on the fractional metal content of a star.
    
    Inputs:
    -------    
    zeta: floating point
    
        Metal content of the star.

    base_zeta: floating point, optional, default = 0.019
    
        The base of the [M/H] scale. Default is 0.019 for a solar base
    
    Output:
    -------
    
    Z: floating point.
    
        The [M/H] value.
    """
    Z = math.log10(zeta / base_zeta)
    return Z

def divideInterval(divisor,reverse=False):
    """
    Takes an interval expressed as a string in the form "low,high,interval" and
        then performs black magic to turn it into a linear array.
        
        In particular, depending on the value of interval:
        
            interval is an integer X:
            interval is 'n' followed by an integer X:
            
                Use numpy.linspace(low,high,interval) to produce X intervals,
                    including both low and high, and evenly divided between them.
            
            interval is 'i' followed by an integer X:
            
                Use numpy.linspace(low,high,((high-low)/interval)+1) to produce
                    intervals *of* X from low to high, and including both. Uses
                    'linspace' rather than 'arange' because 'linspace' does not
                    produce weird effects at start and end points, and will
                    actually make zero 0. rather than 1.58395e-15 or some such.
            
            interval is 'd':
            
                Use numpy.linspace() called multiple times to produce decades.
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
            
                Use numpy.linspace() called multiple times to produce half-decades.
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

def getStars(age, zeta):
    """
    Returns a list of stars (for use by e.g. phoenix grid) with necessary parameters to determine
    their characteristics and magnitude.
    
    Inputs:
    -------
    
    age: float
    
        The age of the stellar population, in years
    
    zeta: float
    
        The metal content of the stars, where the Sun has 0.019 as its metal content.
    
    Output:
    -------
    
    star_list: list of lists
    
        The output stars. The list has the following columns:
            Z: [M/H]
            Age (Gyr)
            M_ini (M_\odot)
            M_act (M_\odot)
            Te
            log(g)
            int_IMF
            Johnson,I
    """
    result_str = re.compile(r"The results are available at <a href=(.*?)>output\d*\.dat</a>")
    request = mechanize.Request("http://stev.oapd.inaf.it/cgi-bin/cmd_2.5")
    response = mechanize.urlopen(request)
    forms = mechanize.ParseResponse(response,backwards_compat=False)
    response.close()
    form = forms[0]
    #The reasoning here is that I can *get* Johnson filters in JWST pysynphot, but can't figure
    #out how to do the Spitzer equivalents.
#    form["photsys_file"] = ["tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat"]
    form["photsys_file"] = ["tab_mag_odfnew/tab_mag_ubvrijhk.dat"]
    #0 = single isochrone, single metallicity.
    #1 = sequence of isochrones, different ages, constant metallicity
    #2 = sequence of isochrones, constant age, different metallicities
    form["isoc_val"] = ["0"]
    #Age for single-single
    form["isoc_age"] = '%g' % (age)
    form["isoc_zeta"] = '%g' % (zeta)
    request2 = form.click()
    response2 = mechanize.urlopen(request2)
    response_value = response2.read()
    response_url = response2.geturl()
    match = result_str.search(response_value)
    star_list = []
    if match is not None:
        output_url = match.group(1)
        response_result = mechanize.urlopen(mechanize.urljoin(response_url,output_url))
        output_lines = response_result.read().split("\n")
        output_lines = output_lines[13:]
        for line in output_lines:
            if line != "":
                #Z, log(age/yr), M_ini, M_act, logL/Lo, logTe, logG, mbol, U, B, V, R, I, J, H, K, int_IMF, stage
                items = line.split()
                star = [None]
                star.append(getZ(float(items[0])))
                star.append(10**float(items[1]))
                star.append(float(items[2]))
                star.append(float(items[3]))
                star.append(10**float(items[5]))
                star.append(float(items[6]))
                star.append(float(items[6]))
                star.append(float(items[16]))
                star.append(float(items[12]))
                star_list.append(star)
    return star_list

def makeDB(dbname):
    """
    Makes a database to hold information about the various stars and their JWST count rates
    
    Inputs:
    -------
    
    dbname: string
    
        Name of the database to create
    """
    if os.path.exists(dbname): #Remove database if it already exists
        os.remove(dbname)
    db = sqlite3.connect(dbname)
    c = db.cursor()
    c.execute("""SELECT name FROM sqlite_master WHERE name = 'star_info'""")
    results = c.fetchall()
    if len(results) == 0:
        c.execute("""CREATE TABLE star_info
                        (id integer primary key,
                         Z real,
                         age real,
                         m_ini real,
                         m_act real,
                         Te real,
                         log_g real,
                         log_g_orig real,
                         int_imf real,
                         johnson_i_abs real)""")
    db.commit()
    db.close()

def makeStars(dbname,Zs,ages,null,stdout,verbose=False,logfile="",override=False,erase_existing=False,reverse=False):
    """
    Stores stellar count rate in all JWST filters in a database.
    
    Inputs:
    -------
    
    dbname: string
    
        Path to the stellar database
    
    Zs: string
    
        The divideInterval() string of metallicities for the grid.
    
    ages: string
    
        The divideInterval() string of ages for the grid.
    
    verbose: boolean, optional, default=False
    
        Provide additional logging information
    
    logfile: string, optional, default=""
    
        File to store the logging information in (in addition to printing it)
    
    override: boolean, optional, default=False
    
        If true, will store entries in database even if they exist already.
    
    erase_existing: boolean, optional, default=False
    
        If true, will delete existing duplicate entries before storing the new version. If
        override = False, delete_existing will be ignored.
        
    reverse: boolean, optional, default=False
    
    	If true, stellar grid will be run in reverse order
    """
    if not os.path.exists(dbname):
        makeDB(dbname)
    db = sqlite3.connect(dbname)
    if verbose:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(message)s')
        stdout_handler = logging.StreamHandler()
        logger.addHandler(stdout_handler)
        stdout_handler.setFormatter(formatter)
        if logfile != "":
            file_handler = logging.FileHandler(logfile)
            logger.addHandler(file_handler)
            file_handler.setFormatter(formatter)
    for Z in divideInterval(Zs,reverse):
        if verbose:
            logger.info("Handling Z = %g",Z)
        for age in divideInterval(ages,reverse):
            if verbose:
                logger.info("\tHandling age = %g",age)
            star_list = getStars(age,getZeta(Z))
            if verbose:
                logger.info("\tGot list of stars for %g %g",Z,age)
#                          0     1    2      3      4   5       6       7        8          9
            # star_list = [None, Z, Age, M_ini, M_act, Te, log(g), log(g), int_IMF, Johnson,I]
            for star in star_list:
                star[1] = float("{:.2e}".format(star[1]))
                star[2] = float("{:.2e}".format(star[2]))
                c = db.cursor()
                
                c.execute("""SELECT * FROM star_info WHERE Z = ? AND age = ? AND m_ini = ? AND
                                m_act = ? AND Te = ? AND log_g = ? AND log_g_orig = ? AND
                                int_imf = ? AND johnson_i_abs = ?""",tuple(star[1:]))
                results = c.fetchall()
                if len(results) == 0 or override: #Add to database if not there already
                
                    if len(results) != 0 and erase_existing: #Delete existing entries to avoid confusion.
                        for result in results:
                            idnum = result[0]
                            c.execute("""DELETE FROM star_info WHERE id = ?""",(idnum,))
                            logger.info("\t\tREMOVING FROM DATABSE INDEX %d: star %g %g %g %g %g %g %g",idnum,*tuple(star[3:]))
                        db.commit()

                    # ***** THIS IS A HACK *****
                    # The Phoenix grid can't handle log(g) < 0
                    # As such, I'm altering negative values to be equal to zero.
                    # I'm also storing them in the database that way, so that taking something out of
                    #   the database and running it through the phoenix grid will get exactly the same
                    #   results that were stored.
                    if star[6] < 0:
                        star[6] = 0.

                    # 0                       1  2    3      4      5   6      7           8        9
                    # id integer primary key, Z, age, m_ini, m_act, Te, log_g, log_g_orig, int_imf, johnson_i_abs
                    c.execute("""INSERT INTO star_info VALUES (?,?,?,?,?,?,?,?,?,?)""", tuple(star))
                    db.commit()
                else:
                    logger.info("\t\tALREADY IN DATABASE: Star %g %g %g %g %g %g %g",*tuple(star[3:]))
    db.close()
    if verbose:
        logf.close()

if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with open(os.devnull, 'w') as f:
            stdout = sys.stdout
        
            parser = ArgumentParser()
            parser.add_argument("-d","--database",action="store",dest="db",default="IsochroneGrid.db",
                                metavar="DATABASE",help="Database Name")
            parser.add_argument("-m","--metallicity",action="store",dest="metallicity",default="-2.5,0.5,i0.1",
                                metavar="MIN,MAX,NSTEPS",help="Metallicity axis")
            parser.add_argument("-a","--age",action="store",dest="age",default="1.0e6,1.35e10,d5",
                                metavar="MIN,MAX,NSTEPS",help="Age axis")
            parser.add_argument("-l","--logfile",action="store",dest="logfile",default="",
                                metavar="FILE",help="Output logging file")
            parser.add_argument("-o","--override",action="store_true",dest="override",default=False,
                                help="Store values in database even if they already exist")
            parser.add_argument("-e","--erase_existing",action="store_true",dest="erase_existing",default=False,
                                help="If override is active, delete duplicate entries from database before storing new entries")
            parser.add_argument("-r","--reverse_order",action="store_true",dest="reverse",default=False,
    		    		    	help="Reverse the order of traversing the stellar grid")
            parser.add_argument("-v","--verbosity",action="store_true",dest="verbosity",default=False,
                                help="Verbose output")
            results = parser.parse_args()
    
            makeStars(results.db,results.metallicity,results.age,f,stdout,verbose=results.verbosity,
                      logfile=results.logfile,override=results.override, erase_existing=results.erase_existing,
                      reverse=results.reverse)

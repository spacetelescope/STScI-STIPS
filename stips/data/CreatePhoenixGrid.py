from __future__ import print_function

import glob
import os
import sys
import time

# Set plot backend. Must happen before importing pylab.
import matplotlib
matplotlib.use('Agg')

pandeia_version_file = os.path.join(os.environ["pandeia_refdata"], "VERSION_PSF")
with open(pandeia_version_file, 'r') as inf:
    pandeia_version_info = inf.readline().strip()
print("Pandeia Version: {}".format(pandeia_version_info))

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), "..", "..")))

import numpy as np
import pysynphot as ps
from stips.utilities import InstrumentList
from .. import __version__ as stips_version_info
from astropy.io import fits as pyfits

print("STIPS Version: {}".format(stips_version_info))

COMPFILES =  sorted(glob.glob(os.path.join(os.environ["PYSYN_CDBS"],"mtab","*tmc.fits")))
GRAPHFILES = sorted(glob.glob(os.path.join(os.environ["PYSYN_CDBS"],"mtab","*tmg.fits")))
THERMFILES = sorted(glob.glob(os.path.join(os.environ["PYSYN_CDBS"],"mtab","*tmt.fits")))
wfc3_refs = {
                'comptable': COMPFILES[-1],
                'graphtable': GRAPHFILES[-1],
                'thermtable': THERMFILES[-1],
                'area': 45238.93416,
                'waveset': (500,26000,10000.,'log')
            }

area =  {
            'nircamlong':   253260.0,
            'nircamshort':  253260.0,
            'miri':         253260.0,
            'wfc3ir':       45238.93416,
            'wfi':          45238.93416
        }
modes = {
            'nircamshort':  ['sw_imaging'],
            'nircamlong':   ['lw_imaging'],
            'miri':         ['imaging'],
            'wfc3ir':       ['imaging'],
            'wfi':          ['imaging']
        }
filters = {
            'nircamshort':  {'sw_imaging': ["f070w", "f090w", "f115w", "f140m", "f150w", "f162m", "f164n", "f182m", "f187n", "f200w", "f210m", "f212n"]},
            'nircamlong':   {'lw_imaging': ["f250m", "f277w", "f300m", "f323n", "f335m", "f356w", "f360m", "f405n", "f410m", "f430m", "f444w", "f460m", "f466n", "f470n", 
                                            "f480m"
                                            ]
                            },
            'miri':         {'imaging': ["f560w", "f770w", "f1000w", "f1130w", "f1280w", "f1500w", "f1800w", "f2100w", "f2550w"]},
            'wfc3ir':       {'imaging': ["f110w", "f160w"]},
            'wfi':          {'imaging': ['f062', 'z087', 'y106', 'w149', 'j129', 'h158', 'f184']}
          }
apertures = {
                'nircamshort':  {'sw_imaging':  "sw"},
                'nircamlong':   {'lw_imaging':  "lw"},
                'miri':         {'imaging':     "imager"},
                'wfc3ir':       {'imaging':     "default"},
                'wfi':          {'imaging':     "any"}
            }

def get_grid_points():
    grid_file = os.path.abspath(os.path.join(os.environ['PYSYN_CDBS'], "grid", "phoenix", "catalog.fits"))
    teff, Z, logg = np.array(()), np.array(()), np.array(())
    with pyfits.open(grid_file) as inf:
        indices = inf[1].data.field('INDEX')
        for row in indices:
            items = row.split(",")
            teff = np.append(teff, float(items[0]))
            Z = np.append(Z, float(items[1]))
            logg = np.append(logg, float(items[2]))
    return np.array((np.unique(Z), np.unique(logg), np.unique(teff), np.arange(-5.5, 16.0)))


if __name__ == '__main__':
    norm_bandpass = ps.ObsBandpass('johnson,i')
    coords = get_grid_points()
    print(coords)
    bandpasses = {}
    result_arrays = {}

    instruments = InstrumentList()
    print("{}: Making Bandpasses...".format(time.ctime()))
    for instrument in instruments:
        my_instrument = instruments[instrument](log_level="WARNING")
        bandpasses[instrument.lower()] = {}
        result_arrays[instrument.lower()] = {}
        for mode in modes[instrument.lower()]:
            for filter in filters[instrument.lower()][mode]:
                print("\t{}: {},{},{},{}".format(time.ctime(), instrument, mode, filter, apertures[instrument.lower()][mode]))
                my_instrument.reset(0., 0., 0., filter.upper(), 0., psf=False, detectors=False)
                bandpasses[instrument.lower()][filter] = my_instrument.bandpass
                result_arrays[instrument.lower()][filter] = np.empty((len(coords[0]), len(coords[1]), len(coords[2]), len(coords[3])))
    print("Done\n")

    total = len(coords[0]) * len(coords[1]) * len(coords[2]) * len(coords[3])
    n = 0
    for i, Z in enumerate(coords[0]):
        print("{}: Starting Z = {}".format(time.ctime(), Z))
        for j, logg in enumerate(coords[1]):
            print("\t{}: Starting log(g) = {}".format(time.ctime(), logg))
            for k, teff in enumerate(coords[2]):
                print("\t\t{}: Starting Teff = {}".format(time.ctime(), teff))
                try:
                    spec = ps.Icat('phoenix', teff, Z, logg)
                    counts = False
                    if sum(spec.flux) > 0:
                        counts = True
                except ps.exceptions.ParameterOutOfBounds:
                    counts = False
                for l, mag in enumerate(coords[3]):
                    print("\t\t\t{} of {}: {}: Starting Z = {}, log(g) = {}, Teff = {}, Mabs = {:>4}".format(n+1, total, time.ctime(), Z, logg, teff, mag), end='')
                    if counts:
                        spec_norm = spec.renorm(mag, 'vegamag', norm_bandpass)
                    for instrument in instruments:
                        ps.setref(area=area[instrument.lower()], waveset=(500, 260000, 10000., 'log'))
                        for mode in modes[instrument.lower()]:
                            for filter in filters[instrument.lower()][mode]:
                                if counts:
                                    obs = ps.Observation(spec_norm, bandpasses[instrument.lower()][filter], binset=spec_norm.wave)
                                    result_arrays[instrument.lower()][filter][i,j,k,l] = obs.countrate()
                                    print(".", end='')
                                else:
                                    result_arrays[instrument.lower()][filter][i,j,k,l] = 0.
                                    print("x", end='')
                    print("")
                    n += 1


    print("{}: Saving files...".format(time.ctime()), end='')
    with open(os.path.join(os.getcwd(), "grid", "VERSION.txt"), "wt") as outf:
        outf.write("Pandeia: {}\n".format(pandeia_version_info))
        outf.write("STIPS: {}\n".format(stips_version_info))
    np.save(os.path.join(os.getcwd(), 'grid', 'input.npy'), coords)
    for instrument in instruments:
        for mode in modes[instrument.lower()]:
            for filter in filters[instrument.lower()][mode]:
                np.save(os.path.join(os.getcwd(), 'grid', 'result_{}_{}.npy'.format(instrument.lower(), filter)),
                        result_arrays[instrument.lower()][filter])
    print("done")

from __future__ import division

import inspect,numpy,os,pytest,sys
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
topdir = os.path.dirname(parentdir)
sys.path.insert(0,parentdir)

import astropy.io.fits as pyfits
from astropy.table import Table,Column

from SceneModule import SceneModule

def verifyParameters(scene,results):
    assert scene.in_path == results['in_path']
    assert scene.out_path == results['out_path']
    assert scene.prefix == results['prefix']
    assert scene.seed == results['seed']
    assert scene.ra == results['ra']
    assert scene.dec == results['dec']
    assert scene.params == results['params']
    assert scene.catalogues == results['catalogues']

creation_data =         [
                            (
                                {},
                                {'in_path':os.getcwd(),'out_path':os.getcwd(),'prefix':'sim',
                                 'seed':0, 'ra': 0., 'dec': 0., 'catalogues':{},
                                 'params':['Random seed: 0','Centre (RA,DEC) = (0.000000,0.000000)']}
                            ),
                            (
                                {'in_path':'/path/to/somewhere','out_path':'/path/to/somewhere/else',
                                 'out_prefix':'sim123','general':{'seed':4321,'ra':11.5,'dec':-8.7}},
                                {'in_path':'/path/to/somewhere','out_path':'/path/to/somewhere/else',
                                 'prefix':'sim123','seed':4321, 'ra':11.5, 'dec':-8.7, 'catalogues':{},
                                 'params':['Random seed: 4321','Centre (RA,DEC) = (11.500000,-8.700000)']}
                            )
                        ]

@pytest.mark.parametrize(("inputs","results"), creation_data)
def test_creation(inputs,results):
    scene = SceneModule(**inputs)
    verifyParameters(scene,results)

createPopulation_data = [
                            (
                                {'n_stars':500,'age_low':1.e9,'age_high':1.e9,'z_low':0.,'z_high':0.,
                                 'imf':'kroupa','alpha':-2.5,'distribution':'uniform',
                                 'clustered':False,'radius':88.5,'radius_units':'pc',
                                 'distance_low':5.,'distance_high':15.,'binary_fraction':0.2,
                                 'offset_ra':0.,'offset_dec':0.},
                                3333,
                                os.path.join(currentdir,"SceneModule_createPopulation_01")
                            ),
                            (
                                {'n_stars':50000,'age_low':1.e11,'age_high':1.e12,'z_low':-2.5,
                                 'z_high':-2.0,'imf':'salpeter','alpha':-2.5,
                                 'distribution':'invpow','clustered':True,'radius':34.5,
                                 'radius_units':'arcsec','distance_low':18.,'distance_high':18.,
                                 'binary_fraction':0.05,'offset_ra':0.15,'offset_dec':-0.17},
                                7876,
                                os.path.join(currentdir,"SceneModule_createPopulation_02")
                            )
                        ]

@pytest.mark.parametrize(("pop","seed","fname"), createPopulation_data)
def test_createPopulation(pop,seed,fname):
    scene = SceneModule(seed=seed,in_path=os.path.join(topdir,"sim_input"))
    pop_file = scene.CreatePopulation(pop)
    created_pop = open(pop_file,mode='r')
    created_lines = created_pop.readlines()
    created_pop.close()
    reference_pop = open(fname,mode='r')
    reference_lines = reference_pop.readlines()
    reference_pop.close()
    for i in range(len(created_lines)):
        assert created_lines[i] == reference_lines[i]
    os.remove(pop_file)

createGalaxies_data =         [
                            (
                                {'n_gals':500,'z_low':0.,'z_high':1.,'rad_low':1.,'rad_high':5.,
                                 'vmag_low':23.,'vmag_high':28.,'distribution':'uniform',
                                 'clustered':False,'radius':88.5,'radius_units':'pc','offset_ra':0.,
                                 'offset_dec':0.},
                                3333,
                                os.path.join(currentdir,"SceneModule_createGalaxies_01")
                            ),
                            (
                                {'n_gals':500,'z_low':2.,'z_high':5.,'rad_low':0.5,'rad_high':3.5,
                                 'vmag_low':21.,'vmag_high':29.,'distribution':'invpow',
                                 'clustered':True,'radius':32.6,'radius_units':'arcsec',
                                 'offset_ra':2.5,'offset_dec':1.3},
                                7876,
                                os.path.join(currentdir,"SceneModule_createGalaxies_02")
                            )
                        ]

@pytest.mark.parametrize(("gal","seed","fname"), createGalaxies_data)
def test_createGalaxies(gal,seed,fname):
    scene = SceneModule(seed=seed,in_path=os.path.join(topdir,"sim_input"))
    gal_file = scene.CreateGalaxies(gal)
    created_gal = open(gal_file,mode='r')
    created_lines = created_gal.readlines()
    created_gal.close()
    reference_gal = open(fname,mode='r')
    reference_lines = reference_gal.readlines()
    reference_gal.close()
    for i in range(len(created_lines)):
        assert created_lines[i] == reference_lines[i]
    os.remove(gal_file)

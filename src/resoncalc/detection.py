# -*- coding: utf-8 -*
"""Detection

.. module:: detection
   :synopsis: Module for detection of eigenstates

"""

# imports
import numpy as np
import resoncalc.femdvr as femdvr
import resoncalc.potential as potential
import resoncalc.output as output
from time import time
from os import path
from scipy.signal import find_peaks

# globals
nquad = 15            # order of quadrature
emax = 1.0            # maximum energy in atomic units
mu = 1.0              # reduced mass in atomic units
l = 1                 # secondar quantum number
x0 = 0.0              # ecs point x0
phases = [40.0, 30.0] # ecs phases
prec = 1e-8           # ecs precision

def perform_detection_loop(cfg):
    """Perform eigenstate detection in loop using potential parameters variation 

    Args:
        cfg (dict): configuration

    sample
    {
     'title': 'test',
     'potential': 'gaussian',
     'nquad' : 15,
     'x0' : 0.0,
     'phases' : [40.0, 30.0],
     'prec' : 1e-8,
     'emax' : 1.0,
     'mu' : 1.0,
     'l' : 1,
     'intervals' : [
                    {'a': 0.0,   'b': 9.0,     'elems': 15, 'type': 'equidistant'},
                    {'a': 9.0,   'b': 100.0,   'elems': 15, 'type': 'progressive', 'len': 0.6},
                    {'a': 100.0, 'b': 10000.0, 'elems': 15, 'type': 'progressive', 'len': 6.0}
                   ],
     'params' : [
                 {'a': -0.620, 'b': -0.560, 'cnt': 13},
                 {'a': 0.1,    'b': 0.1,    'cnt': 1}
                ],
     'params2' : [0.7416],
     'outstates' : ['bound', 'resonance'],
     'outfiles' : ['eigenstates', 'eigenvalues', 'potential_grid', 'spectrum', 'log', 'test']
    }

    """     

    # ranges
    rng1 = np.linspace(cfg['params'][0]['a'], cfg['params'][0]['b'], cfg['params'][0]['cnt'])
    rng2 = np.linspace(cfg['params'][1]['a'], cfg['params'][1]['b'], cfg['params'][1]['cnt'])
    total = len(rng1) * len(rng2)

    # optional parameters
    title = cfg['title'] if ('title' in cfg) else cfg['potential']
    globals()['nquad'] = cfg['nquad'] if ('nquad' in cfg) else globals()['nquad']
    globals()['emax'] = cfg['emax'] if ('emax' in cfg) else globals()['emax']
    globals()['mu'] = cfg['mu'] if ('mu' in cfg) else globals()['mu']
    globals()['l'] = cfg['l'] if ('l' in cfg) else globals()['l']
    globals()['x0'] = cfg['x0'] if ('x0' in cfg) else globals()['x0']
    globals()['phases'] = cfg['phases'] if ('phases' in cfg) else globals()['phases'] 
    globals()['prec'] = cfg['prec'] if ('prec' in cfg) else globals()['prec']
    params2 = cfg['params2'] if ('params2' in cfg) else []    
    
    if ('outstates' in cfg):
        output.outstates = cfg['outstates']
    if ('outfiles' in cfg):
        output.outfiles = cfg['outfiles']

    # use all combinations of potential parameters
    start = time()
    output.set_outdir(title)
    output.export_test(cfg)
    output.info('Running test: {0} for potential: {1}, total tests: {2}'.format(title, cfg['potential'], total))
    test = 1
    eigenstates = []
    params = []

    try:        

        # create common grid
        grid = create_grid(nquad, cfg['intervals'])
        
        # param1 loop
        for param1 in rng1:
            # param2 loop
            for param2 in rng2:

                output.info('Test: {0}/{1}, parameters: {2}, {3}'.format(test, total, round(param1, 5), round(param2, 5)))

                # detection
                try:
                    states = perform_detection(potential.mapping[cfg['potential']], nquad, emax, mu, cfg['intervals'], grid,
                                               param1, param2, l, mu, *params2)

                # exceptions handled outside
                except KeyboardInterrupt as ex:
                    raise KeyboardInterrupt
                except FileNotFoundError as ex:
                    raise FileNotFoundError
                
                # logged exception
                except Exception as ex:
                    states = None
                    output.error('Exception {0}'.format(ex))

                # states detected
                if (states is not None):
                    for state in states:
                        eigenstates.append(state)
                        params.append((param1, param2, l))
                
                test += 1

        # export eigenstates
        output.export_eigenstates(params, eigenstates)
        output.info('Tests finished after {0}s'.format(round(time() - start, 1)))

    # export already detected eigenstates    
    except KeyboardInterrupt as ex:
        output.export_eigenstates(params, eigenstates)
        raise KeyboardInterrupt

    # output directory deleted
    except FileNotFoundError as ex:
        print('Output directory {0} not found'.format(output.outdir))

def perform_detection(potential, n, emax, mu, intervals, grid, *params):
    """Perform eigenstate detection for one potential configuration

    Args:
        potential (func): potential function
        n (int): order of quadrature
        emax (float): maximum energy in atomic units
        mu (float): reduced mass in atomic units
        intervals (dict): interval configuration             
        params (args): potential specific parameters

    Returns:
        list: eigenstates

    """    

    # spectrum for first ECS phase
    phase = phases[0]
    output.debug('Using ECS x0: {0}, phase: {1}'.format(x0, phase))
    output.debug('Creating Hamiltonian matrix')
    matrix, potential_grid = femdvr.hamiltonian(n, grid['endpoints1'], potential, grid['points1'], grid['stiffness1'], mu, *params)
    output.export_potential_grid(grid['points1'], potential_grid, phase, *params)
    output.debug('Calculating eigenvalues')
    eigenvalues1 = femdvr.spectrum(matrix)
    output.export_eigenvalues(eigenvalues1, phase, *params)

    # spectrum for second ECS phase
    phase = phases[1]
    output.debug('Using ECS x0: {0}, phase: {1}'.format(x0, phase))
    output.debug('Creating Hamiltonian matrix')
    matrix, potential_grid = femdvr.hamiltonian(n, grid['endpoints2'], potential, grid['points2'], grid['stiffness2'], mu, *params)
    output.export_potential_grid(grid['points2'], potential_grid, phase, *params)
    output.debug('Calculating eigenvalues')
    eigenvalues2 = femdvr.spectrum(matrix)
    output.export_eigenvalues(eigenvalues2, phase, *params)

    # detect eigenstate
    output.debug('Detecting eigenstate')
    init_point = grid['endpoints1'][1]
    potmax = potential_maximum(potential, np.real(init_point), np.real(init_point)+20.0, *params)
    states = detect_states(emax, potmax, eigenvalues1, eigenvalues2)
    
    if (states is not None):
        for state in states:
            output.info('Detected eigenstate real: {0:e}, imag: {1:e}, type: {2}'.format(np.real(state), np.imag(state), get_state_type(state)))

        # plot spectrum
        output.export_complex_spectrum_fig(eigenvalues1, eigenvalues2, states, *params)

        # plot eigenstates 
        output.export_eigenstates_fig(potential, states, np.real(init_point), np.real(init_point)+20.0, emax, *params)
    else:
        output.info('Eigenstate not detected')

    return states

def create_grid(n, intervals):
    """Create common grid and stiffness matrix for both ECS phases

    Args:
        n (int): order of quadrature
        intervals (dict): interval configuration

    Returns:
        dict: endpoints1, points1, stiffness1
              endpoints2, points2, stiffness2

    """     

    output.debug('Creating grid')
    grid = {}

    # first ECS phase
    endpoints = create_endpoints(intervals)
    endpoints = femdvr.exterior_complex_scaling(endpoints, x0, phases[0])
    points, weights = femdvr.fem_points_weights(n, endpoints)
    stiffness = femdvr.stiffness_matrix(n, endpoints)
    grid['endpoints1'] = endpoints
    grid['points1'] = points
    grid['stiffness1'] = stiffness

    # second ECS phase
    endpoints = create_endpoints(intervals)
    endpoints = femdvr.exterior_complex_scaling(endpoints, x0, phases[1])
    points, weights = femdvr.fem_points_weights(n, endpoints)
    stiffness = femdvr.stiffness_matrix(n, endpoints)
    grid['endpoints2'] = endpoints
    grid['points2'] = points
    grid['stiffness2'] = stiffness

    return grid
                
def create_endpoints(intervals):
    """Create endpoints from given intervals

    Args:
        intervals (dict): interval configuration

    Returns:
        list: endpoints

    """      

    endpoints = []
    for interval in intervals:

        # equistant endpoints
        if (interval['type'] == 'equidistant'):
            endpoints.append(femdvr.equidistant_endpoints(interval['a'], interval['b'], interval['elems']))
        # progressive endpoints
        elif (interval['type'] == 'progressive'):
            endpoints.append(femdvr.progressive_endpoints(interval['a'], interval['b'], interval['elems'], interval['len']))

    # merge intervals
    endpoints = femdvr.merge_endpoints(endpoints)
    
    return endpoints

def detect_states(emax, potmax, eigenvalues1, eigenvalues2):
    """Detect eigenstates using ECS method executed twice

    Args:
        emax (float): maximum energy in atomic units
        potmax (float): potential maximum
        eigenvalues1 (list): eigenvalues for first phase
        eigenvalues2 (list): eigenvalues for second phase

    Returns:
        list: eigenstates

    """    

    # select eigenvalues up to maximum energy
    vals1, vals2 = [], []
    maximum = min(emax, potmax)
    for i in range(len(eigenvalues1)):
        if (np.real(eigenvalues1[i]) <= maximum):
            vals1.append(eigenvalues1[i])
        if (np.real(eigenvalues2[i]) <= maximum):
            vals2.append(eigenvalues2[i])

    # find eigenstate independent on phase    
    states = []
    for i in range(len(vals1)):
        for j in range(len(vals2)):
            if (np.abs(vals1[i] - vals2[j]) < prec):
                states.append((vals1[i] + vals2[j])/2.0)

    # not detected
    if (len(states) == 0):
        states = None

    # sort states
    else:
        states = sorted(states)

    return states

def get_state_type(state):
    """Get eigenstate type

    Args:
        state (complex): eigenstate

    Returns:
        str

    """     

    stype = 'bound' if (np.real(state) < 0.0) else 'resonance'

    return stype

def potential_maximum(potential, a, b, *params):
    """Get potential maximum

    Args:
        potential (func): potential function
        a (float): left boundary of interval
        b (float): right boundary of interval
        params (args): potential specific parameters

    Returns:
        float

    """     

    n = 500
    x = np.linspace(a, b, n)
    y = np.zeros(n, dtype=float)
    
    for i in range(n):
        y[i] = potential(x[i], *params)
        
    peaks, _ = find_peaks(y, height=0)
    maximum = y[peaks[0]]

    return maximum

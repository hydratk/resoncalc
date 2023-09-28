# -*- coding: utf-8 -*
"""Output

.. module:: output
   :synopsis: Module for generating output

"""

# imports
import numpy as np
import resoncalc.graphics as graphics
import resoncalc.detection as detection
from datetime import datetime
from os import getcwd, path, mkdir
from json import dumps

# globals
outdir = ''   # output directory
logfile = ''  # log filename
log_level = 1 # log level 0 - ERROR, 1 - INFO, 2 - DEBUG
stdout = True # write to standard output
outstates = ['bound', 'resonance'] # exported eigenstates
outfiles = ['eigenstates', 'eigenvalues', 'potential_grid', 'spectrum', 'log', 'test'] # generated output files

def set_outdir(dirname):
    """Set output directory

    Args:
        dirname (str): directory name

    """    

    globals()['outdir'] = path.join(getcwd() if (len(outdir) == 0) else outdir, '{0}_{1}'.format(dirname, datetime.now().strftime('%Y%m%d%H%M%S')))
    mkdir(outdir)
    globals()['logfile'] = path.join(outdir, 'log.log')
 
def log(msg, level):
    """Log message

    Args:
        msg (str): message
        level (int): message level

    """       

    # message levels
    levels = ['ERROR', 'INFO', 'DEBUG']

    # check level
    if (level <= log_level):
        msg_formatted = '{0} [{1}] {2}'.format(datetime.now(), levels[level], msg)

        # write to stdout
        if (stdout):
            print(msg_formatted)

        # output enabled
        if ('log' in outfiles):
            with open(logfile, 'a+') as f:
                f.write(msg_formatted + '\n')

def error(msg):
    """Log error message

    Args:
        msg (str): message

    """     

    log(msg, 0)                

def info(msg):
    """Log info message

    Args:
        msg (str): message

    """       

    log(msg, 1)

def debug(msg):
    """Log debug message

    Args:
        msg (str): message

    """     

    log(msg, 2)

def format_params(*params):
    """Format parameters to be used in filename

    Args:
        params (args): parameters

    """     

    params_formatted = '_'.join([str(round(param, 5)).replace('.', '') for param in params[:-1]])

    return params_formatted

def export_test(cfg):
    """Export test configuration to file

    Args:
        cfg (dict): test configuration

    """  

    # output enabled
    if ('test' in outfiles):

        fname = path.join(outdir, '{0}.json'.format(cfg['potential']))
        debug('Exporting test configuration')

        with open(fname, 'a+') as f:
            f.write(dumps(cfg, indent=2))

def export_eigenstates(params, eigenstates):
    """Export eigenstates to file

    Args:
        params (list): potential parameters
        states (list): states

    """  

    fname = path.join(outdir, 'eigenstates.csv')
    debug('Exporting eigenstates')
    header = 'param1,param2,l,real,imag,type\n'

    with open(fname, 'a+') as f:
        f.write(header)
        
        for i in range(len(eigenstates)):
            
            # state enabled
            state_type = detection.get_state_type(eigenstates[i])
            if (state_type in outstates):
                f.write('{0},{1},{2},{3:e},{4:e},{5}\n'.format(round(params[i][0], 5), round(params[i][1], 5), params[i][2], np.real(eigenstates[i]),
                                                               np.imag(eigenstates[i]), state_type))

def export_eigenvalues(eigenvalues, phase, *params):
    """Export eigenvalues to file

    Args:
        eigenvalues (list): eigenvalues
        phase (float): ECS phase
        params (args): potential parameters

    """  

    # output enabled
    if ('eigenvalues' in outfiles):
        fname = path.join(outdir, 'eigenvalues_{0}_{1}.csv'.format(format_params(*params), phase))
        debug('Exporting eigenvalues')
        header = 'real,imag\n'

        with open(fname, 'a+') as f:
            f.write(header)
        
            for val in eigenvalues:
                f.write('{0:e},{1:e}\n'.format(np.real(val), np.imag(val)))            

def export_potential_grid(points, potential, phase, *params):
    """Export potential grid to file

    Args:
        points (list): points
        potential (list): potential values
        phase (float): ECS phase
        params (args): potential parameters

    """  

    # output enabled
    if ('potential_grid' in outfiles):
        fname = path.join(outdir, 'potential_grid_{0}_{1}.csv'.format(format_params(*params), phase))
        debug('Exporting potential grid')
        header = 'point real,point imag,potential real,potential imag\n'
        points = points[1:-1]

        with open(fname, 'a+') as f:
            f.write(header)
        
            for i in range(len(points)):
                f.write('{0:e},{1:e},{2:e},{3:e}\n'.format(np.real(points[i]), np.imag(points[i]), np.real(potential[i]), np.imag(potential[i]))) 

def export_complex_spectrum_fig(eigenvalues1, eigenvalues2, states, *params):
    """Export complex spectrum to figure

    Args:
        eigenvalues (list): eigenvalues for first angle or rotation
        eigenvalues2 (list): eigenvalues for second angle of rotation
        states (list): highlighted eigenstates        
        params (args): potential parameters

    """  

    # output enabled
    if ('spectrum' in outfiles):
        fname = path.join(outdir, 'spectrum_{0}.png'.format(format_params(*params)))
        debug('Exporting complex spectrum figure')
        graphics.plot_complex_spectrum(eigenvalues1, eigenvalues2, states, fname)

def export_eigenstates_fig(potential, states, a, b, emax, *params):
    """Export eigenstates to figure

    Args:
        potential (func): potential function
        states (list): eigenstates
        a (float): left boundary of interval
        b (float): right boundary of interval
        emax (float): maximum energy in atomic units
        params (args): potential parameters

    """ 

    # output enabled
    if ('eigenstates' in outfiles):
        fname = path.join(outdir, 'eigenstates_{0}.png'.format(format_params(*params)))
        debug('Exporting eigenstates figure')
        graphics.plot_eigenstates(potential, states, a, b, emax, fname, *params)    

def generate_graphs(infile, title=None):
    """Generate graphs from input

    Args:
        infile (str): input filename
        title (str): title, default empty

    """ 

    title = '_{0}'.format(title) if (title is not None) else ''
    graphics.plot_resonances_complex(infile, path.join(outdir, 'resonances_complex{0}.png'.format(title)))
    graphics.plot_resonances_params(infile, True, path.join(outdir, 'resonances_energy{0}.png'.format(title)))
    graphics.plot_resonances_params(infile, False, path.join(outdir, 'resonances_width{0}.png'.format(title)))

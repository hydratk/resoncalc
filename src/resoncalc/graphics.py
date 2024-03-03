# -*- coding: utf-8 -*
"""Graphics

.. module:: graphics
   :synopsis: Module for creating graphs 

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import resoncalc.output as output
from math import ceil
from csv import DictReader

# globals
n_points = 1000 # count of points in graph
n_rows = 20.0   # count of rows in graph legend

def plot_fem_base_function(func, a, b):
    """Plot base function

    Args:
        func (func): base function
        a (float): left boundary of interval
        b (float): right boundary of interval

    """    

    # coordinates init 
    x = np.linspace(np.real(a), np.real(b), n_points)
    n = len(x)
    y = np.zeros(n, dtype=float)

    # calculate function values
    for i in range(n):
        y[i] = np.real(func(x[i]))

    # create graph
    plt.plot(x, y)
    plt.title('FEM base function')
    plt.tight_layout()
    plt.show()

def plot_fem_base(funcs, a, b):
    """Plot base functions

    Args:
        funcs (list): base functions
        a (float): left boundary of interval
        b (float): right boundary of interval

    """    

    # coordinate init
    x = np.linspace(np.real(a), np.real(b), n_points)
    n = len(x)    

    # function loop
    for i in range(len(funcs)):
        # coordinate init
        y = np.zeros(n, dtype=float)

        # bridging function
        if (type(funcs[i]) is list):
            for j in range(n):
                y[j] = np.real(funcs[i][0](x[j]) + funcs[i][1](x[j]))

        # function within element        
        else:
            for j in range(n):
                y[j] = np.real(funcs[i](x[j]))

        # plot function     
        plt.plot(x, y, label=i+2)

    # create graph
    plt.title('FEM base functions')
    plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0.0, ncol=ceil(len(funcs)/n_rows))
    plt.tight_layout()
    plt.show()

def plot_potential(potential, a, b, *params):
    """Plot potential

    Args:
        potential (func): potential function
        a (float): left boundary of interval
        b (float): right boundary of interval
        params (args): potential specific parameters

    """    

    # coordinates init
    x = np.linspace(np.real(a), np.real(b), n_points)
    n = len(x)
    y = np.zeros(n, dtype=float)

    # plot potential
    for i in range(n):
        y[i] = potential(x[i], *params)
    plt.plot(x, y)


    # create graph
    plt.title('Potential')
    plt.xlabel('r')
    plt.ylabel('E')
    plt.grid()
    plt.tight_layout()
    plt.show()

def plot_states(potential, states, a, b, emax, fname='', *params):
    """Plot states for given potential

    Args:
        potential (func): potential function
        eigenvalues (list): states
        a (float): left boundary of interval
        b (float): right boundary of interval
        emax (float): maximum energy in atomic units
        fname (str): export filename
        params (args): potential specific parameters

    """    

    # coordinates init
    states = [np.real(v['value']) for v in states]
    x = np.linspace(a, b, n_points)
    n = len(x)
    y = np.zeros(n, dtype=float)

    # plot potential
    for i in range(n):
        y[i] = potential(x[i], *params)
    plt.plot(x, y)

    # plot states
    for i in range(len(states)):
        z = np.zeros(n, dtype=float)
        val = states[i]
        cross = 0

        # line segment representing state
        for j in range(n):
            
            # intersection points with potential
            if (val > y[j]):
                if (cross == 0):
                    cross += 1
            elif (cross == 1):
                cross += 1

            # plot line
            if (cross == 1):
                z[j] = val
            else:
                z[j] = np.nan

        plt.plot(x, z, 'r', label='E{0} = {1:e}'.format(i+1, val))

    # create graph
    plt.xlim(left=a, right=b)
    plt.ylim(bottom=None, top=emax)
    plt.title('Potential')
    plt.xlabel('r')
    plt.ylabel('E')
    plt.grid()
    plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0.0, ncol=ceil(len(states)/n_rows))
    plt.tight_layout()

    # export graph
    if (len(fname) > 0):
        plt.savefig(fname)
        plt.close()
        
    # display graph
    else:
        plt.show()

def plot_complex_spectrum(eigenvalues, eigenvalues2=[], states=None, fname=''):
    """Plot complex spectrum used in ECS method

    Args:
        eigenvalues (list): eigenvalues for first angle or rotation
        eigenvalues2 (list): eigenvalues for second angle of rotation, default empty
        states (list): highlighted states
        fname (str): export filename

    """    
    
    # plot first eigenvalues    
    x = []
    y = []
    if (states is None):        
        x = np.real(eigenvalues)
        y = np.imag(eigenvalues)
    else:
        limit = np.abs(states[0]['value'])
        for val in eigenvalues:
            if (np.real(val) <= limit):
                x.append(np.real(val))
                y.append(np.imag(val))
        
    plt.plot(x, y, 'b.')

    # plot second eigenvalues
    x = []
    y = []
    if (len(eigenvalues2) > 0):
        if (states is None):
            x = np.real(eigenvalues2)
            y = np.imag(eigenvalues2)
        else:
            limit = np.abs(states[0]['value'])
            for val in eigenvalues2:
                if (np.real(val) <= limit):
                    x.append(np.real(val))
                    y.append(np.imag(val))

    plt.plot(x, y, 'g.')

    # highlight states
    if (states is not None):
        for state in states:
            value = state['value']
            plt.plot(np.real(value), np.imag(value), 'ro', label='real={0:e}, imag={1:e}'.format(np.real(value), np.imag(value)))
        plt.legend(loc='lower left')

    # create graph    
    plt.title('Complex spectrum')
    plt.xlabel('Re')
    plt.ylabel('Im')    
    plt.grid()    

    # export graph
    if (len(fname) > 0):
        plt.savefig(fname)
        plt.close()
        
    # display graph
    else:
        plt.show()

def plot_resonances_complex(infile, outfile=''):
    """Plot resonances in complex plane

    Args:
        infile (str): input filename
        outfile (str): output filename, default empty

    """   

    # parse csv file
    with open(infile, 'r') as f:
        reader = DictReader(f)
        states = {}

        # get resonance states, group by param2
        for row in reader:                        
            if (row['type'] == 'resonance'):
                param2 = row['param2']
                if (not param2 in states):
                    states[param2] = []
                    
                states[param2].append(float(row['real']) + 1j*float(row['imag']))

    # plot state trajectories
    for k, v in states.items():
        plt.plot(np.real(v), np.imag(v), 'b-')

    # create graph    
    plt.title('Resonances')
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.grid()

    # export graph
    if (len(outfile) > 0):
        plt.savefig(outfile)
        plt.close()
        
    # display graph
    else:
        plt.show()

def plot_resonances_params(infile, energy=True, outfile=''):
    """Plot resonances according to parameters

    Args:
        infile (str): input filename
        energy (bool): energy or width, default energy
        outfile (str): output filename, default empty

    """   

    # parse csv file
    with open(infile, 'r') as f:
        reader = DictReader(f)
        param1, param2, value = [], [], []

        # get resonance states
        for row in reader:                        
            if (row['type'] == 'resonance'):
                param1.append(float(row['param1']))
                param2.append(float(row['param2']))
                value.append(float(row['real']) if (energy) else -0.5 * float(row['imag']))

    # create graph
    plt.scatter(param1, param2, c=value, cmap=plt.cm.rainbow)    
    plt.colorbar()
    plt.title('Resonance {0}'.format('energy' if (energy) else 'width'))
    plt.xlabel('param a')
    plt.ylabel('param b')
    plt.grid()

    # export graph
    if (len(outfile) > 0):
        plt.savefig(outfile)
        plt.close()
        
    # display graph
    else:
        plt.show()        

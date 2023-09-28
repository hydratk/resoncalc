# -*- coding: utf-8 -*
"""FEM DVR

.. module:: femdvr
   :synopsis: Module for calculating Schrodinger equation spectrum using FEM-DVR method

"""

# imports
import numpy as np
import resoncalc.output as output
from sympy.integrals.quadrature import gauss_lobatto
from scipy.linalg import eigvals, eig
from functools import reduce
from operator import mul, add

def equidistant_endpoints(a, b, n):
    """Create equidistant endpoints for FEM method

    Args:
        a (float): left boundary of interval
        b (float): right boundary of interval
        n (int): count of elements

    Returns:
        list: endpoints 

    """

    # list init
    endpoints = np.zeros(n+1, dtype=np.complex_)

    # equidistant division
    d = (b-a)/n
    for i in range(n+1):
        endpoints[i] = a + i*d

    return endpoints

def progressive_endpoints(a, b, n, l):
    """Create progressive endpoints for FEM method

    Args:
        a (float): left boundary of interval
        b (float): right boundary of interval
        n (int): count of elements
        l (float): length of first element

    Returns:
        list: endpoints 

    """

    # list init
    endpoints = np.zeros(n+1, dtype=np.complex_)

    # progressive division
    d = 2 * (b - a - n*l) / (n*(n-1))
    endpoints[0] = a
    
    for i in range(1, n+1):
        endpoints[i] = endpoints[i-1] + l + (i-1)*d

    return endpoints

def merge_endpoints(lists):
    """Merge multiple endpoints

    Args:
        lists (list): endpoints

    Returns:
        list: endpoints 

    """

    endpoints = lists[0]
    for i in range(1, len(lists)):

        # common boundary stored once
        if (endpoints[-1] == lists[i][0]):
            endpoints = np.concatenate((endpoints, lists[i][1:]))
        else:
            endpoints = np.concatenate((endpoints, lists[i]))

    return endpoints

def exterior_complex_scaling(endpoints, x0, phase):
    """Apply coordinate transformation using ECS method

    Args:
        endpoints (list): endpoints coordinates
        x0 (float): centre of rotation
        phase (int): angle of rotation in degrees

    Returns:
        list: endpoints 

    """    
    
    for i in range(len(endpoints)):
        point = endpoints[i]

        # apply ECS
        if (point >= x0):
            point = x0 + (point-x0) * np.exp(1j*phase*np.pi/180)
        endpoints[i] = point

    return endpoints

def gauss_lobatto_points_weights(n, a=-1.0, b=1.0):
    """Calculate Gauss-Lobatto quadrature

    Args:
        n (int): order of quadrature 
        a (float): left boundary of interval, default -1
        b (float): right boundary of interval, default 1

    Returns:
        tuple: points (list), weights (list) 

    """    

    # quadrature calculation on [-1, 1] using Sympy 
    points, weights = gauss_lobatto(n, 16)

    # rescale to [a, b]
    c1 = (b-a)/2.0
    c2 = (b+a)/2.0
    
    for i in range(n):
        points[i] = c1 * points[i] + c2
        weights[i] = c1 * weights[i]
        
    return (points, weights)
 

def fem_points_weights(n, endpoints):
    """Calculate points and weights for all FEM elements

    Args:
        n (int): order of quadrature 
        endpoints (list): endpoints coordinates

    Returns:
        tuple: points (list), weights (list) 

    """    

    n_elems = len(endpoints) - 1   # count of elements
    n_points = n_elems * (n-1) + 1 # count of points

    # list init
    points = np.zeros(n_points, dtype=np.complex_)
    weights = np.zeros(n_points, dtype=np.complex_)

    # element loop
    for i in range(n_elems):

        # quadrature for given element
        el_points, el_weights = gauss_lobatto_points_weights(n, endpoints[i], endpoints[i+1])

        # store points and weights
        idx = i * (n-1)
        points[idx : idx + n] = el_points
        weights[idx : idx + n] = weights[idx : idx + n] + el_weights # bridging on element boundary

    return (points, weights)

def fem_base_function(idx, points, weight):
    """Calculate base function for FEM-DVR method

    Args:
        idx (int): function index
        points (list): element endpoints
        weight (float): weight

    Returns:
        func: base function

    """    

    n = len(points)

    # piecewise definition of function using product of polynomials
    func = (lambda x: 0 if (x < points[0] or x > points[n-1])
                        else 1.0/np.sqrt(weight) * reduce(mul, [(x - points[i]) / (points[idx] - points[i]) for i in range(n) if i != idx]))
    
    return func

def fem_base(n, endpoints):
    """Calculate function base for FEM-DVR method

    Args:
        n (int): order of quadrature
        endpoints (list): element endpoints

    Returns:
        list: base functions

    """    

    # quadrature for all elements
    points, weights = fem_points_weights(n, endpoints)
    funcs = []

    # element loop
    for i in range(len(endpoints)-1):
        idx = i * (n-1)

        # within element
        for j in range(n):

            # calculate base function
            func = fem_base_function(j, points[idx : idx + n], weights[idx + j])

            # bridging function composed of two base functions
            if (j == 0 and i > 0):
                funcs[-1] = [funcs[-1], func]
            else:
                funcs.append(func)

    # first and last functions not used to satisfy boundary conditions
    return funcs[1 : len(funcs)-1]

def lagrange_derivative_matrix(n):
    """Calculate matrix of derivatives of Lagrange polynomials used in stiffness matrix

    Args:
        n (int): order of quadrature

    Returns:
        2d array: matrix

    """    

    # matrix init
    matrix = np.zeros([n, n], dtype=np.complex_)

    # quadrature
    points, weights = gauss_lobatto_points_weights(n)

    # row loop
    for i in range(n):

        # column loop
        for j in range(n):

            # diagonal element     
            if (i == j):
                total = 0.0
                for k in [a for a in range(n) if a != i]:
                    total += 1.0/(points[i] - points[k])
                    
            # non-diagonal element
            else:
                total = 1.0/(points[i] - points[j])
                for k in [a for a in range(n) if (a != i and a != j)]:
                    total *= (points[j] - points[k]) / (points[i] - points[k])

            # store element
            matrix[i, j] = total

    return matrix

def stiffness_matrix(n, endpoints):
    """Calculate stiffness matrix

    Args:
        n (int): order of quadrature
        endpoints (list): endpoints coordinates

    Returns:
        2d array: matrix

    """    

    n_elems = len(endpoints) - 1   # count of elements
    n_points = n_elems * (n-1) + 1 # count of points
    overlap = 0j

    # matrix init
    matrix = np.zeros([n_points, n_points], dtype=np.complex_)

    # quadrature for all elements
    points, weights = fem_points_weights(n, endpoints)

    # matrix of derivative of lagrange polynomials
    lagrange_matrix = lagrange_derivative_matrix(n)

    # element loop
    for i in range(n_elems):

        # rescale matrix to given interval
        scaled = 2.0/(endpoints[i+1] - endpoints[i]) * lagrange_matrix

        # quadrature for given element
        el_points, el_weights = gauss_lobatto_points_weights(n, endpoints[i], endpoints[i+1])        
        idx = i * (n-1)

        # apply weight
        for j in range(n):
            scaled[j, :] = scaled[j, :] / np.sqrt(weights[idx + j])
            el_weights[j] = complex(el_weights[j])

        # within element
        for j in range(n):
            row = idx + j

            # matrix is symetric, only part is calculated
            for k in range(j+1):                
                col = idx + k
                total = 0j
                
                for l in range(n):                    
                    total -= el_weights[l] * scaled[j, l] * scaled[k, l]

                # store data
                matrix[row, col] = total
                matrix[col, row] = total                

        # apply overlap on element block matrices
        matrix[idx, idx] += overlap
        overlap = matrix[idx+n-1, idx+n-1]

    return matrix

def hamiltonian(n, endpoints, potential, points, mu=1, *params):
    """Calculate hamiltonian matrix

    Args:
        n (int): order of quadrature
        endpoints (list): endpoints coordinates
        potential (func): potential function
        points (list): points where potential is calculated
        mu (float): reduced mass in atomic units, default 1
        params (args): potential specific parameters

    Returns:
        2d array: matrix
        list: potential on grid

    """
    
    n_elems = len(endpoints) - 1   # count of elements
    n_points = n_elems * (n-1) + 1 # count of points
    h = 1                          # Planck's constant

    # kinetic term from stiffness matrix, boundary rows and column not used to satisfy boundary conditions
    matrix = -h**2/(2*mu) * stiffness_matrix(n, endpoints)[1 : n_points-1, 1 : n_points-1]    

    # potential term, operator is diagonal
    potential_grid = np.zeros(n_points-2, dtype=np.complex_)
    
    for i in range(0, n_points-2):
        pot = potential(points[i+1], *params)
        matrix[i, i] += pot
        potential_grid[i] = pot
    
    return (matrix, potential_grid)

def spectrum(hamiltonian):
    """Calculate spectrum of Hamiltonian matrix

    Args:
        hamiltonian (2d array): Hamiltonian matrix

    Returns:
        list: eigenvalues

    """    

    eigenvalues = []

    # calculate eigenvalues using Scipy
    vals = eigvals(hamiltonian)

    # sort by real part ascending
    vals = np.flip(np.sort_complex(vals))

    return vals

# -*- coding: utf-8 -*
"""Potential

.. module:: potential
   :synopsis: Module for potential functions definition

"""

# imports
import numpy as np

def centrifugal(x, l=1, mu=1.0):
    """Common centrifugal term

    Args:
        x (float): variable 
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1

    Returns:
        function

    """    

    potential = l*(l+1)/(2*mu*x**2)

    return potential

def morse(x, a=4.7446, b=1.440558, l=1.0, mu=1.0, c=0.7416):
    """Morse potential

    Args:
        x (float): variable
        a (float): parameter, default 4.7446
        b (float): parameter, default 1.440558
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1
        c (float): parameter, default 0.7416

    Returns:
        function

    """    

    potential = a*(1.0 - np.exp(-b * (x-c)))**2 - a + centrifugal(x, l, mu)

    return potential

def parabolic_gaussian(x, a=0.028, b=0.028, l=1.0, mu=1.0, c=1.0):
    """Parabolic Gaussian potential

    Args:
        x (float): variable
        a (float): parameter, default 0.028        
        b (float): parameter, default 0.028
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1
        c (float): parameter, default 1.0

    Returns:
        function

    """    

    potential = (a*x**2 - c) * np.exp(-b*x**2) + centrifugal(x, l, mu)
    
    return potential

def parabolic_gaussian2(x, a=0.028, b=0.028, l=1.0, mu=1.0, c=5.0):
    """Parabolic Gaussian potential 2

    Args:
        x (float): variable
        a (float): parameter, default 0.028        
        b (float): parameter, default 0.028
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1
        c (float): parameter, default 5.0

    Returns:
        function

    """    

    potential = a*(x**2 - c) * np.exp(-b*x**2) + centrifugal(x, l, mu)
    
    return potential

def gaussian(x, a, b, l=1.0, mu=1.0):
    """Gaussian potential

    Args:
        x (float): variable
        a (float): parameter
        b (float): parameter
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1

    Returns:
        function

    """    

    potential = a * np.exp(-b*x**2) + centrifugal(x, l, mu)
    
    return potential

def exponential(x, a, b, l=1.0, mu=1.0):
    """Exponential potential

    Args:
        x (float): variable
        a (float): parameter
        b (float): parameter
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1

    Returns:
        function

    """    

    potential = a * np.exp(-b*x) + centrifugal(x, l, mu)

    return potential

# map title to function
mapping = {
           'morse' : morse,
           'parabolic_gaussian' : parabolic_gaussian,
           'parabolic_gaussian2' : parabolic_gaussian2,
           'gaussian' : gaussian,
           'exponential' : exponential
          }

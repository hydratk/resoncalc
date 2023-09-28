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

def morse(x, De=4.7446, a=1.440558, l=1.0, mu=1.0, re=0.7416):
    """Morse potential

    Args:
        x (float): variable
        De (float): parameter, default 4.7446
        a (float): parameter, default 1.440558
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default 1
        r (float): parameter, default 0.7416

    Returns:
        function

    """    

    potential = De*(1.0 - np.exp(-a * (x-re)))**2 - De + centrifugal(x, l, mu)

    return potential

def harmonic_gaussian(x, a=0.028, c=0.028, l=1.0, mu=1.0, b=1.0):
    """Harmonic Gaussian potential

    Args:
        x (float): variable
        a (float): parameter, default 0.028        
        c (float): parameter, default 0.028
        l (int): secondary quantum number, default 1
        mu (float): reduced mass in atomic units, default
        b (float): parameter, default 1.01

    Returns:
        function

    """    

    potential = (a*x**2 - b) * np.exp(-c*x**2) + centrifugal(x, l, mu)
    
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
           'harmonic_gaussian' : harmonic_gaussian,
           'gaussian' : gaussian,
           'exponential' : exponential
          }

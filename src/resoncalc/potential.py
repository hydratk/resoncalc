# -*- coding: utf-8 -*
"""Potential

.. module:: potential
   :synopsis: Module for potential functions definition

"""

# imports
import numpy as np

def centrifugal(r, l=1, mu=1.0):
    """Common centrifugal term

    Args:
        r (float): variable 
        l (int): angular momentum, default 1
        mu (float): reduced mass in atomic units, default 1

    Returns:
        function

    """    

    potential = l*(l+1)/(2*mu*r**2)

    return potential

def morse(r, a=4.7446, b=1.440558, l=1.0, mu=1.0, c=0.7416):
    """Morse potential

    Args:
        r (float): variable
        a (float): parameter, default 4.7446
        b (float): parameter, default 1.440558
        l (int): angular momentum, default 1
        mu (float): reduced mass in atomic units, default 1
        c (float): parameter, default 0.7416

    Returns:
        function

    """    

    potential = a*(1.0 - np.exp(-b * (r-c)))**2 - a + centrifugal(r, l, mu)

    return potential

def parabolic_gaussian(r, a=0.028, b=0.028, l=1.0, mu=1.0, c=1.0):
    """Parabolic Gaussian potential

    Args:
        r (float): variable
        a (float): parameter, default 0.028        
        b (float): parameter, default 0.028
        l (int): angular momentum, default 1
        mu (float): reduced mass in atomic units, default 1
        c (float): parameter, default 1.0

    Returns:
        function

    """    

    potential = (a*r**2 - c) * np.exp(-b*r**2) + centrifugal(r, l, mu)
    
    return potential

def parabolic_gaussian2(r, a=0.028, b=0.028, l=1.0, mu=1.0, c=5.0):
    """Parabolic Gaussian potential 2

    Args:
        r (float): variable
        a (float): parameter, default 0.028        
        b (float): parameter, default 0.028
        l (int): angular momentum, default 1
        mu (float): reduced mass in atomic units, default 1
        c (float): parameter, default 5.0

    Returns:
        function

    """    

    potential = a*(r**2 - c) * np.exp(-b*r**2) + centrifugal(r, l, mu)
    
    return potential

def gaussian(r, a, b, l=1.0, mu=1.0):
    """Gaussian potential

    Args:
        r (float): variable
        a (float): parameter
        b (float): parameter
        l (int): angular momentum, default 1
        mu (float): reduced mass in atomic units, default 1

    Returns:
        function

    """    

    potential = a * np.exp(-b*r**2) + centrifugal(r, l, mu)
    
    return potential

def exponential(r, a, b, l=1.0, mu=1.0):
    """Exponential potential

    Args:
        r (float): variable
        a (float): parameter
        b (float): parameter
        l (int): angular momentum, default 1
        mu (float): reduced mass in atomic units, default 1

    Returns:
        function

    """    

    potential = a * np.exp(-b*r) + centrifugal(r, l, mu)

    return potential

# map title to function
mapping = {
           'morse' : morse,
           'parabolic_gaussian' : parabolic_gaussian,
           'parabolic_gaussian2' : parabolic_gaussian2,
           'gaussian' : gaussian,
           'exponential' : exponential
          }

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 06:30:00 2019

@author: Shibabrat Naik, shiba@vt.edu
"""

import numpy as np
from scipy.optimize import fsolve


def potential_energy(x, y, params_pe):
    """
    Potential energy function for the 2 DOF model
    x, y: N x N array of position values as meshgrid data
    params: 1 x 4 array of parameters
    """
    if np.size(x) == 1:
        nX = 1
        nY = 1
    else:
        nX = np.size(x,1)
        nY = np.size(x,0)
        x = np.ravel(x, order = 'C')
        y = np.ravel(y, order = 'C')
    
    mu, alpha, omega, epsilon  = params_pe
    
    vx = -np.sqrt(mu)*x**2 + (alpha/3.0)*x**3
    vy = 0.5*omega**2*y**2
    vxy = 0.5*epsilon*(y - x)**2
    
    pe = np.reshape(vx + vy + vxy, (nY, nX))
#     print(np.size(vy))
#     print(x,y,vyx)
    
    return pe


def vector_field(params, t, states):
    
    massA, massB, mu, alpha, omega, epsilon = params
    q, x, p, px = states
    
    q1Dot = p/massA
    q2Dot = px/massB
    p1Dot = -( -2*np.sqrt(mu)*q + alpha*q**2 - epsilon*(x - q) )
    p2Dot = -( omega**2*x + epsilon*x - epsilon*q )
    
    return np.array([q1Dot, q2Dot, p1Dot, p2Dot])

def cross_sos_center_eqpt(params, t, states):
    """
    Event function to catch trajectories crossing the SOS located at the bottom of the well, that is the center equilibrium point
    """

    massA, massB, mu, alpha, omega, epsilon = params
    x_eqpt = (1/alpha)*(2.0*np.sqrt(mu) \
        - (omega**2*epsilon)/(omega**2 + epsilon)) 

    # cross_sos_center_eqpt.terminal = True
    # cross_sos_center_eqpt.direction = 0

    return states[0] - x_eqpt

# 5 is arbitrary
# just to ensure the trajectory is not going to turn back
def cross_xneg5(params, t, states):
    return states[0] - (-5) 


# Obtain one of the momenta coordinate based on input values of other coordinates on the fixed energy surface
def momentum_fixed_energy(x, y, px, params, E):

    massA, massB, mu, alpha, omega, epsilon = params
    py = 0
    # print(params)

    potential_energy_val = potential_energy(x, y, params[2:]) 

    if (E >= (potential_energy_val + (1/(2.0*massA))*(px**2.0))):
        py = np.sqrt( 2.0*massB*(E - (potential_energy_val \
            + (1/(2.0*massA))*(px**2.0)) ) )
    else:
        py = np.NaN
        # print("Momentum isn't real!")
    
    return py 


def energysurface_intersect_sos(params, deltaEnergy, y_guess, numpts_boundary = 1000):
    
    massA, massB, mu, alpha, omega, epsilon = params 
    saddleEnergy = 0

    x_e = (1/alpha)*(2.0*np.sqrt(mu) \
        - (omega**2*epsilon)/(omega**2 + epsilon))

    totalEnergy = saddleEnergy + deltaEnergy
    
    H_xe = lambda y: totalEnergy - (-np.sqrt(mu)*x_e**2 + (alpha/3)*x_e**3 \
        + ((omega**2)/2)*y**2 + 0.5*epsilon*(x_e - y)**2) 
        
    # minimum and maximum of the energy boundary's y coordinate on the sos
    yMax = fsolve(H_xe,  y_guess)
    yMin = fsolve(H_xe,  -y_guess)
    # print(yMin, yMax)
    
    yGrid = np.linspace(yMin, yMax, numpts_boundary)
    
    pyGridPos = np.zeros(np.size(yGrid))
    pyGridNeg = np.zeros(np.size(yGrid))
    
    for i in range(np.size(yGrid)):
        if H_xe(yGrid[i]) >= 0:
            pyGridPos[i] = np.sqrt((2*massB)*(H_xe(yGrid[i])))
            
            pyGridNeg[i] = -np.sqrt((2*massB)*(H_xe(yGrid[i])))
            
        
    yGrid = np.append(yGrid, np.flip(yGrid,0))
    pyGrid = np.append(pyGridPos, np.flip(pyGridNeg,0))
    
#     print(np.shape(xGrid))
        
    return np.array([yGrid, pyGrid]).T







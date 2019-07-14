#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 06:30:00 2019

@author: Shibabrat Naik, shiba@vt.edu
"""
import numpy as np

def V_sn2dof(x, y, params_pe):
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


def vec_field_sn2dof(states, *params):
    
    massA, massB, mu, alpha, omega, epsilon = params
    q, x, p, px = states
    
    q1Dot = p/massA
    q2Dot = px/massB
    p1Dot = -( -2*np.sqrt(mu)*q + alpha*q**2 - epsilon*(x - q) )
    p2Dot = -( omega**2*x + epsilon*x - epsilon*q )
    
    return np.array([q1Dot, q2Dot, p1Dot, p2Dot])















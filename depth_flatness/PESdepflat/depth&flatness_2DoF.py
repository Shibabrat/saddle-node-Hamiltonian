# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 11:20:06 2019

@author: wl16298
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint,quad,trapz,solve_ivp
import math
from IPython.display import Image # for Notebook
from IPython.core.display import HTML
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
import scipy.linalg as linalg
from scipy.optimize import fsolve
import time
from functools import partial
#import saddlenode_tpcd ### import module xxx where xxx is the name of the python file xxx.py 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib import cm
from pylab import rcParams
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
%matplotlib

#%%
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
#alpha=1.0
mu=0.1
omega=1.0
epsilon=0
axis_fs=20

#%% Depth for the 2DoF (un)coupled system
def depth(alpha,mu,omega,epsilon):
    """ This function returns the value of depth for a given set of parameters
    """
    depth = -(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)
    return depth

alpha = np.linspace(0,5,1000)
epsilon = np.array([0,0.01,0.05,0.1,0.5,1.0,1.5])
ax = plt.gca()
# we want to use non zero alpha as depeth is not defined for alpha=0
plot1 = ax.plot(alpha[1:],depth(alpha[1:],mu,omega,epsilon[0]),label=r'$\mathcal{D},\epsilon=0$')
plot2 = ax.plot(alpha[1:],depth(alpha[1:],mu,omega,epsilon[2]),label=r'$\mathcal{D},\epsilon=0.05$')
plot3 = ax.plot(alpha[1:],depth(alpha[1:],mu,omega,epsilon[3]),label=r'$\mathcal{D},\epsilon=0.1$')
plot4 = ax.plot(alpha[1:],depth(alpha[1:],mu,omega,epsilon[4]),label=r'$\mathcal{D},\epsilon=0.5$')
plot5 = ax.plot(alpha[1:],depth(alpha[1:],mu,omega,epsilon[5]),label=r'$\mathcal{D},\epsilon=1.0$')
plot6 = ax.plot(alpha[1:],depth(alpha[1:],mu,omega,epsilon[6]),label=r'$\mathcal{D},\epsilon=1.5$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_xlim(0, 1)
ax.set_ylim(-5, 0)

plt.grid()
plt.show()

#%% Flatness for the 2DoF (un)coupled system
alpha = 1
mu=0.1
omega=1.0
epsilon=0.1
axis_fs=20
e=0.01

def dVdx(x,y,alpha,mu,epsilon):
    dVdx = alpha*x**2+ (epsilon-2*np.sqrt(mu))*x - epsilon*y
    return dVdx

def V_2dof(x,y,alpha,mu,omega,epsilon,V):
    """This function returns the function of the potential energy 
    """
#    f = (alpha/3)*x**3 +(0.5*epsilon-np.sqrt(mu))*x**2-epsilon*y*x +0.5*(omega**2+epsilon)*y**2-V
    f = (alpha/3)*x**3 -np.sqrt(mu)*x**2 +0.5*epsilon*x**2-V
    return f

def direct_flatness(y,alpha,mu,omega,epsilon):
    """This function returns the flatness in x, y direction for a given value of y. 
    Note that for the saddle node Hamiltonian, the flatness in y direction does not depend on x coordinate.
    However, the flatness in x direction does depend on y coordinate.
    x_sad: x coordinate of the saddle equilibrium pt.
    x_cen: x coordinate of the centre equilibrium pt.
    x_int: x coordinate of the pt with V=0, y=0.
    (We define this pt as the intersection pt between V=0 and the plane y=0 
      y-axis
        |  xxxx                  os = x_sad, 
        | x    x                  
        |os-----oi---->  x-axis  oi = x_int            
        |  x   x
        |   xxx                  ox..xox..x represents the projection of V=0 in x-y plane.
        |)
    x_guess is a guess value of x_int.
    flatness_x: flatness in x direction(requires a y value) . 
    flatness_y: flatness in y direction.
    """
    x_sad = 0.
    x_cen = (2*np.sqrt(mu) - epsilon + np.sqrt((epsilon- 2*np.sqrt(mu))**2 + 4*alpha*epsilon* y))/(2*alpha)
    x_guess = (2*np.sqrt(mu)/alpha- (omega**2*epsilon)/(alpha*(omega**2+epsilon)))*4/3
#    x_guess = x_cen+(1/3)*np.sqrt((epsilon- 2*np.sqrt(mu))**2 + 4*alpha*epsilon* y)/alpha
    x_int = fsolve(lambda x: V_2dof(x,y,alpha,mu,omega,epsilon,0),x_guess)
    flatness_x = 0.5*(0.25*np.sqrt((2*np.sqrt(mu)-epsilon)**2+4*epsilon*alpha*y) + abs(dVdx(x_int,y,alpha,mu,epsilon)/(x_int-x_cen)))
    flatness_y = omega**2+epsilon
    return flatness_x,flatness_y

#%% Compute flattness over a domain in x-y plane
num_alp=100 # number of alphas we want to calculate
F = np.zeros((num_alp,7))
epsilon = np.array([0,0.01,0.05,0.1,0.5,1.0,1.5])
alpha = np.linspace(0,5,num_alp)
for i in range(100):
    for j in range(7):
        y = np.linspace(-5,5,1000)  # defined interval of y coordinate
        flat_vec = np.zeros(1000)
        for k in range(1000):
            flat_x, flat_y= direct_flatness(y[k],alpha[i],mu,omega,epsilon[j])
            flat_vec[k] = np.sqrt(flat_x**2+flat_y**2)
        F[i,j] = np.nanmean(flat_vec)
        # We can calculate flatness values for a given set of parameters of a domain in x-y plane 
        # We then take the mean of the non nan values and define this number as the flatness over the domain in x-y plane.
        print(F[i,j])

flat_group = open("flatness_minalpha=%smaxalpha=%s_epsilon=%s.txt" %(alpha[0],alpha[-1],epsilon),'a+')
np.savetxt(flat_group.name,F,fmt='%1.16e')
flat_group.close()
#%% Plot flattness as a function of alpha
flat_group = open("flatness_minalpha=%smaxalpha=%s_epsilon=%s.txt" %(alpha[0],alpha[-1],epsilon),'a+')
F = np.loadtxt(flat_group.name)
flat_group.close()        
ax = plt.gca()
plot1 = ax.plot(alpha[1:],F[1:,0],label=r'$\mathcal{F},\epsilon=0$')
plot2 = ax.plot(alpha[1:],F[1:,1],label=r'$\mathcal{F},\epsilon=0.01$')
plot3 = ax.plot(alpha[1:],F[1:,2],label=r'$\mathcal{F},\epsilon=0.05$')
plot4 = ax.plot(alpha[1:],F[1:,3],label=r'$\mathcal{F},\epsilon=0.1$')
plot5 = ax.plot(alpha[1:],F[1:,4],label=r'$\mathcal{F},\epsilon=0.5$')
plot6 = ax.plot(alpha[1:],F[1:,5],label=r'$\mathcal{F},\epsilon=1.0$')
plot7 = ax.plot(alpha[1:],F[1:,6],label=r'$\mathcal{F},\epsilon=1.5$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 5)
ax.set_ylim(-0.01, 4)

plt.grid()
plt.show()
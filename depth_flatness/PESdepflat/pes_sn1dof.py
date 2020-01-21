# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 16:32:46 2019

@author: wl16298
"""

#%% plot of 1 DoF PES with alpha and mu, keep alpha as a constant

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pylab import rcParams
    
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

plt.style.use('seaborn') # use sans-serif fonts

mpl.rcParams['axes.spines.left'] = True   ## display axis spines
mpl.rcParams['axes.spines.bottom'] = True
mpl.rcParams['axes.spines.top'] = True
mpl.rcParams['axes.spines.right'] = True
# mpl.rcParams['xtick.top'] = True
# mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0

rcParams['figure.figsize'] = 8, 8 


# axis_fs=15

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal
#m = cm.ScalarMappable(cmap=cm.jet)
mpl.rcParams['font.weight'] = 'normal'

lw = 2.0 # linewidth for line plots
axis_fs=15


# Definition of the PES
def V(alpha,mu,x):
    """
    Potential energy function for the 1 dimensional saddle-node 
    Hamiltonian

    Parameters
    ----------
    alpha: float
    Parameter of the potential energy function

    mu: float
    Parameter of the potential energy function (must be > 0)

    x: float
    Position coordinate of the system

    Returns
    -------
    pot: float
    Potential energy value at the input position coordinate
    """

    pot = (1/3)*alpha*x**3 - math.sqrt(mu)*x**2
    return pot

# Depth definition
# def depth_sn1dof(alpha, mu):

#     # Get the saddle equilibrium point
#     x_saddle = 

#     # Get the center equilibrium point
#     x_center = 


#     depth = V(alpha, mu, x_saddle) - V(alpha, mu, x_center)

#     return depth


x = np.linspace(-5, 10, 1000)

#% plot of 1 DoF PES with alpha and mu, keep alpha as a constant
fig = plt.figure()
ax = plt.gca()

plot1 = ax.plot(x,V(1,0.5,x), linewidth = lw, \
                label=r'$\alpha=1,\mu=0.5$')

# plot2 = ax.plot(x,V(1,0.5,x),label=r'$\alpha=1,\mu=0.5$')
 
# plot3 = ax.plot(x,V(1,1.0,x),label=r'$\alpha=1,\mu=1.0$')

plot4 = ax.plot(x,V(1,2.0,x), linewidth = lw, \
                label=r'$\alpha=1,\mu=2.0$')
 
ax.scatter(0, 0, s = 100, c = 'r', marker = 'X')
#ax.scatter(0, 1, s = 100, c = 'r', marker = 'X')
ax.set_xlabel(r'$x$', fontsize=fal)
ax.set_ylabel(r'$V(x)$', fontsize=fal)

legend = ax.legend(loc='best', fontsize = ftl)

ax.set_xlim(-5, 15)
ax.set_ylim(-10, 5)
plt.grid()
plt.draw()
plt.pause(0.01)


#% plot of 1 DoF PES with alpha and mu, keep mu as a constant

ax = plt.gca()

# plot1 = ax.plot(x,V(0.5,1,x), linewidth = 1.0, \
                # label=r'$\alpha=0.5,\mu=1$')

plot2 = ax.plot(x,V(0.5,1.0,x), linewidth = lw, \
                label=r'$\alpha=0.5,\mu=1$')
 
# plot3 = ax.plot(x,V(1.5,1.0,x),label=r'$\alpha=1.5,\mu=1$')

plot4 = ax.plot(x,V(2.0,1.0,x), linewidth = lw, \
                label=r'$\alpha=2.0,\mu=1$')
 
ax.scatter(0, 0, s = 200, c = 'r', marker = 'X')
#ax.scatter(0, 1, s = 100, c = 'r', marker = 'X')
ax.set_xlabel('$x$', fontsize=fal)
ax.set_ylabel('$V(x)$', fontsize=fal)

legend = ax.legend(loc='best', fontsize = ftl)

ax.set_xlim(-5, 10)
ax.set_ylim(-10, 5)

plt.grid()
plt.draw()
# plt.pause(0.01)
plt.savefig('temp.pdf', bbox_inches = 'tight')
plt.show() #block = False



# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:36:40 2019

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
#mu=alpha**(3/4)
mu=0.1
omega=1.0
epsilon=0
axis_fs=20
#parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega]);
#eqNum = 1;  
#eqPt = saddlenode_tpcd.get_eq_pts_saddlenode(eqNum, parameters)

#depth = -(2*math.sqrt(parameters[4])- (parameters[6]**2*parameters[5])/(parameters[6]**2+parameters[5]))**3/(6*parameters[3]**2)
#%%
def depth(alpha,mu,omega,epsilon):
    depth = -(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)
    return depth
#%% definition of flatness(1st attempt)
alpha = np.linspace(1.e-5, 5, 1000)

def flatness(mu,omega):
    flatness = math.sqrt((7/4)**2*mu+ omega**4) 
    return flatness
ax = plt.gca()
plot1 = ax.plot(alpha,depth(alpha,mu,omega,epsilon),label=r'$\mathcal{D}, \mu=0.1,\omega=1$') # plot depth as a function of alpha
plot2 = ax.plot(alpha,np.ones(1000)*flatness(mu,omega),label=r'$\mathcal{F}, \mu=0.1,\omega=1$') # plot depth as a function of alpha


ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F},\mathcal{D}$', fontsize=axis_fs)

legend = ax.legend(loc='best')
ax.set_xlim(0, 2)
ax.set_ylim(-4, 2)
plt.grid()
plt.show()

#%%
plt.close('all')
# ratio of width of the bottleneck/width of the well
#widratio = lambda alpha: math.sqrt(e/(e+(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)))
alpha = np.linspace(1.e-5, 10, 100)
def widthratio(alpha,mu,omega,epsilon,e):
    ratio = np.sqrt(e/(e+(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)))
    return ratio

ax = plt.gca()
e= 0.01
plot1 = ax.plot(alpha,widthratio(alpha,mu,omega,epsilon,e),label=r'$R_{bw},\mu=0.1,\omega=1,H=0.01$') # plot width ratio as a function of alpha
e= 0.05
plot2 = ax.plot(alpha,widthratio(alpha,mu,omega,epsilon,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=0.05$') # plot width ratio as a function of alpha
e= 0.1
plot3 = ax.plot(alpha,widthratio(alpha,mu,omega,epsilon,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=0.1$') # plot width ratio as a function of alpha
e= 0.5
plot4 = ax.plot(alpha,widthratio(alpha,mu,omega,epsilon,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=0.5$') # plot width ratio as a function of alpha
e=1.0
plot6 = ax.plot(alpha,widthratio(alpha,mu,omega,epsilon,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=1.0$') # plot width ratio as a function of alpha


ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$R_{bw}$', fontsize=axis_fs)

legend = ax.legend(loc='best')
ax.set_xlim(0, 10)
ax.set_ylim(0, 1)
plt.grid()
plt.show()

#%% plot the ratio as a function of depth
plt.close('all')
# ratio of width of the bottleneck/width of the well
#widratio = lambda alpha: math.sqrt(e/(e+(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)))
depth = np.linspace(-10, -1.e-5, 100)
def widthratio(depth,e):
    ratio = np.sqrt(e/(e-depth))
    return ratio

ax = plt.gca()
e= 0.01
plot1 = ax.plot(depth,widthratio(depth,e),label=r'$R_{bw},\mu=0.1,\omega=1,H=0.01$') # plot width ratio as a function of alpha
e= 0.05
plot2 = ax.plot(depth,widthratio(depth,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=0.05$') # plot width ratio as a function of alpha
e= 0.1
plot3 = ax.plot(depth,widthratio(depth,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=0.1$') # plot width ratio as a function of alpha
e= 0.5
plot4 = ax.plot(depth,widthratio(depth,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=0.5$') # plot width ratio as a function of alpha
e=1.0
plot6 = ax.plot(depth,widthratio(depth,e),label=r'$R_{bw}, \mu=0.1,\omega=1,H=1.0$') # plot width ratio as a function of alpha



ax.set_xlabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_ylabel(r'$R_{bw}(\mathcal{D})$', fontsize=axis_fs)

legend = ax.legend(loc='best')
ax.set_xlim(-10, 0)
ax.set_ylim(0, 1)
plt.grid()
plt.show()
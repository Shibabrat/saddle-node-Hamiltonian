# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 15:46:40 2019

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

#%% parameters

resX=100

alpha = 1
mu=4
axis_fs=20

def depth(alpha, mu):
     """ definition of depth for 1dof system
     
     Parameters: 
     alpha
     mu
     
     Returns:
     float:
     depth
     """
     depth = -4*mu*(3/2)/(3*alpha**2)
     return depth


# definition of flatness (2nd version, i.e. mean of norm(dVdx) over some domain \Omega)
#%% 1dof flatness calculation
def flatness(a,b, alpha, mu):
    """This function returns the flatness for 1 dof system. flatness is defined as 
    the mean of norm(dVdx) over some domain \Omega \in [a, b].
    
    Parameters:
    a is left end of the domain \Omega
    b is right end of the domain \Omega
    alpha
    mu
    
    Returns:
    float:
        flatness
    """     
    flatness = 1/(b-a) * (-np.sqrt(mu)*(b**2-a**2)+ alpha/3*(b**3-a**3))
    
    return flatness


#%%1dof depth, flatness plottings
ax = plt.gca()
num_alp=100
alpha = np.linspace(0,10,num_alp)
F = flatness(-1,10,alpha, mu)
D = depth(alpha, mu)
plot1 = ax.plot(alpha[1:],abs(F[1:]),label=r'$\mathcal{F},\mu = 4$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-20, 300)

plt.grid()
plt.show()
#%%
plt.close('all')
ax = plt.gca()
plot1 = ax.plot(alpha[1:],D[1:],label=r'$\mathcal{D},\mu = 4$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-5, 0)

plt.grid()
plt.show()
#%% flatness aginist Rbw plottings
plt.close('all')

def widthratio(alpha,mu,e):
    """ This function returns the Rbw = wb/ww for 1 dof system .
    
    Parameters:
    alpha
    mu
    e is the total energy of the system
    
    Returns:
    float:
        Rbw
    
    """
    Rbw = np.sqrt(e/(e+4*mu**1.5/(3*alpha**2)))
    return Rbw

ax = plt.gca()
e= 0.005
plot1 = ax.plot(F[:],widthratio(alpha,mu,e),label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.005$') 
e= 0.01
plot2 = ax.plot(F[:],widthratio(alpha,mu,e),label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.01$') 
e= 0.05
plot3 = ax.plot(F[:],widthratio(alpha,mu,e),label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.05$') 
e= 0.1
plot4 = ax.plot(F[:],widthratio(alpha,mu,e),label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.1$') 
e=0.5
plot6 = ax.plot(F[:],widthratio(alpha,mu,e),label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.5$') 


ax.set_xlabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_ylabel('$R_{bw}(\mathcal{F})$', fontsize=axis_fs)

legend = ax.legend(loc='best')
ax.set_xlim(0, 300)
ax.set_ylim(0, 1)
plt.grid()
plt.show()
#%% depth aginist Rbw plottings
plt.close('all')

ax = plt.gca()
e= 0.005
plot1 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:],label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.005$') 
e= 0.01
plot2 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:],label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.01$') 
e= 0.05
plot3 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:],label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.05$') 
e= 0.1
plot4 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:],label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.1$') 
e=0.5
plot6 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:],label=r'$R_{bw}$ for the 1 dof system, $\mu=4,e=0.5$') 


ax.set_xlabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_ylabel('$R_{bw}(\mathcal{D})$', fontsize=axis_fs)

legend = ax.legend(loc='best')
ax.set_xlim(-20, 0)
ax.set_ylim(0, 1)
plt.grid()
plt.show()
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

import seaborn as sns
fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal
#m = cm.ScalarMappable(cmap=cm.jet)
mpl.rcParams['font.weight'] = 'normal'

lw = 3.0 # linewidth for line plots
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

plt.style.use('seaborn')

rcParams['figure.figsize'] = 8, 8 
axis_fs = 30
savefig_flag = True

#%% parameters

resX = 100
alpha = 1
mu = 4


def depth(alpha, mu):
     """ definition of depth for 1dof system
     
     Parameters: 
        alpha
        
        mu
     
     Returns:
        depth : float

     """
     depth = -4*mu*(3/2)/(3*alpha**2)
     return depth


# definition of flatness (2nd version, i.e. mean of norm(dVdx) over some domain \Omega), 1dof flatness calculation
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


#%%1dof depth, flatness plottings

sns.set(font_scale = 2)
figH = plt.figure()
ax = figH.gca()
num_alp=100
alpha = np.linspace(1e-10,10,num_alp)
F = flatness(-1,10,alpha, mu)
D = depth(alpha, mu)
plot1 = ax.plot(alpha[1:],abs(F[1:]),label=r'$\mathcal{F},\mu = 4$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-20, 300)

plt.grid('on')
plt.show()




#%%

plt.close('all')
sns.set(font_scale = 2)
figH = plt.figure()
ax = figH.gca()
plot1 = ax.plot(alpha[1:],D[1:],label=r'$\mathcal{D},\mu = 4$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-5, 0)

plt.grid('on')
plt.show()

#%% flatness against Rbw plottings
plt.close('all')

sns.set(font_scale = 2)
figH = plt.figure()
ax = figH.gca()

#$R_{bw}$ for the 1 dof system, $\mu=4,
e= 0.005
plot1 = ax.plot(abs(F[:]),widthratio(alpha,mu,e), lw=lw, label=r'$e=0.005$') 
e= 0.01
plot2 = ax.plot(abs(F[:]),widthratio(alpha,mu,e), lw=lw, label=r'$e=0.01$') 
e= 0.05
plot3 = ax.plot(abs(F[:]),widthratio(alpha,mu,e), lw=lw, label=r'$e=0.05$') 
e= 0.1
plot4 = ax.plot(abs(F[:]),widthratio(alpha,mu,e), lw=lw, label=r'$e=0.1$') 
e=0.5
plot6 = ax.plot(abs(F[:]),widthratio(alpha,mu,e), lw=lw, label=r'$e=0.5$') 


ax.set_xlabel(r'Flatness, $\mathcal{F}$', fontsize=axis_fs)
# ax.set_ylabel('Ratio of bottleneck to well-width, $R_{bw}$', fontsize=axis_fs) # (\mathcal{F})

legend = ax.legend(loc='best')
ax.set_xlim(0, 300)
ax.set_ylim(0, 1)
plt.grid('on')
plt.show()

if savefig_flag:
    figH.savefig('Rbw_1dof_flat_mu=4.pdf', bbox_inches = 'tight')


#%% depth against Rbw plottings
plt.close('all')

sns.set(font_scale = 2)
figH = plt.figure()
ax = figH.gca()

e= 0.005
plot1 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:], lw=lw, label=r'$e=0.005$') 
e= 0.01
plot2 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:], lw=lw, label=r'$e=0.01$') 
e= 0.05
plot3 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:], lw=lw, label=r'$e=0.05$') 
e= 0.1
plot4 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:], lw=lw, label=r'$e=0.1$') 
e=0.5
plot6 = ax.plot(D[1:],widthratio(alpha,mu,e)[1:], lw=lw, label=r'$e=0.5$') 


ax.set_xlabel(r'Depth, $\mathcal{D}$', fontsize = axis_fs)
ax.set_ylabel('Ratio of bottleneck to well-width, $R_{bw}$', fontsize = axis_fs) #(\mathcal{D})

legend = ax.legend(loc='best')
ax.set_xlim(-20, 0)
ax.set_ylim(0, 1)
plt.grid('on')
plt.show()


if savefig_flag:
    figH.savefig('Rbw_1dof_dep_mu=4.pdf', bbox_inches = 'tight')




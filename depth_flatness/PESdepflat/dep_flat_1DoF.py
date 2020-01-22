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
#%%
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=1.0
mu=4.0
axis_fs=30
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu])
#%% Depth for the 1DoF system
def depth_1dof_sn(par):
    """ This function returns the value of depth for a given set of parameters
    """
    depth = -4*(math.sqrt(par[4]))**3/(3*par[3]**2)
    return depth

def flatness_1dof_sn(x_min,x_max,num_pts,par):
    """Returns the value of flatness for a given domain [x_min,x_max] where we discretise this domain
    
       with num_pts number of points. definition of flatness is defined as the mean(norm) of nonnan values
       
       of dV/dx over some domain in x coordinates.
       
       
        Parameters
    
        ----------
    
        x_min : float
    
            = min value of the boundary of the domain we want to define our flatness on 
    
        x_max : float
    
            = max value of the boundary of the domain we want to define our flatness on 
            
        num_pts : int
            
            = number of points we want to discretise our domain with this number
              
              of points
    
        parameters : float (list)
    
            model parameters
    
    
    
        Returns
    
        -------
    
        F : float
    
            value of flatness over this particular domain with model parameters
        
    """
    def grad_pot_saddlenode(x, par):
        """This function returns the gradient of the potential energy function V(x,y)
        """     
    
        dVdx = par[3]*x**2-2*np.sqrt(par[4])*x
        
        return abs(dVdx)
    x = np.linspace(x_min,x_max,num_pts)
    normF = np.zeros((num_pts))
    for i in range(num_pts):
        normF[i] = grad_pot_saddlenode(x[i],par)
    F = np.nanmean(normF)
    return F

#%% 1dof depth plottings
sns.set(font_scale = 2)
num_alp = 1000
alpha = np.linspace(0,10,num_alp)
ax = plt.gca()
# we want to use non zero alpha as depeth is not defined for alpha=0
plot1 = ax.plot(alpha[1:],depth_1dof_sn([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu]),label=r'$\mathcal{D},\mu =%s$' %(mu))
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-5, 0)

#plt.grid()
plt.show()
plt.savefig('dep_mu=4_1dof.pdf', format='pdf', \

            bbox_inches='tight')

#%% 1dof flatness calculation
num_alp = 1000
alpha = np.linspace(0,10,num_alp)
F = np.zeros(num_alp)
for i in range(num_alp):
    parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i], mu])
    F[i] = flatness_1dof_sn(-1,10,500,parameters)
    print(F[i])


#%% 1dof flatness plottings
ax = plt.gca()
plot1 = ax.plot(alpha[1:],F[1:],label=r'$\mathcal{F},\mu=%s$' %(mu))

legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-0.01, 400)

#plt.grid()
plt.show()
#plt.savefig('flat_2ndtry_mu=4_1dof.pdf', format='pdf', \
#
#            bbox_inches='tight')
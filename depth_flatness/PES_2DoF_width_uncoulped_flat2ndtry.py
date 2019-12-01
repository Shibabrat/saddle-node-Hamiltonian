# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 16:43:25 2019

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
resX=100

alpha = 1
mu=4
omega=3.0
epsilon=0.0
axis_fs=20
e=0.01
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega])
#eqNum = 1;  
#eqPt = saddlenode_tpcd.get_eq_pts_saddlenode(eqNum, parameters)

#depth = -(2*math.sqrt(parameters[4])- (parameters[6]**2*parameters[5])/(parameters[6]**2+parameters[5]))**3/(6*parameters[3]**2)

def dVdx(x,y,alpha,mu,epsilon):
    dVdx = alpha*x**2+ (epsilon-2*np.sqrt(mu))*x - epsilon*y
    return dVdx

def V_2dof(x,y,alpha,mu,omega,epsilon,V):
#    f = (alpha/3)*x**3 +(0.5*epsilon-np.sqrt(mu))*x**2-epsilon*y*x +0.5*(omega**2+epsilon)*y**2-V
    f = (alpha/3)*x**3 -np.sqrt(mu)*x**2 +0.5*epsilon*x**2-V
    return f


#%% definition of flatness (2nd version, i.e. mean of norm(dVdx) over some domain \Omega)
#%% 2dof flatness calculation
def grad_pot_saddlenode(x, par):
    """This function returns the gradient of the potential energy function V(x,y)
    """     

    dVdx = par[3]*x[0]**2-2*np.sqrt(par[4])*x[0]+par[5]*(x[0]-x[1])
    dVdy = par[6]**2*x[1]-par[5]*(x[0]-x[1])
    
    normF = np.sqrt(dVdx**2 + dVdy**2)
    
    return normF

x = np.linspace(-1,10,500)
y = np.linspace(-5,5,500)
normF = np.zeros((500,500))
for i in range(500):
    for j in range(500):
        normF[i,j] = grad_pot_saddlenode([x[i],y[j]], parameters)

num_alp=100 # number of alphas we want to calculate
F = np.zeros((num_alp,8))
epsilon = np.array([0,1,2,3,4,5,6,7])
#epsilon = np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75])
alpha = np.linspace(0,10,num_alp)
#%% Perform the calculation
for i in range(100):
    for j in range(8):
        parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i], mu, epsilon[j], omega])
        normF = np.zeros((500,500))
        for k1 in range(500):
            for k2 in range(500):
                normF[k1,k2] = grad_pot_saddlenode([x[k1],y[k2]], parameters)
            
        F[i,j] = np.nanmean(normF)
        # We can calculate flatness values for a given set of parameters of a domain in x-y plane 
        # We then take the mean of the non nan values and define this number as the flatness over the domain in x-y plane.
        print(F[i,j])
        
flat_group = open("flatness_2ndtry_minalpha=%smaxalpha=%s_epsilon=%s.txt" %(alpha[0],alpha[-1],epsilon),'a+')
np.savetxt(flat_group.name,F,fmt='%1.16e')
flat_group.close()
#%% load data

flat_group = open("flatness_2ndtry_minalpha=%smaxalpha=%s_epsilon=%s.txt" %(alpha[0],alpha[-1],epsilon),'a+')
F = np.loadtxt(flat_group.name)
flat_group.close()

#%%2dof flatness plottings
ax = plt.gca()
plot1 = ax.plot(alpha[1:],F[1:,0],label=r'$\mathcal{F},\epsilon=0$')
#plot2 = ax.plot(alpha[1:],F[1:,1],label=r'$\mathcal{F},\epsilon=1$')
plot3 = ax.plot(alpha[1:],F[1:,2],label=r'$\mathcal{F},\epsilon=2$')
#plot4 = ax.plot(alpha[1:],F[1:,3],label=r'$\mathcal{F},\epsilon=3$')
plot5 = ax.plot(alpha[1:],F[1:,4],label=r'$\mathcal{F},\epsilon=4$')
#plot6 = ax.plot(alpha[1:],F[1:,5],label=r'$\mathcal{F},\epsilon=5$')
plot7 = ax.plot(alpha[1:],F[1:,6],label=r'$\mathcal{F},\epsilon=6$')
#plot8 = ax.plot(alpha[1:],F[1:,7],label=r'$\mathcal{F},\epsilon=7$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-0.01, 350)

plt.grid()
plt.show()


#%% flatness aginist widthratio plottings
plt.close('all')
# ratio of width of the bottleneck/width of the well
#widratio = lambda alpha: math.sqrt(e/(e+(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)))
alpha = np.linspace(0,10,num_alp)
epsilon = 0
def saddle_pt(alpha,mu,omega,epsilon):
    xe = 2*math.sqrt(parameters[4])/parameters[3] - (parameters[6]**2*parameters[5])/(parameters[3]*(parameters[6]**2+parameters[5]))
    ye = xe*parameters[5]/(parameters[6]**2+parameters[5])
    return xe,ye

def widthratio(alpha,mu,omega,epsilon,e):
    ratio = np.sqrt(e/(e+(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)))
    return ratio

ax = plt.gca()
e= 0.005
plot1 = ax.plot(F[:,0],widthratio(alpha,mu,omega,epsilon,e),label=r'width ratio of the uncoupled system, $\mu=4,\omega=3,e=0.005$') # plot width ratio as a function of alpha
e= 0.01
plot2 = ax.plot(F[:,0],widthratio(alpha,mu,omega,epsilon,e),label=r'width ratio of the uncoupled system, $\mu=4,\omega=3,e=0.01$') # plot width ratio as a function of alpha
e= 0.05
plot3 = ax.plot(F[:,0],widthratio(alpha,mu,omega,epsilon,e),label=r'width ratio of the uncoupled system, $\mu=4,\omega=3,e=0.05$') # plot width ratio as a function of alpha
e= 0.1
plot4 = ax.plot(F[:,0],widthratio(alpha,mu,omega,epsilon,e),label=r'width ratio of the uncoupled system, $\mu=4,\omega=3,e=0.1$') # plot width ratio as a function of alpha
e=0.5
plot6 = ax.plot(F[:,0],widthratio(alpha,mu,omega,epsilon,e),label=r'width ratio of the uncoupled system, $\mu=4,\omega=3,e=0.5$') # plot width ratio as a function of alpha


ax.set_xlabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_ylabel('$R_{bw}(\mathcal{F})$', fontsize=axis_fs)

legend = ax.legend(loc='best')
ax.set_xlim(0, 300)
ax.set_ylim(0, 1)
plt.grid()
plt.show()
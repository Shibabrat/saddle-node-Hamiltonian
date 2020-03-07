# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:12:09 2019

@author: wl16298
"""

"""We define the flatness as the mean || dVdx,dVdy || over some domain 
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
from scipy import optimize
from functools import partial

from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
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
resX=100
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=5
#mu=alpha**(3/4)
mu=4
omega=3.0
epsilon=0.05
axis_fs=30
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega])

#%% 2dof depth & flatness calculation
def depth_2dof_sn(par):
    """ This function returns the value of depth for a given set of parameters
    """
    depth = (2*math.sqrt(par[4])- (par[6]**2*par[5])/(par[6]**2+par[5]))**3/(6*par[3]**2)
    return depth

def flatness_2dof_sn(x_min,x_max,y_min,y_max,num_pts,par):
    """Returns the value of flatness for a given domain [x_min,x_max] x [y_min, y_max] where we discretise this domain
    
       with num_pts number of points. definition of flatness is defined as the mean(norm) of nonnan values
       
       of ||(dV/dx,dV/dy)||= \sqrt((dV/dx)**2+(dV/dy)**2) over some domain in x,y plane.
       
       
        Parameters
    
        ----------
    
        x_min : float
    
            = min x value of the boundary of the domain we want to define our flatness on 
    
        x_max : float
    
            = max x value of the boundary of the domain we want to define our flatness on 

        y_min : float
    
            = min y value of the boundary of the domain we want to define our flatness on 
    
        y_max : float
    
            = max y value of the boundary of the domain we want to define our flatness on 
            
        num_pts : int
            
            = number of points we want to discretise our domain with this number
              
              of points
    
        par : float (list)
    
            model parameters
    
    
    
        Returns
    
        -------
    
        F : float
    
            value of flatness over this particular domain with model parameters
        
    """
    def grad_pot_saddlenode2(x, par):
        """This function returns the gradient of the potential energy function V(x,y)
        """     
    
        dVdx = par[3]*x[0]**2-2*np.sqrt(par[4])*x[0]+par[5]*(x[0]-x[1])
        dVdy = par[6]**2*x[1]-par[5]*(x[0]-x[1])
        
        normF = np.sqrt(dVdx**2 + dVdy**2)
        
        return normF
    x = np.linspace(x_min,x_max,num_pts)
    y = np.linspace(y_min,y_max,num_pts)
    normF = np.zeros((num_pts,num_pts))
    for i in range(num_pts):
        for j in range(num_pts):
            normF[i,j] = grad_pot_saddlenode2([x[i],y[j]],par)
    F = np.nanmean(normF)
    return F

#%%Depth as a function of alpha
sns.set(font_scale = 2)
num_alp=1000
alpha = np.linspace(0,10,num_alp)
epsilon = np.array([0,1,2,3,4,5,6,7])
ax = plt.gca()
# we want to use non zero alpha as depeth is not defined for alpha=0
plot1 = ax.plot(alpha[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu,omega,epsilon[0]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[0]))
plot2 = ax.plot(alpha[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu,omega,epsilon[2]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[2]))
#plot3 = ax.plot(alpha[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu,omega,epsilon[3]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[3]))
plot4 = ax.plot(alpha[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu,omega,epsilon[4]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[4]))
#plot5 = ax.plot(alpha[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu,omega,epsilon[5]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[5]))
plot6 = ax.plot(alpha[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha[1:],mu,omega,epsilon[6]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[6]))
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
#plt.grid()
plt.show()
plt.savefig('dep_mu=4_omega=3_2dof.pdf', format='pdf', \

            bbox_inches='tight')

#%%Depth as a function of omega
sns.set(font_scale = 2)
num_ome=1000
alpha = 1
omega = np.linspace(0,10,num_ome)
epsilon = np.array([0,1,2,3,4,5,6,7])
ax = plt.gca()
# we want to use non zero alpha as depeth is not defined for alpha=0
plot1 = ax.plot(omega[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha,mu,omega[1:],epsilon[0]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[0]))
plot2 = ax.plot(omega[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha,mu,omega[1:],epsilon[2]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[2]))
#plot3 = ax.plot(omega[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha,mu,omega[1:],epsilon[3]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[3]))
plot4 = ax.plot(omega[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha,mu,omega[1:],epsilon[4]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[4]))
#plot5 = ax.plot(omega[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha,mu,omega[1:],epsilon[5]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[5]))
plot6 = ax.plot(omega[1:],depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S,alpha,mu,omega[1:],epsilon[6]])),lw=lw,label=r'$\mathcal{D},\epsilon=%s$' %(epsilon[6]))
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\omega$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(0, 12)
#plt.grid()
plt.show()
plt.savefig('dep_mu=4_alpha=1_2dof.pdf', format='pdf', \

            bbox_inches='tight')

#%%
num_alp=100 # number of alphas we want to calculate
epsilon = np.array([0,1,2,3,4,5,6,7])
#epsilon = np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75])
F = np.zeros((num_alp,len(epsilon)))
alpha = np.linspace(0,10,num_alp)

#%% Perform the calculation
for i in range(num_alp):
    for j in range(len(epsilon)):
        parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i], mu, epsilon[j], omega])
        
        F[i,j] = flatness_2dof_sn(-1,10,-5,5,500,parameters)
        print(F[i,j])
        
with open("flatness_2ndtry_minalpha=%smaxalpha=%s_epsilon=%s.txt" %(alpha[0],alpha[-1],epsilon),'a+') as flat_group:
    np.savetxt(flat_group.name,F,fmt='%1.16e')
#%% 2dof flatness plottings
sns.set(font_scale = 2)
with open("flatness_2ndtry_minalpha=%smaxalpha=%s_epsilon=%s.txt" %(alpha[0],alpha[-1],epsilon),'a+') as flat_group:
    F = np.loadtxt(flat_group.name)

ax = plt.gca()
plot1 = ax.plot(alpha[1:],F[1:,0],lw=lw,label=r'$\mathcal{F},\epsilon=0$')
#plot2 = ax.plot(alpha[1:],F[1:,1],lw=lw,label=r'$\mathcal{F},\epsilon=1$')
plot3 = ax.plot(alpha[1:],F[1:,2],lw=lw,label=r'$\mathcal{F},\epsilon=2$')
#plot4 = ax.plot(alpha[1:],F[1:,3],lw=lw,label=r'$\mathcal{F},\epsilon=3$')
plot5 = ax.plot(alpha[1:],F[1:,4],lw=lw,label=r'$\mathcal{F},\epsilon=4$')
#plot6 = ax.plot(alpha[1:],F[1:,5],lw=lw,label=r'$\mathcal{F},\epsilon=5$')
plot7 = ax.plot(alpha[1:],F[1:,6],lw=lw,label=r'$\mathcal{F},\epsilon=6$')
#plot8 = ax.plot(alpha[1:],F[1:,7],lw=lw,label=r'$\mathcal{F},\epsilon=7$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-0.01, 350)

#plt.grid()
plt.show()
plt.savefig('flat_2ndtry_mu=4_omega=3_2dof.pdf', format='pdf', \

            bbox_inches='tight')

#%% flatness aginist omega
num_ome=100 # number of alphas we want to calculate
alpha=1
epsilon = np.array([0,1,2,3,4,5,6,7])
#epsilon = np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75])
F = np.zeros((num_ome,len(epsilon)))
omega = np.linspace(0,10,num_ome)

#%% Perform the calculation
for i in range(num_ome):
    for j in range(len(epsilon)):
        parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon[j], omega[i]])
        
        F[i,j] = flatness_2dof_sn(-1,10,-5,5,500,parameters)
        print(F[i,j])
        
with open("flatness_2ndtry_minalpha=%smaxomega=%s_epsilon=%s.txt" %(omega[0],omega[-1],epsilon),'a+') as flat_group:
    np.savetxt(flat_group.name,F,fmt='%1.16e')
#%% 2dof flatness plottings
sns.set(font_scale = 2)
with open("flatness_2ndtry_minalpha=%smaxomega=%s_epsilon=%s.txt" %(omega[0],omega[-1],epsilon),'a+') as flat_group:
    F = np.loadtxt(flat_group.name)

ax = plt.gca()
plot1 = ax.plot(omega[1:],F[1:,0],lw=lw,label=r'$\mathcal{F},\epsilon=0$')
#plot2 = ax.plot(omega[1:],F[1:,1],lw=lw,label=r'$\mathcal{F},\epsilon=1$')
plot3 = ax.plot(omega[1:],F[1:,2],lw=lw,label=r'$\mathcal{F},\epsilon=2$')
#plot4 = ax.plot(omega[1:],F[1:,3],lw=lw,label=r'$\mathcal{F},\epsilon=3$')
plot5 = ax.plot(omega[1:],F[1:,4],lw=lw,label=r'$\mathcal{F},\epsilon=4$')
#plot6 = ax.plot(omega[1:],F[1:,5],lw=lw,label=r'$\mathcal{F},\epsilon=5$')
plot7 = ax.plot(omega[1:],F[1:,6],lw=lw,label=r'$\mathcal{F},\epsilon=6$')
#plot8 = ax.plot(omega[1:],F[1:,7],lw=lw,label=r'$\mathcal{F},\epsilon=7$')
legend = ax.legend(loc='best')
ax.set_xlabel(r'$\omega$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{F}$', fontsize=axis_fs)
ax.set_xlim(0, 10)
ax.set_ylim(-0.01, 350)

#plt.grid()
plt.show()
plt.savefig('flat_2ndtry_mu=4_alpha=1_2dof.pdf', format='pdf', \

            bbox_inches='tight')
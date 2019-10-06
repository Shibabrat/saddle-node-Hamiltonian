# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:11:01 2019

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
from scipy import optimize
from functools import partial
import saddlenode_tpcd ### import module xxx where xxx is the name of the python file xxx.py 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib import cm
from pylab import rcParams
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
%matplotlib

#%%
resX=100
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=0.5
#mu=alpha**(3/4)
mu=0.1
omega=1.0
epsilon=0.1

parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega])
eqNum = 1;  
eqPt = saddlenode_tpcd.get_eq_pts_saddlenode(eqNum, parameters)

#%%
"""Import initial conditions of the NHIM
"""
e=0.2
deltaE = 0.2
parameters = np.array([MASS_A, MASS_B, EPSILON_S, 0.5, mu, epsilon, omega])
po_fam_file = open("1111x0_tpcd_deltaE%s_saddlenode_parameters%s.txt" %(deltaE,parameters),'a+');
print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
x0podata = np.loadtxt(po_fam_file.name)
po_fam_file.close()
x0po_1 = x0podata[-1,0:4]

e=0.2
deltaE = 0.2
parameters = np.array([MASS_A, MASS_B, EPSILON_S, 1, mu, epsilon, omega])
po_fam_file = open("1111x0_tpcd_deltaE%s_saddlenode_parameters%s.txt" %(deltaE,parameters),'a+');
print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
x0podata = np.loadtxt(po_fam_file.name)
po_fam_file.close()
x0po_2 = x0podata[-1,0:4]

e=0.2
deltaE = 0.2
parameters = np.array([MASS_A, MASS_B, EPSILON_S, 2, mu, epsilon, omega])
po_fam_file = open("1111x0_tpcd_deltaE%s_saddlenode_parameters%s.txt" %(deltaE,parameters),'a+');
print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
x0podata = np.loadtxt(po_fam_file.name)
po_fam_file.close()
x0po_3 = x0podata[-1,0:4]
#%% 
"""Integrate the initial conditions to get data for the full period 
"""
TSPAN = [0,50]
plt.close('all')
axis_fs = 15
RelTol = 3.e-10
AbsTol = 1.e-10

f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
soln = solve_ivp(f, TSPAN, x0po_1,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
te = soln.t_events[0]
tt = [0,te[2]]
t1,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_1,parameters)
x_po1 = x[:,0]
y_po1 = x[:,1]
px_po1 = x[:,2]
py_po1 = x[:,3]

f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
soln = solve_ivp(f, TSPAN, x0po_2,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
te = soln.t_events[0]
tt = [0,te[2]]
t2,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_2,parameters)
x_po2 = x[:,0]
y_po2 = x[:,1]
px_po2 = x[:,2]
py_po2 = x[:,3]
f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
soln = solve_ivp(f, TSPAN, x0po_3,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
te = soln.t_events[0]
tt = [0,te[2]]
t3,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_3,parameters)
x_po3 = x[:,0]
y_po3 = x[:,1]
px_po3 = x[:,2]
py_po3 = x[:,3]

#%% Calculate reaction flux 
"""Calculate reaction flux
"""
q1=0
q2=0
q3=0
for i in range(len(px_po1)-1):
    q1 = q1+ 0.5*(t1[i+1]-t1[i])*(py_po1[i+1]**2+py_po1[i]**2) + 0.5*(t1[i+1]-t1[i])*(px_po1[i+1]**2+px_po1[i]**2)
for i in range(len(px_po2)-1):
    q2 = q2+ 0.5*(t2[i+1]-t2[i])*(py_po2[i+1]**2+py_po2[i]**2) + 0.5*(t2[i+1]-t2[i])*(px_po2[i+1]**2+px_po2[i]**2)
for i in range(len(px_po3)-1):
    q3 = q3+ 0.5*(t3[i+1]-t3[i])*(py_po3[i+1]**2+py_po3[i]**2) + 0.5*(t3[i+1]-t3[i])*(px_po3[i+1]**2+px_po3[i]**2)

Q = np.array([q1,q2,q3])    
print("flux is (%s,%s,%s)" % (q1,q2,q3))

#%% Plot reaction flux as a function of Depth, Flatness

def depth(alpha,mu,omega,epsilon):
    depth = -(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)
    return depth


epsilon = np.array([0,0.05,0.1,0.2])
alpha = np.array([0.5,1,2])   
D = np.zeros((3,4))
for i in range(len(alpha)):
    for j in range(len(epsilon)):
        D[i,j] = depth(alpha[i],mu,omega,epsilon[j])
ax = plt.gca()
plot1 = ax.plot(alpha[0:],Q,label=r'$\mathcal{F},\epsilon=0$')
plot2 = ax.plot(D[:,1],Q,label=r'$\mathcal{F},\epsilon=0$')


legend = ax.legend(loc='best')
ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{Q}$', fontsize=axis_fs)
ax.set_xlim(-1, 2)
ax.set_ylim(1, 1.5)

plt.grid()
plt.show()
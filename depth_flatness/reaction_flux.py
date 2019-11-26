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
import diffcorr_UPOsHam2dof ### import module xxx where xxx is the name of the python file xxx.py 
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
epsilon=0.05

parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega])
eqNum = 1
#model = 'saddlenode'
#eqPt = diffcorr_UPOsHam2dof.get_eq_pts(eqNum,model, parameters)
eqPt = diffcorr_UPOsHam2dof.get_eq_pts(eqNum, init_guess_eqpt_saddlenode, \
                                       grad_pot_saddlenode, parameters)

#energy of the saddle eq pt
eSaddle = diffcorr_UPOsHam2dof.get_total_energy([eqPt[0],eqPt[1],0,0], pot_energy_saddlenode, \
                                                parameters) 

#%%
def varEqns_saddlenode(t,PHI,par):
    
    """
    PHIdot = varEqns_saddlenode(t,PHI) 
    
    This here is a preliminary state transition, PHI(t,t0),
    matrix equation attempt for a ball rolling on the surface, based on...
    
    d PHI(t, t0)
    ------------ =  Df(t) * PHI(t, t0)
        dt
    
    """
    
    phi = PHI[0:16]
    phimatrix  = np.reshape(PHI[0:16],(4,4))
    x,y,px,py = PHI[16:20]
    
    
    # The first order derivative of the potential energy.
    dVdx = par[3]*x**2-2*np.sqrt(par[4])*x+par[5]*(x-y)
    dVdy = par[6]**2*y-par[5]*(x-y)

    # The second order derivative of the potential energy. 
    d2Vdx2 = 2*par[3]*x-2*np.sqrt(par[4])+par[5]
            
    d2Vdy2 = par[6]**2+par[5]
    
    d2Vdydx = -par[5]
    
    d2Vdxdy = d2Vdydx    

    Df    = np.array([[  0,     0,    par[0],    0],
              [0,     0,    0,    par[1]],
              [-d2Vdx2,  -d2Vdydx,   0,    0],
              [-d2Vdxdy, -d2Vdy2,    0,    0]])

    
    phidot = np.matmul(Df, phimatrix) # variational equation

    PHIdot        = np.zeros(20)
    PHIdot[0:16]  = np.reshape(phidot,(1,16)) 
    PHIdot[16]    = px/par[0]
    PHIdot[17]    = py/par[1]
    PHIdot[18]    = -dVdx 
    PHIdot[19]    = -dVdy
    
    return list(PHIdot)

def ham2dof_saddlenode(t, x, par):
    """
    Hamilton's equations of motion
    """
    
    xDot = np.zeros(4)
    
    dVdx = par[3]*x[0]**2-2*np.sqrt(par[4])*x[0]+par[5]*(x[0]-x[1])
    dVdy = par[6]**2*x[1]-par[5]*(x[0]-x[1])
        
    xDot[0] = x[2]/par[0]
    xDot[1] = x[3]/par[1]
    xDot[2] = -dVdx 
    xDot[3] = -dVdy
    
    return list(xDot)    

def half_period_saddlenode(t, x, par):
    """
    Return the turning point where we want to stop the integration                           
    
    pxDot = x[0]
    pyDot = x[1]
    xDot = x[2]
    yDot = x[3]
    """
    
    terminal = True
    # The zero can be approached from either direction
    direction = 0 #0: all directions of crossing
    
    return x[3]
#%% Calculate reaction flux
"""Import initial conditions of the NHIM
"""
e=0.2
deltaE = 0.2
num_alpha = 10
num_ep = 4
alpha = np.zeros((num_alpha,num_ep))
alpha[:,0] = np.linspace(0.1,5,num_alpha)
alpha[:,1] = np.linspace(0.03,1,num_alpha)
alpha[:,2] = np.linspace(0.01,0.5,num_alpha)
alpha[:,3] = np.linspace(0.001,0.1,num_alpha)
epsilon = np.array([1e-20,0.5,1.0,1.5])

"""Save reaction flux in a matrix Q using a for loop.
For each value of alpha, epsilon, import data from the data folder, 
calculate the reaction flux and save it in the matrix Q.
"""
Q = np.zeros((num_alpha,num_ep))
Q[0,0] = 0.4*np.pi
Q[1,0] = 0.4*np.pi
Q[2,0] = 0.4*np.pi
for i in range(num_alpha):
    for j in range(num_ep):
        parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i,j], mu, epsilon[j], omega])
        po_fam_file = open("x0_diffcorr_alpha%s_ep%s_deltaE%s_saddlenode.txt" %(alpha[i,j],epsilon[j],deltaE),'a+')
        print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
        x0podata = np.loadtxt(po_fam_file.name)
        po_fam_file.close()
        x0po_1 = x0podata[0:4]
        TSPAN = [0,50]
        plt.close('all')
        axis_fs = 15
        RelTol = 3.e-10
        AbsTol = 1.e-10
        f= lambda t,x: ham2dof_saddlenode(t,x,parameters)
        soln = solve_ivp(f, TSPAN, x0po_1,method='RK45',dense_output=True, \
                     events = lambda t,x : half_period_saddlenode(t,x,parameters), \
                     rtol=RelTol, atol=AbsTol)
#        f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
#        soln = solve_ivp(f, TSPAN, x0po_1,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
        te = soln.t_events[0]
        tt = [0,te[2]]
#        print(tt)
#        t1,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_1,parameters)
        t1,x,phi_t1,PHI = diffcorr_UPOsHam2dof.stateTransitMat(tt,x0po_1,parameters,varEqns_saddlenode)
        x_po1 = x[:,0]
        y_po1 = x[:,1]
        px_po1 = x[:,2]
        py_po1 = x[:,3]
        
        q1=0
        for k in range(len(px_po1)-1):
            q1 = q1+ 0.5*(t1[k+1]-t1[k])*(py_po1[k+1]**2+py_po1[k]**2) + 0.5*(t1[k+1]-t1[k])*(px_po1[k+1]**2+px_po1[k]**2)
        Q[i,j] = q1


#%%
#"""Import initial conditions of the NHIM
#"""
#epsilon=0.2
#e=0.2
#deltaE = 0.2
#parameters = np.array([MASS_A, MASS_B, EPSILON_S, 0.5, mu, epsilon, omega])
#po_fam_file = open("1111x0_tpcd_deltaE%s_saddlenode_parameters%s.txt" %(deltaE,parameters),'a+');
#print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
#x0podata = np.loadtxt(po_fam_file.name)
#po_fam_file.close()
#x0po_1 = x0podata[-1,0:4]
#
#e=0.2
#deltaE = 0.2
#parameters = np.array([MASS_A, MASS_B, EPSILON_S, 1, mu, epsilon, omega])
#po_fam_file = open("1111x0_tpcd_deltaE%s_saddlenode_parameters%s.txt" %(deltaE,parameters),'a+');
#print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
#x0podata = np.loadtxt(po_fam_file.name)
#po_fam_file.close()
#x0po_2 = x0podata[-1,0:4]
#
#e=0.2
#deltaE = 0.2
#parameters = np.array([MASS_A, MASS_B, EPSILON_S, 2, mu, epsilon, omega])
#po_fam_file = open("1111x0_tpcd_deltaE%s_saddlenode_parameters%s.txt" %(deltaE,parameters),'a+');
#print('Loading the periodic orbit family from data file',po_fam_file.name,'\n'); 
#x0podata = np.loadtxt(po_fam_file.name)
#po_fam_file.close()
#x0po_3 = x0podata[-1,0:4]
##%% 
#"""Integrate the initial conditions to get data for the full period 
#"""
#TSPAN = [0,50]
#plt.close('all')
#axis_fs = 15
#RelTol = 3.e-10
#AbsTol = 1.e-10
#
#f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
#soln = solve_ivp(f, TSPAN, x0po_1,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
#te = soln.t_events[0]
#tt = [0,te[2]]
#t1,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_1,parameters)
#x_po1 = x[:,0]
#y_po1 = x[:,1]
#px_po1 = x[:,2]
#py_po1 = x[:,3]
#
#f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
#soln = solve_ivp(f, TSPAN, x0po_2,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
#te = soln.t_events[0]
#tt = [0,te[2]]
#t2,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_2,parameters)
#x_po2 = x[:,0]
#y_po2 = x[:,1]
#px_po2 = x[:,2]
#py_po2 = x[:,3]
#f = partial(saddlenode_tpcd.saddlenode2dof, par=parameters) 
#soln = solve_ivp(f, TSPAN, x0po_3,method='RK45',dense_output=True, events = saddlenode_tpcd.half_period,rtol=RelTol, atol=AbsTol)
#te = soln.t_events[0]
#tt = [0,te[2]]
#t3,x,phi_t1,PHI = saddlenode_tpcd.stateTransitMat_saddlenode(tt,x0po_3,parameters)
#x_po3 = x[:,0]
#y_po3 = x[:,1]
#px_po3 = x[:,2]
#py_po3 = x[:,3]
#
##%% Calculate reaction flux 
#"""Calculate reaction flux
#"""
#q1=0
#q2=0
#q3=0
#for i in range(len(px_po1)-1):
#    q1 = q1+ 0.5*(t1[i+1]-t1[i])*(py_po1[i+1]**2+py_po1[i]**2) + 0.5*(t1[i+1]-t1[i])*(px_po1[i+1]**2+px_po1[i]**2)
#for j in range(len(px_po2)-1):
#    q2 = q2+ 0.5*(t2[j+1]-t2[j])*(py_po2[j+1]**2+py_po2[j]**2) + 0.5*(t2[j+1]-t2[j])*(px_po2[j+1]**2+px_po2[j]**2)
#for k in range(len(px_po3)-1):
#    q3 = q3+ 0.5*(t3[k+1]-t3[k])*(py_po3[k+1]**2+py_po3[k]**2) + 0.5*(t3[k+1]-t3[k])*(px_po3[k+1]**2+px_po3[k]**2)
#
#Q = np.array([q1,q2,q3])    
#print("flux is (%s,%s,%s)" % (q1,q2,q3))


#%% Plot reaction flux as a function of Depth, Flatness
def depth(alpha,mu,omega,epsilon):
    """Definition of depth
    """
    depth = -(2*math.sqrt(mu)- (omega**2*epsilon)/(omega**2+epsilon))**3/(6*alpha**2)
    return depth


e=0.2
deltaE = 0.2
num_alpha = 10
num_ep = 4
alpha = np.zeros((num_alpha,num_ep))
alpha[:,0] = np.linspace(0.1,5,num_alpha)
alpha[:,1] = np.linspace(0.03,1,num_alpha)
alpha[:,2] = np.linspace(0.01,0.5,num_alpha)
alpha[:,3] = np.linspace(0.001,0.1,num_alpha)
epsilon = np.array([1e-20,0.5,1.0,1.5]) 

#file = open("reaction_flux_e=0.2_alpha=%s_ep=%s.txt" %(alpha,epsilon),'a+');
#np.savetxt(file.name,reaction_flux,fmt='%1.16e')
#po_fam_file.close()
 
D = np.zeros((num_alpha,num_ep))
for i in range(num_alpha):
    for j in range(num_ep):
        D[i,j] = depth(alpha[i,j],mu,omega,epsilon[j])
ax = plt.gca()
plot1 = ax.scatter(D[:,0],Q[:,0],label=r'$\mathcal{Q},\epsilon=0$')
plot2 = ax.scatter(D[:,1],Q[:,1],label=r'$\mathcal{Q},\epsilon=0.5$')
plot3 = ax.scatter(D[:,2],Q[:,2],label=r'$\mathcal{Q},\epsilon=1$')
plot4 = ax.scatter(D[:,3],Q[:,3],label=r'$\mathcal{Q},\epsilon=1.5$')

#plot1 = ax.plot(alpha[:,0],Q[:,0],label=r'$\mathcal{Q},\epsilon=0$')
#plot2 = ax.plot(alpha[:,1],Q[:,1],label=r'$\mathcal{Q},\epsilon=0.5$')
#plot3 = ax.plot(alpha[:,2],Q[:,2],label=r'$\mathcal{Q},\epsilon=1$')
#plot4 = ax.plot(alpha[:,3],Q[:,3],label=r'$\mathcal{Q},\epsilon=1.5$')

legend = ax.legend(loc='best')
ax.set_xlabel(r'$\mathcal{D}$', fontsize=axis_fs)
ax.set_ylabel(r'$\mathcal{Q}$', fontsize=axis_fs)
ax.set_xlim(-0.2, 0)
ax.set_ylim(0.1, 1.3)

#legend = ax.legend(loc='best')
#ax.set_xlabel(r'$\alpha$', fontsize=axis_fs)
#ax.set_ylabel(r'$\mathcal{Q}$', fontsize=axis_fs)
#ax.set_xlim(0.0, 2)
#ax.set_ylim(0.1, 1.3)

plt.grid()
plt.show()
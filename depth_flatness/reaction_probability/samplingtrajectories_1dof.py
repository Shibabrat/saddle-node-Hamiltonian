# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 19:20:52 2019

@author: wl16298
"""

#%%
from scipy.integrate import odeint,quad,trapz,solve_ivp, ode
from IPython.display import Image # for Notebook
from IPython.core.display import HTML
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
import scipy.linalg as linalg
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import math
from matplotlib import cm
from scipy import optimize
from scipy.optimize import fsolve
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
#%% illustration of the choice of the initional conditions on the configuration space and the evolution of the trajectories.
# 1dof
resX=200
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=2.0
mu = 4

parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu]);

def H_SN1dof(x,px,par):
    H = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2 + 0.5*px**2
    return H

def V_SN_1dof(x,par):
    V = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2
    return V

def get_H_surf_proj(xVec, yVec,par):            

    resX = np.size(xVec)
    resY = np.size(xVec)
    surfProj = np.zeros([resX, resY])
    for i in range(len(xVec)):
        for j in range(len(yVec)):
            surfProj[i,j] = H_SN1dof(xVec[j], yVec[i],par)

    return surfProj

def saddlenode1dof(t,x, par):
    # Hamilton's equations, used for ploting trajectories later.
    xDot = np.zeros(2)

    dVdx = par[3]*(x[0])**2-2*math.sqrt(par[4])*(x[0])
    
    xDot[0] = x[1]*par[0]
    xDot[1] = -dVdx 
    return list(xDot)

def reaction_event(t,x,par):
# reaction is defined when the trajectory reaches the DS(x=0)                            

    terminal = True
    # The zero can be approached from either direction
    direction = 0; #0: all directions of crossing
    return x[0]

reaction_event.terminal = True 

# The zero can be approached from either direction

reaction_event.direction=0#0: all directions of crossing

def reaction_probability(x,px,e,par):
    """
    This function records the time when the trajectory reaches the DS(x=0) 
    from some initial conditions in the reactants region(x>0).
    num_sim is the number of simulations we want to perform at a particular point on the configuration space 
    Parameters:
        x: initial condition of x coordinate
        px: initial condition of px coordinate
        e: total energy of the system
        par: parameter values
        num_sim: number of simulation we want to perform for given values of x, y.
    Returns:
        event_t: time when the trajecotry reaches the DS(x=0)
    """
    TSPAN = [0,40]
    RelTol = 3.e-10
    AbsTol = 1.e-10

    f1 = lambda t,x: saddlenode1dof(t,x,par) 
    soln1 = solve_ivp(f1, TSPAN, [x,px], \
                                  method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,par),rtol=RelTol, atol=AbsTol)
    te1 = soln1.t_events[0]
    event_t = te1[0]
    print("exits products region at time%s" %(te1))
    return event_t

def sampling_phase(num_pts,par,e):
    """This function returns the sampling points for calculating reaction probability in configuration space x-y.
       Returns all the points in the domain x \in [x_min, x_max], y \in [y_min, y_max] s.t. V(x,y) <=e
       
       Variables:
           num_pts is the number of points in each x, y direction.(In total there are num_pts**2 points in the domain)
           par is the parameter values
           e is the total energy of the system
           (all the points in the configuration space have energy <= e)
       
       Returns:
           accept_pts are the positions of all possible sampling points 
    """
    TSPAN = [0,50]
    RelTol = 3.e-10
    AbsTol = 1.e-10 
    x_max_guess = 2*np.sqrt(par[4])/par[3]*4/3
    x_max = fsolve(lambda x: H_SN1dof(x,0,par)-e,x_max_guess)
    x_max=x_max[0]
    f1 = lambda t,x: saddlenode1dof(t,x,par) 
    soln1 = solve_ivp(f1, TSPAN, [x_max,0],\
                              method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,par),rtol=RelTol, atol=AbsTol)
    coordinates = soln1.y
    
#
#    x = np.linspace(0,x_max,num_pts)
#    px = -np.sqrt(2*e-2*V_SN_1dof(x,par))
#    accept_pts=np.reshape((x,px),(2,num_pts))
#    accept_pts=np.transpose(accept_pts)
    boundary_value = np.where(coordinates[0,:]>0)[0][-1]
    accept_pts=coordinates[:,:boundary_value]
    return accept_pts

#%% depth and flatness
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

num_alp=3 # number of alphas we want to calculate
alpha = np.array([1,2,5])
D = depth_1dof_sn([MASS_A, MASS_B, EPSILON_S,alpha,mu])
F = np.zeros(num_alp)
for i in range(num_alp):
    parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i], mu])
    F[i] = flatness_1dof_sn(-1,10,500,parameters)
    print(F[i])

#%% Illustration of the initial conditions in the phase space
xVec = np.linspace(-2,6,resX)
pxVec = np.linspace(-5,5,resX)
xMat, yMat = np.meshgrid(xVec, pxVec)
axis_fs=30
alpha=2.0
mu = 4
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu]);
sns.set(font_scale = 2)
e=0.5
ax = plt.gca()
# projection of the PES on the x-y configuration space
cset1 = ax.contour(xMat, yMat, get_H_surf_proj(xVec, pxVec,parameters), [-0.5],\
                   linewidths = 1.0, cmap=cm.viridis, alpha = 0.8)
cset2 = ax.contour(xMat, yMat, get_H_surf_proj(xVec, pxVec,parameters), [0.0],\
                   linewidths = 1.0, cmap=cm.viridis, alpha = 0.8)
cset3 = ax.contour(xMat, yMat, get_H_surf_proj(xVec, pxVec,parameters), [0.5],\
                   linewidths = 1.0, colors = 'b', alpha = 0.8)
cset5 = ax.contour(xMat, yMat, get_H_surf_proj(xVec, pxVec,parameters), [1.0],\
                   linewidths = 1.0, cmap=cm.viridis, alpha = 0.8)
        
    
sampling_points = sampling_phase(100,parameters,0.5)
for i in range(len(sampling_points[0,:])):
    # plot the sampling points in the configuration space
    ax.scatter(sampling_points[0,i],sampling_points[1,i], s = 5, c = 'k', marker = '.')
#TSPAN = [0,10]
#RelTol = 3.e-10
#AbsTol = 1.e-10 
#f1 = lambda t,x: saddlenode1dof(t,x,parameters) 
## randomly pick a x and px values on the configuration space.
#pick_x = 0
#e=0.5
#pick_px = np.sqrt(2*e-2*pick_x)
## integrate the Hamilton's equation and plot the trajectory, with negative px
#soln1 = solve_ivp(f1, TSPAN, [pick_x,pick_px],\
#                              method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,parameters),rtol=RelTol, atol=AbsTol)
#te1 = soln1.t_events[0]
#coordinates = soln1.y
#ax.plot(coordinates[0,:],coordinates[1,:],c = 'g')
x_max = fsolve(lambda x: H_SN1dof(x,0,parameters)-e,2*np.sqrt(mu)/alpha*4/3)
def equ_pos(par):
    """Returns the position of the centre equilibrium point
    """
    xe = 2*np.sqrt(par[4])/par[3]
    pxe = 0
    return xe, pxe

ax.set_xlim(-2, 6)
ax.set_ylim(-5, 5)
ax.set_xlabel('$x$', fontsize=axis_fs)
ax.set_ylabel('$p_x$', fontsize=axis_fs)
#ax.scatter(x_max,0,s = 100, c = 'r', marker = 'o')
ax.scatter(0,0, s = 100, c = 'r', marker = 'x')
ax.scatter(equ_pos(parameters)[0],equ_pos(parameters)[1], s = 100, c = 'r', marker = 'o')
#ax.set_title(r'$\alpha=%s, \mu=%s$' %(alpha, mu))
plt.show()
plt.savefig('samp_pos_mu=4_omega=3_1dof.pdf', format='pdf', \

            bbox_inches='tight')
#%% Reaction probability calculation
e=0.5
alpha1=1
alpha2=2
alpha3=5
parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu])

sampling_pts = sampling_phase(32,parameters,0.5)
sampling_pts = np.transpose(sampling_pts)
num_pts = len(sampling_pts)
with open("e=%s_par=%s_1dof.txt" %(e,parameters[3:]),'a+') as react_time_file:
    event_times = np.zeros(num_pts)
    for i in range(num_pts):
        print(i)
        event_t = reaction_probability(sampling_pts[i,0],sampling_pts[i,1],e,parameters)
        event_times[i] = event_t
    np.savetxt(react_time_file.name,event_times,fmt='%1.16e')

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu])
sampling_pts = sampling_phase(50,parameters,0.5)
sampling_pts = np.transpose(sampling_pts)
num_pts = len(sampling_pts)
with open("e=%s_par=%s_1dof.txt" %(e,parameters[3:]),'a+') as react_time_file:
    event_times = np.zeros(num_pts)
    for i in range(num_pts):
        print(i)
        event_t = reaction_probability(sampling_pts[i,0],sampling_pts[i,1],e,parameters)
        event_times[i] = event_t
    np.savetxt(react_time_file.name,event_times,fmt='%1.16e')

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu])
sampling_pts = sampling_phase(75,parameters,0.5)
sampling_pts = np.transpose(sampling_pts)
num_pts = len(sampling_pts)
with open("e=%s_par=%s_1dof.txt" %(e,parameters[3:]),'a+') as react_time_file:
    event_times = np.zeros(num_pts)
    for i in range(num_pts):
        print(i)
        event_t = reaction_probability(sampling_pts[i,0],sampling_pts[i,1],e,parameters)
        event_times[i] = event_t
    np.savetxt(react_time_file.name,event_times,fmt='%1.16e')
    
#%% load data
e=0.5
alpha1=1
alpha2=2
alpha3=5

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu])
with open("e=%s_par=%s_1dof.txt" %(e,parameters[3:]),'a+') as react_time_file:
    print('Loading the reaction time from data file',react_time_file.name,'\n') 
    event_t1 = np.loadtxt(react_time_file.name)


parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu])
with open("e=%s_par=%s_1dof.txt" %(e,parameters[3:]),'a+') as react_time_file:
    print('Loading the reaction time from data file',react_time_file.name,'\n') 
    event_t2 = np.loadtxt(react_time_file.name)


parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu])
with open("e=%s_par=%s_1dof.txt" %(e,parameters[3:]),'a+') as react_time_file:
    print('Loading the reaction time from data file',react_time_file.name,'\n') 
    event_t3 = np.loadtxt(react_time_file.name)

#%% Plot the reaction probability as a function of time
axis_fs=30
sns.set(font_scale = 2)
event_t1 = np.sort(event_t1,axis=None)
event_t2 = np.sort(event_t2,axis=None)
event_t3 = np.sort(event_t3,axis=None)
event_time = [event_t1, event_t2, event_t3]
alpha = [1,2,5]
ax = plt.gca()

for i in range(3):
    event_t = event_time[i]

    react_prob = np.zeros((len(event_t),1))
    for j in range(len(event_t)):
        react_prob[j] = (j+1)/len(event_t)
    ax.plot(event_t,react_prob, label=r'$\alpha=%s,\mathcal{D}=%.2f,\mathcal{F}=%.2f$'\
            %(alpha[i],D[i],F[i]))

ax.set_xlabel('$t$', fontsize=axis_fs)
ax.set_ylabel('reaction probability', fontsize=axis_fs)
ax.set_xlim(0, 5)
ax.set_ylim(0, 1)
legend = ax.legend(loc='best')
plt.show()
plt.savefig('reactprob2_alp=1_2_5_mu4_e=0dot5_1of.pdf', format='pdf', \

            bbox_inches='tight')
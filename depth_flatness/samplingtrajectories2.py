# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 12:29:15 2019

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
fal = 30 # fontsize axis labels
ftl = 30  #fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal
#m = cm.ScalarMappable(cmap=cm.jet)

#%%
resX=100
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=1.0
#mu=alpha**(3/4)
mu = 0.1
omega=1.0
epsilon=0.0

parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega]);


#%%
def V_SN2dof(x,y,par):
    pot = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2 + 0.5*par[6]*y**2+0.5*par[5]*(x-y)**2
    return pot
    

def saddlenode2dof(t,x, par):
    # Hamilton's equations, used for ploting trajectories later.
    xDot = np.zeros(4)

    dVdx = par[3]*(x[0])**2-2*math.sqrt(par[4])*(x[0])+par[5]*(x[0]-x[1])
    dVdy = par[6]**2*x[1]-par[5]*(x[0]-x[1])
    
    xDot[0] = x[2]*par[0]
    xDot[1] = x[3]*par[0]
    xDot[2] = -dVdx 
    xDot[3] = -dVdy
    return list(xDot)

def reaction_event(t,x,par):
# reaction is defined when the trajectory reaches the DS(x=0)                            
    #pxDot = x[0]
    #pyDot = x[1]
    #xDot = x[2]
    #yDot = x[3]
    #xe = 2*math.sqrt(par[4])/par[3] - (par[6]**2*par[5])/(par[3]*(par[6]**2+par[5]))
    
    terminal = True
    # The zero can be approached from either direction
    direction = 0; #0: all directions of crossing
    return x[0]

def get_pot_surf_proj(xVec, yVec,par):            

    resX = np.size(xVec)
    resY = np.size(xVec)
    surfProj = np.zeros([resX, resY])
    for i in range(len(xVec)):
        for j in range(len(yVec)):
            surfProj[i,j] = V_SN2dof(xVec[j], yVec[i],par)

    return surfProj

def deter_py(x,y,px,e,par):
    """ this function determines the postive value of py for a given x,y,px and total energy e.
    
    Parameters:
        x: initial condition of x coordinate
        y: initial condition of y coordinate
        px: initial condition of px coordinate
        e: total energy of the system
        par: parameter values
    """
    py = np.sqrt(e-V_SN2dof(x,y,par)-0.5*px**2)
    return py

#%% defines function for calculating reaction probability, sampling points in the configuration space
def reaction_probability(x,y,e,par,num_sim):
    """
    This function records the time when the trajectory reaches the DS(x=0) 
    from some initial conditions in the reactants region(x>0).
    num_sim is the number of simulations we want to perform at a particular point on the configuration space 
    Parameters:
        x: initial condition of x coordinate
        y: initial condition of y coordinate
        e: total energy of the system
        par: parameter values
        num_sim: number of simulation we want to perform for given values of x, y.
    """
    TSPAN = [0,40]
    RelTol = 3.e-10
    AbsTol = 1.e-10
    px_max = deter_py(x,y,0,e,par) # maximum value of px we can get for a fixed x, y, total energy e and parameter values.
#    px = np.linspace(0,px_max,num_sim) # sampling positive px
    px = np.linspace(-px_max,0,num_sim) # sampling negative px
    event_t = np.zeros((num_sim,1))
    for i in range(num_sim):
        print(i)
        f1 = partial(saddlenode2dof, par=parameters) 
        reaction = partial(reaction_event, par=parameters)
        soln1 = solve_ivp(f1, TSPAN, [x,y,px[i],deter_py(x,y,px[i],e,par)], \
                                      method='RK45',dense_output=True, events = reaction,rtol=RelTol, atol=AbsTol)
        te1 = soln1.t_events[0]
        react_number=0
        if len(te1) > 0: # there is a possibility that the trajectory won't react,
                         # i.e. it will not reach the DS(x=0), which means te1 is [](empty).
            event_t[i] = te1[0]
            react_number = react_number+1
            
        print(te1)
        

    print("exits products region at time%s" %(te1))
    return event_t

def sampling_configuration(x_min,x_max,y_min,y_max,num_pts,par,e):
    """This function returns the sampling points for calculating reaction probability in configuration space x-y.
       Returns all the points in the domain x \in [x_min, x_max], y \in [y_min, y_max] s.t. V(x,y) <=e
       
       Variables:
           x_min, x_max, y_min, y_max define a domain in the configuration space x-y
           num_pts is the number of points in each x, y direction.(In total there are num_pts**2 points in the domain)
           par is the parameter values
           e is the total energy of the system
           (all the points in the configuration space have energy <= e)
       
       Returns:
           accept_pts[:num_accept_pts,:] are the positions of all possible sampling points 
           accept_pts[:num_accept_pts,0] is the x coordinate of the sampling points
           accept_pts[:num_accept_pts,1] is the y coordinate of the sampling points
           num_accept_pts is the total number of accepted points(out of num_pts**2 points).
    """
    x = np.linspace(x_min,x_max,num_pts)
    y = np.linspace(y_min,y_max,num_pts)
    num_accept_pts = 0
    accept_pts = np.zeros((num_pts**2,2))
    for i in range(num_pts):
        for j in range(num_pts):
            if V_SN2dof(x[i],y[j],par) <= e:
                num_accept_pts = num_accept_pts+1
                accept_pts[num_accept_pts,:] = [x[i],y[j]]
    return accept_pts[:num_accept_pts,:]
#%%
num_sim=100 # The expected time to run the simulation for 3 different alphas, when num_sim=100 is approx 100 min.
e=0.5
alpha1=1
alpha2=2
alpha3=5
epsilon= 0.0
parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu, epsilon, omega])

DS = np.linspace(-np.sqrt(2*e/(omega**2+epsilon)),np.sqrt(2*e/(omega**2+epsilon)),12)
num_pts=10 #num_pts is the number of points we want to sample in one direction on the configuration space.
chosen_pts = np.zeros((num_pts,num_pts))
# Record the choice of initial conditions on the configuration space.
for i in range(10):
    V_max= fsolve(lambda x: V_SN2dof(x,DS[1:11][i],parameters)-e,1.5)
    V_x_range= np.linspace(0,V_max,12)
    chosen_pts[i] = V_x_range[1:11]
    
    
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
event_times = np.zeros((num_sim,1))
for i in range(10):
    print(i)
    for j in range(10):
        event_t = reaction_probability(chosen_pts[i,j],DS[1:11][i],e,parameters,num_sim,react_time_file)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
np.savetxt(react_time_file.name,event_times,fmt='%1.16e')
react_time_file.close()
        
parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu, epsilon, omega])

DS = np.linspace(-np.sqrt(2*e/(omega**2+epsilon)),np.sqrt(2*e/(omega**2+epsilon)),12)
num_pts=10
chosen_pts = np.zeros((num_pts,num_pts))
for i in range(10):
    V_max= fsolve(lambda x: V_SN2dof(x,DS[1:11][i],parameters)-e,1.5)
    V_x_range= np.linspace(0,V_max,12)
    chosen_pts[i] = V_x_range[1:11]
    
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
event_times = np.zeros((num_sim,1))
for i in range(10):
    print(i)
    for j in range(10):
        event_t = reaction_probability(chosen_pts[i,j],DS[1:11][i],e,parameters,num_sim,react_time_file)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
np.savetxt(react_time_file.name,event_times,fmt='%1.16e')
react_time_file.close()

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu, epsilon, omega])

DS = np.linspace(-np.sqrt(2*e/(omega**2+epsilon)),np.sqrt(2*e/(omega**2+epsilon)),12)
num_pts=10
chosen_pts = np.zeros((num_pts,num_pts))
for i in range(10):
    V_max= fsolve(lambda x: V_SN2dof(x,DS[1:11][i],parameters)-e,1.5)
    V_x_range= np.linspace(0,V_max,12)
    chosen_pts[i] = V_x_range[1:11]
    

react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
event_times = np.zeros((num_sim,1))
for i in range(10):
    print(i)
    for j in range(10):
        event_t = reaction_probability(chosen_pts[i,j],DS[1:11][i],e,parameters,num_sim,react_time_file)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
np.savetxt(react_time_file.name,event_times,fmt='%1.16e')
react_time_file.close()

#%%
num_sim=100
x=-0.1
y=0.5
e=0.5
alpha1=1
alpha2=2
alpha3=5
epsilon= 0.0
parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu, epsilon, omega])
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
print('Loading the reaction time from data file',react_time_file.name,'\n') 
event_t1 = np.loadtxt(react_time_file.name)
event_t1 = event_t1[num_sim:]
react_time_file.close()

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu, epsilon, omega])
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
print('Loading the reaction time from data file',react_time_file.name,'\n') 
event_t2 = np.loadtxt(react_time_file.name)
event_t2 = event_t2[num_sim:]
react_time_file.close()

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu, epsilon, omega])
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
print('Loading the reaction time from data file',react_time_file.name,'\n') 
event_t3 = np.loadtxt(react_time_file.name)
event_t3 = event_t3[num_sim:]
react_time_file.close()

#%% Plot the reaction probability as a function of time
axis_fs=30
event_t1 = np.sort(event_t1,axis=None)
event_t2 = np.sort(event_t2,axis=None)
event_t3 = np.sort(event_t3,axis=None)
event_time = [event_t1, event_t2, event_t3]
alpha = [1,2,5]
ax = plt.gca()

for i in range(3):
    event_t = event_time[i]
    non_react = len(np.where(event_t==0)[0])
    react_time=event_t[non_react:]
    react_prob = np.zeros((len(react_time),1))
    for j in range(len(react_time)):
        react_prob[j] = (j+1)/len(event_t)
    ax.plot(react_time,react_prob, label=r'$\alpha=%s,\mu=%s,\epsilon=%s,\omega=%s,e=%s$'\
            %(alpha[i],parameters[4],epsilon,parameters[6],e))

ax.set_xlabel('$t$', fontsize=axis_fs)
ax.set_ylabel('reaction probability', fontsize=axis_fs)
ax.set_xlim(0, 5)
ax.set_ylim(0, 1)
legend = ax.legend(loc='best')
plt.show()
#%% illustration of the choice of the initional conditions on the configuration space and the evolution of the trajectories.
resX=100
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=1.0
#mu=alpha**(3/4)
mu = 4
omega=3.0
epsilon=0

parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega]);

def V_SN2dof(x,y,par):
    pot = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2 + 0.5*par[6]*y**2+0.5*par[5]*(x-y)**2
    return pot

xVec = np.linspace(-2,6,resX)
yVec = np.linspace(-3,3,resX)
xMat, yMat = np.meshgrid(xVec, yVec)

ax = plt.gca()
# projection of the PES on the x-y configuration space
cset1 = ax.contour(xMat, yMat, get_pot_surf_proj(xVec, yVec,parameters), [-0.5],\
                   linewidths = 1.0, cmap=cm.viridis, alpha = 0.8)
cset2 = ax.contour(xMat, yMat, get_pot_surf_proj(xVec, yVec,parameters), [0.0],\
                   linewidths = 1.0, cmap=cm.viridis, alpha = 0.8)
cset3 = ax.contour(xMat, yMat, get_pot_surf_proj(xVec, yVec,parameters), [0.5],\
                   linewidths = 1.0, colors = 'b', alpha = 0.8)
cset5 = ax.contour(xMat, yMat, get_pot_surf_proj(xVec, yVec,parameters), [1.0],\
                   linewidths = 1.0, cmap=cm.viridis, alpha = 0.8)
        
def sampling_configuration(x_min,x_max,y_min,y_max,num_pts,par,e):
    """This function returns the sampling points for calculating reaction probability in configuration space x-y.
       Returns all the points in the domain x \in [x_min, x_max], y \in [y_min, y_max] s.t. V(x,y) <=e
       
       Variables:
           x_min, x_max, y_min, y_max define a domain in the configuration space x-y
           num_pts is the number of points in each x, y direction.(In total there are num_pts**2 points in the domain)
           par is the parameter values
           e is the total energy of the system
           (all the points in the configuration space have energy <= e)
       
       Returns:
           accept_pts[:num_accept_pts,:] are the positions of all possible sampling points 
           accept_pts[:num_accept_pts,0] is the x coordinate of the sampling points
           accept_pts[:num_accept_pts,1] is the y coordinate of the sampling points
           num_accept_pts is the total number of accepted points(out of num_pts**2 points).
    """
    x = np.linspace(x_min,x_max,num_pts)
    y = np.linspace(y_min,y_max,num_pts)
    num_accept_pts = 0
    accept_pts = np.zeros((num_pts**2,2))
    for i in range(num_pts):
        for j in range(num_pts):
            if V_SN2dof(x[i],y[j],par) <= e:
                num_accept_pts = num_accept_pts+1
                accept_pts[num_accept_pts,:] = [x[i],y[j]]
    return accept_pts[:num_accept_pts,:]

sampling_points = sampling_configuration(0,6,-3,3,100,parameters,0.5)
for i in range(len(sampling_points)):
    ax.scatter(sampling_points[i,0],sampling_points[i,1], s = 5, c = 'b', marker = '.')
TSPAN = [0,40]
RelTol = 3.e-10
AbsTol = 1.e-10 
f1 = partial(saddlenode2dof, par=parameters) 
reaction = partial(reaction_event, par=parameters)
# randomly pick a x and y values on the configuration space.
pick_x = sampling_points[int(len(sampling_points)/102),0]
pick_y = sampling_points[int(len(sampling_points)/102),1]
# integrate the Hamilton's equation and plot the trajectory, with negative px
soln1 = solve_ivp(f1, TSPAN, [pick_x,pick_y,-0.1,-deter_py(pick_x,pick_y,-0.1,0.5,parameters)],\
                              method='RK45',dense_output=True, events = reaction,rtol=RelTol, atol=AbsTol)
te1 = soln1.t_events[0]
coordinates = soln1.y
ax.plot(coordinates[0,:],coordinates[1,:])

def equ_pos(par):
    """Returns the position of the centre equilibrium point
    """
    xe = 2*math.sqrt(par[4])/par[3] - (par[6]**2*par[5])/(par[3]*(par[6]**2+par[5]))
    ye = xe*par[5]/(par[6]**2+par[5])
    return xe, ye

ax.set_xlim(-2, 6)
ax.set_ylim(-3, 3)
ax.scatter(0,0, s = 100, c = 'r', marker = 'x')
ax.scatter(equ_pos(parameters)[0],equ_pos(parameters)[1], s = 100, c = 'r', marker = 'o')
ax.set_title(r'$\alpha=%s, \mu=%s$' %(alpha, mu))
plt.show()


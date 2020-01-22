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
alpha=1.0
#mu=alpha**(3/4)
mu = 4
omega=3.0
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
# reaction is defined when the trajectory reaches the x=-5, to make sure the reaction occurs                          

    terminal = True
    # The zero can be approached from either direction
    direction = 0; #0: all directions of crossing
    return x[0]+5

reaction_event.terminal = True 

# The zero can be approached from either direction

reaction_event.direction=0#0: all directions of crossing


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

#%% defines function for calculating reaction probability, procedure for sampling points in the configuration space
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
    Returns:
        te1: time when the trajecotry reaches the DS(x=0),
             empty if the trajecotory does not reach the DS(x=0).
    """
    TSPAN = [0,40]
    RelTol = 3.e-10
    AbsTol = 1.e-10
    px_max = deter_py(x,y,0,e,par) # maximum value of px we can get for a fixed x, y, total energy e and parameter values.
#    px = np.linspace(0,px_max,num_sim) # sampling positive px
    px = np.linspace(-px_max,0,num_sim) # sampling negative px
    event_t = np.zeros((num_sim*2,1))
    # sampling postive py
    for i in range(num_sim):
        print(i)
        f1 = lambda t,x: saddlenode2dof(t,x,par) 
        soln1 = solve_ivp(f1, TSPAN, [x,y,px[i],deter_py(x,y,px[i],e,par)], \
                                      method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,par),rtol=RelTol, atol=AbsTol)
        te1 = soln1.t_events[0]
        react_number=0
        if len(te1) > 0: # there is a possibility that the trajectory won't react,
                         # i.e. it will not reach the DS(x=0), which means te1 is [](empty).
            event_t[i] = te1[0]
            react_number = react_number+1
        print(te1)

    # sampling negative py
    for i in range(num_sim,2*num_sim):
        print(i)
        f1 = lambda t,x: saddlenode2dof(t,x,par) 
        soln1 = solve_ivp(f1, TSPAN, [x,y,px[i-num_sim],-deter_py(x,y,px[i-num_sim],e,par)], \
                                      method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,par),rtol=RelTol, atol=AbsTol)
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

#%%Depth&flatness

def depth_2dof_sn(par):
    """ This function returns the value of depth for a given set of parameters
    """
    depth = -(2*math.sqrt(par[4])- (par[6]**2*par[5])/(par[6]**2+par[5]))**3/(6*par[3]**2)
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
    
        parameters : float (list)
    
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
  
#%% old version
num_sim=100 # The expected time to run the simulation for 3 different alphas, num_sim=100 takes approx 100 min.
e=0.5
alpha1=1
alpha2=2
alpha3=5
epsilon= 0.0
parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu, epsilon, omega])

DS = np.linspace(-np.sqrt(2*e/(omega**2+epsilon)),np.sqrt(2*e/(omega**2+epsilon)),12)
num_pts=10 #num_pts is the number of points we want to sample in each x,y directions on the configuration space.
chosen_pts = np.zeros((num_pts,num_pts))
# Record the choice of initial conditions on the configuration space.
for i in range(10):
    V_max= fsolve(lambda x: V_SN2dof(x,DS[1:11][i],parameters)-e,1.5)
    V_x_range= np.linspace(0,V_max,12)
    chosen_pts[i] = V_x_range[1:11]
    
    
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
event_times = np.zeros((num_sim*2,1))
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
event_times = np.zeros((num_sim*2,1))
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
event_times = np.zeros((num_sim*2,1))
for i in range(10):
    print(i)
    for j in range(10):
        event_t = reaction_probability(chosen_pts[i,j],DS[1:11][i],e,parameters,num_sim,react_time_file)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
np.savetxt(react_time_file.name,event_times,fmt='%1.16e')
react_time_file.close()

#%% updated version
# num_pts in sampling_configuration for epsilon=0, alpha=1,2,5is 14,25,45
# num_pts in sampling_configuration for epsilon=5, alpha=1,2,5is 32,50,75
num_sim=100 # The expected time to run the simulation for 3 different alphas, when num_sim=100 is approx 100 min.
e=0.5
alpha1=1
alpha2=2
alpha3=5
epsilon= 0.0
parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu, epsilon, omega])

sampling_pts = sampling_configuration(0,6,-3,3,14,parameters,0.5)
num_chosen_pts = len(sampling_pts)
with open("checke=%s_par=%s_numsim=%s.txt" %(e,parameters[3:],num_sim),'a+') as react_time_file:
    event_times = np.zeros((num_sim*2,1))
    for i in range(num_chosen_pts):
        print(i)
        event_t = reaction_probability(sampling_pts[i,0],sampling_pts[i,1],e,parameters,num_sim)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
    np.savetxt(react_time_file.name,event_times,fmt='%1.16e')

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu, epsilon, omega])
sampling_pts = sampling_configuration(0,6,-3,3,25,parameters,0.5)
num_chosen_pts = len(sampling_pts)
with open("checke=%s_par=%s_numsim=%s.txt" %(e,parameters[3:],num_sim),'a+') as react_time_file:
    event_times = np.zeros((num_sim*2,1))
    for i in range(num_chosen_pts):
        print(i)
        event_t = reaction_probability(sampling_pts[i,0],sampling_pts[i,1],e,parameters,num_sim)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
    np.savetxt(react_time_file.name,event_times,fmt='%1.16e')

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu, epsilon, omega])
sampling_pts = sampling_configuration(0,6,-3,3,45,parameters,0.5)
num_chosen_pts = len(sampling_pts)
with open("checke=%s_par=%s_numsim=%s.txt" %(e,parameters[3:],num_sim),'a+') as react_time_file:
    event_times = np.zeros((num_sim*2,1))
    for i in range(num_chosen_pts):
        print(i)
        event_t = reaction_probability(sampling_pts[i,0],sampling_pts[i,1],e,parameters,num_sim)
        event_times = np.concatenate((event_times,event_t),axis=0, out=None)
    np.savetxt(react_time_file.name,event_times,fmt='%1.16e')
#%% old version
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
event_t1 = event_t1[num_sim*2:]
react_time_file.close()

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu, epsilon, omega])
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
print('Loading the reaction time from data file',react_time_file.name,'\n') 
event_t2 = np.loadtxt(react_time_file.name)
event_t2 = event_t2[num_sim*2:]
react_time_file.close()

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu, epsilon, omega])
react_time_file = open("numpts=%s_e=%s_par=%s_numsim=%s.txt" %(num_pts,e,parameters[3:],num_sim),'a+')
print('Loading the reaction time from data file',react_time_file.name,'\n') 
event_t3 = np.loadtxt(react_time_file.name)
event_t3 = event_t3[num_sim*2:]
react_time_file.close()

#%% updated version
num_sim=100
x=-0.1
y=0.5
e=0.5
alpha1=1
alpha2=2
alpha3=5
epsilon= 0

parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha1, mu, epsilon, omega])  
with open("checke=%s_par=%s_numsim=%s.txt" %(e,parameters[3:],num_sim),'a+') as react_time_file:
    print('Loading the reaction time from data file',react_time_file.name,'\n') 
    event_t1 = np.loadtxt(react_time_file.name)
    event_t1 = event_t1[num_sim*2:]


parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha2, mu, epsilon, omega])
with open("checke=%s_par=%s_numsim=%s.txt" %(e,parameters[3:],num_sim),'a+') as react_time_file:
    print('Loading the reaction time from data file',react_time_file.name,'\n') 
    event_t2 = np.loadtxt(react_time_file.name)
    event_t2 = event_t2[num_sim*2:]


parameters=np.array([MASS_A, MASS_B, EPSILON_S, alpha3, mu, epsilon, omega])
with open("checke=%s_par=%s_numsim=%s.txt" %(e,parameters[3:],num_sim),'a+') as react_time_file:
    print('Loading the reaction time from data file',react_time_file.name,'\n') 
    event_t3 = np.loadtxt(react_time_file.name)
    event_t3 = event_t3[num_sim*2:]


#%% Plot the reaction probability as a function of time
axis_fs=30
num_alp=3 # number of alphas we want to calculate
F = np.zeros((num_alp))
alpha = np.array([1,2,5])
D = depth_2dof_sn(np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega]))

for i in range(num_alp):
    parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i], mu, epsilon, omega])
        
    F[i] = flatness_2dof_sn(-1,10,-5,5,500,parameters)
    print(F[i])
        
sns.set(font_scale = 2)
event_t1 = np.sort(event_t1,axis=None)
event_t2 = np.sort(event_t2,axis=None)
event_t3 = np.sort(event_t3,axis=None)
event_time = [event_t1, event_t2, event_t3]
alpha = [1,2,5]
#ax = plt.gca()
fig, ax1 = plt.subplots()

# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.22, 0.45, 0.2, 0.2]
#left, bottom, width, height = [0.65, 0.2, 0.2, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])
#plot1
for i in range(3):
    event_t = event_time[i]
    non_react = len(np.where(event_t==0)[0])
    react_time=event_t[non_react:]
    react_prob = np.zeros((len(react_time),1))
    for j in range(len(react_time)):
        react_prob[j] = (j+1)/len(event_t)
    ax1.plot(react_time,react_prob, label=r'$\alpha=%s,\mathcal{D}=%.2f,\mathcal{F}=%.2f$'\
            %(alpha[i],D[i],F[i]))

ax1.set_xlabel('$t$', fontsize=axis_fs)
ax1.set_ylabel('reaction probability', fontsize=axis_fs)
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 1)
legend = ax1.legend(loc='best')
#plot2
for i in range(3):
    event_t = event_time[i]
    non_react = len(np.where(event_t==0)[0])
    react_time=event_t[non_react:]
    react_prob = np.zeros((len(react_time),1))
    for j in range(len(react_time)):
        react_prob[j] = (j+1)/len(event_t)
    # extend the reaction probaility graph for uncoupled system at a later time
    if react_time[-1]<15:
        react_prob=np.append(react_prob,react_prob[-1])
        react_time=np.append(react_time,20)
    ax2.plot(react_time,react_prob, label=r'$\alpha=%s,\mathcal{D}=%.2f,\mathcal{F}=%.2f$'\
            %(alpha[i],D[i],F[i]))

#ax2.set_xlabel('$t$', fontsize=axis_fs)
#ax2.set_ylabel('reaction probability', fontsize=axis_fs)
ax2.set_xlim(5, 20)
ax2.set_ylim(0.0, 0.5)
#ax2.set_ylim(0.8, 1)

plt.show()
plt.savefig('reactprob2_alp=1_2_5_mu4_ome3_e=0dot5_ep=%s.pdf'%(epsilon), format='pdf', \

            bbox_inches='tight')
#%% illustration of the choice of the initional conditions on the configuration space and the evolution of the trajectories.
#2dof
resX=100
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=1.0
#mu=alpha**(3/4)
mu = 4
omega=3.0
epsilon=5
axis_fs=30
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega]);
sns.set(font_scale = 2)

def V_SN2dof(x,y,par):
    pot = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2 + 0.5*par[6]*y**2+0.5*par[5]*(x-y)**2
    return pot

xVec = np.linspace(-2,6.2,resX)
yVec = np.linspace(-3,3,resX)
xMat, yMat = np.meshgrid(xVec, yVec)

#fig = plt.figure(figsize = (5,5))
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

sampling_points = sampling_configuration(0,6,-3,3,40,parameters,0.5)
for i in range(len(sampling_points)):
    # plot the sampling points in the configuration space
    ax.scatter(sampling_points[i,0],sampling_points[i,1], s = 5, c = 'b', marker = '.')
TSPAN = [0,10]
RelTol = 3.e-10
AbsTol = 1.e-10 
f1 = lambda t,x: saddlenode2dof(t,x,parameters)
# randomly pick a x and y values on the configuration space.
#pick_x = sampling_points[int(len(sampling_points)-17),0]
#pick_y = sampling_points[int(len(sampling_points)-17),1]
## integrate the Hamilton's equation and plot the trajectory, with negative px
#soln1 = solve_ivp(f1, TSPAN, [pick_x,pick_y,-0.01,-deter_py(pick_x,pick_y,-0.01,0.5,parameters)],\
#                              method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,parameters),rtol=RelTol, atol=AbsTol)
#te1 = soln1.t_events[0]
#coordinates = soln1.y
#ax.plot(coordinates[0,:],coordinates[1,:],c = 'g')

pick_x = sampling_points[-5,0]
pick_y = sampling_points[-5,1]
# integrate the Hamilton's equation and plot the trajectory, with negative px
soln2 = solve_ivp(f1, TSPAN, [pick_x,pick_y,-0.01,-deter_py(pick_x,pick_y,-0.01,0.5,parameters)],\
                              method='RK45',dense_output=True, events = lambda t,x: reaction_event(t,x,parameters),rtol=RelTol, atol=AbsTol)
te1 = soln2.t_events[0]
coordinates2 = soln2.y
ax.plot(coordinates2[0,:],coordinates2[1,:],c = 'blue')

def equ_pos(par):
    """Returns the position of the centre equilibrium point
    """
    xe = 2*np.sqrt(par[4])/par[3] - (par[6]**2*par[5])/(par[3]*(par[6]**2+par[5]))
    ye = xe*par[5]/(par[6]**2+par[5])
    return xe, ye

ax.set_xlim(-2, 6.2)
ax.set_ylim(-3, 3)
ax.scatter(0,0, s = 100, c = 'r', marker = 'x')
ax.scatter(equ_pos(parameters)[0],equ_pos(parameters)[1], s = 100, c = 'r', marker = 'o')
ax.set_xlabel('$x$', fontsize=axis_fs)
ax.set_ylabel('$y$', fontsize=axis_fs)
#ax.set_title(r'$\alpha=%s, \mu=%s$' %(alpha, mu))
plt.show()
plt.savefig('samp_pos_mu=4_omega=3_ep=%s.pdf'%(epsilon), format='pdf', \

            bbox_inches='tight')

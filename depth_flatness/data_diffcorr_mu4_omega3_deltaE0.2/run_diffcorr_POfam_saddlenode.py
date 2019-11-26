# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 10:02:48 2019

@author: Wenyang Lyu and Shibabrat Naik

Script to compute unstable periodic orbits at specified energies for the saddlenode quartic 
Hamiltonian using differential correction
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time
import diffcorr_UPOsHam2dof ### import module xxx where xxx is the name of the python file xxx.py 
import matplotlib as mpl
from matplotlib import cm
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

#% Begin problem specific functions
def init_guess_eqpt_saddlenode(eqNum, par):
    """This function returns the position of the equilibrium points with 
        Saddle (EQNUM=1)
        Centre (EQNUM=2)
    """
    
    if eqNum == 1:
        x0 = [0, 0]
    elif eqNum == 2:
        xe = 2*np.sqrt(parameters[4])/parameters[3] - (parameters[6]**2*parameters[5])/(parameters[3]*(parameters[6]**2+parameters[5]))
        ye = xe*parameters[5]/(parameters[6]**2+parameters[5])
        x0 = [xe,ye]    # EQNUM = 2, stable
    
    return x0

def grad_pot_saddlenode(x, par):
    """This function returns the gradient of the potential energy function V(x,y)
    """     

    dVdx = par[3]*x[0]**2-2*np.sqrt(par[4])*x[0]+par[5]*(x[0]-x[1])
    dVdy = par[6]**2*x[1]-par[5]*(x[0]-x[1])
    
    F = [-dVdx, -dVdy]
    
    return F

def pot_energy_saddlenode(x, y, par):
    """This function returns the potential energy function V(x,y)
    """
    
    return (1/3)*par[3]*x**3 - np.sqrt(par[4])*x**2 + 0.5*par[6]**2*y**2+0.5*par[5]*(x-y)**2


#%
def eigvector_saddlenode(L,par):
    """
     eigenvectors and eigenvalues of the Jacobian evaluated at the equilibrium point, 
     which is the correction of the initial condition.
     check the result obtained from the Jacobian matches the analytic result.
     """

    correcx = par[5]/(par[5] - ( 2*np.sqrt(par[2]) + L**2))
    correcy = par[1]
    
    return correcx, correcy


def guess_lin_saddlenode(eqPt, Ax,L, par):
    """This function returns an initial guess of the UPO.
    """
    
    correcx, correcy = eigvector_saddlenode(L,par)
    

    return [eqPt[0]-Ax*correcx,eqPt[1]-Ax*correcy,0,0]


def jacobian_saddlenode(eqPt, par):
    
    x,y,px,py = eqPt[0:4]
    
    # The first order derivative of the Hamiltonian.
    dVdx = par[3]*x**2-2*np.sqrt(par[4])*x+par[5]*(x-y)
    dVdy = par[6]**2*y-par[5]*(x-y)

    # The following is the Jacobian matrix 
    d2Vdx2 = 2*par[3]*x-2*np.sqrt(par[4])+par[5]
            
    d2Vdy2 = par[6]**2+par[5]
    
    d2Vdydx = -par[5]
        
    d2Vdxdy = d2Vdydx    

    Df = np.array([[  0,     0,    par[0],    0],
                   [0,     0,    0,    par[1]],
                   [-d2Vdx2,  -d2Vdydx,   0,    0],
                   [-d2Vdxdy, -d2Vdy2,    0,    0]])
    
    return Df


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


def diffcorr_setup_saddlenode():
    
    dydot1 = 1
    correcty0= 0
    MAXdydot1 = 1.e-10
    drdot1 = dydot1
    correctr0 = correcty0
    MAXdrdot1 = MAXdydot1
    
    return [drdot1, correctr0, MAXdrdot1]


def conv_coord_saddlenode(x1, y1, dxdot1, dydot1):
    return dydot1


def diffcorr_acc_corr_saddlenode(coords, phi_t1, x0, par):
    """This function computes the correction terms to the initial guess of the UPO
    where correcx0 is the correction in x coodinate.
    Using correction in x or y coordinates depend on the problem.
    The return x0 is now our new guess initional condition of the UPO.
    """
    
    x1, y1, dxdot1, dydot1 = coords
    
    dVdx = par[3]*x1**2-2*np.sqrt(par[4])*x1+par[5]*(x1-y1)
    dVdy = par[6]**2*y1-par[5]*(x1-y1)
    vxdot1 = -dVdx
    vydot1 = -dVdy
    
    #correction to the initial y0
    correcty0 = 1/(phi_t1[3,1] - phi_t1[2,1]*vydot1*(1/vxdot1))*dydot1
    x0[1] = x0[1] - correcty0
    
    return x0


def plot_iter_orbit_saddlenode(x, ax, e, par):
    
#    label_fs = 10
    axis_fs = 15 # fontsize for publications 
    
    ax.plot(x[:,0],x[:,1],x[:,3],'-')
    ax.plot(x[:,0],x[:,1],-x[:,3],'--')
    ax.scatter(x[0,0],x[0,1],x[0,3],s=20,marker='*')
    ax.scatter(x[-1,0],x[-1,1],x[-1,3],s=20,marker='o')
    ax.set_xlabel(r'$x$', fontsize=axis_fs)
    ax.set_ylabel(r'$y$', fontsize=axis_fs)
    ax.set_zlabel(r'$p_y$', fontsize=axis_fs)
    ax.set_title(r'$\Delta E$ = %e' %(np.mean(e) - par[2]) ,fontsize=axis_fs)
    #par(3) is the energy of the saddle
    ax.set_xlim(-0.1, 0.1)
    
    return 


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


#% End problem specific functions


#%% Setting up parameters and global variables

N = 4               # dimension of phase space
MASS_A=1.0
MASS_B=1.0
EPSILON_S=0.0
alpha=0.5
#mu=alpha**(3/4)
mu=4
omega=3
epsilon=0.001

parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega])
eqNum = 1
#model = 'saddlenode'
#eqPt = diffcorr_UPOsHam2dof.get_eq_pts(eqNum,model, parameters)
eqPt = diffcorr_UPOsHam2dof.get_eq_pts(eqNum, init_guess_eqpt_saddlenode, \
                                       grad_pot_saddlenode, parameters)

#energy of the saddle eq pt
eSaddle = diffcorr_UPOsHam2dof.get_total_energy([eqPt[0],eqPt[1],0,0], pot_energy_saddlenode, \
                                                parameters) 
#If orbit is an n-array, e.g. orbit = [orbit_0, orbit_1, ..., orbit_n]
# 
#e = np.zeros(n,1)    
#for i in range(n):
#    e[i] = get_total_energy_deleonberne(orbit[i], parameters)
#e = np.mean(e)

#%%
nFam = 100 # use nFam = 10 for low energy

# first two amplitudes for continuation procedure to get p.o. family
Ax1  = 2.e-5 # initial amplitude (1 of 2) values to use: 2.e-3
Ax2  = 2*Ax1 # initial amplitude (2 of 2)

t = time.time()

#  get the initial conditions and periods for a family of periodic orbits
num_alpha = 10
num_ep = 4
alpha = np.zeros((num_alpha,num_ep))
alpha[:,0] = np.linspace(3,10,num_alpha)
alpha[:,1] = np.linspace(3,10,num_alpha)
alpha[:,2] = np.linspace(3,10,num_alpha)
alpha[:,3] = np.linspace(3,10,num_alpha)
epsilon = np.array([1e-20,1.5,3.0,4.5])
for i in range(4):
    for j in range(10):
        al = alpha[j,i]
        ep = epsilon[i]
        parameters = np.array([MASS_A, MASS_B, EPSILON_S, al, mu, ep, omega])
        print(parameters)
        po_fam_file = open("x0_diffcorr_fam_eqPt%s_alpha%s_ep%s_saddlenode.txt" %(eqNum,al,ep),'a+')
        [po_x0Fam,po_tpFam] = diffcorr_UPOsHam2dof.get_POFam(eqNum, Ax1, Ax2, nFam, \
                                                            po_fam_file, init_guess_eqpt_saddlenode, \
                                                            grad_pot_saddlenode, jacobian_saddlenode, \
                                                            guess_lin_saddlenode, diffcorr_setup_saddlenode, \
                                                            conv_coord_saddlenode, \
                                                            diffcorr_acc_corr_saddlenode, ham2dof_saddlenode, \
                                                            half_period_saddlenode, pot_energy_saddlenode, \
                                                            varEqns_saddlenode, plot_iter_orbit_saddlenode, \
                                                            parameters)  
        
        poFamRuntime = time.time() - t
        x0podata = np.concatenate((po_x0Fam, po_tpFam),axis=1)
        po_fam_file.close()



#%%
# begins with a family of periodic orbits and steps until crossing the
# initial condition with target energy 
# fileName = 'x0po_T_energy_case1_L41.txt'
# fileName = 'x0po_T.txt'

#deltaE_vals = [0.01, 0.1, 1.00, 2.0, 4.0]
#linecolor = ['b','r','g','m','c']
deltaE=0.2
#linecolor = ['b','r']

for i in range(10):
    for j in range(4):
        parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i,j], mu, epsilon[j], omega])
        po_fam_file = open("x0_diffcorr_fam_eqPt%s_alpha%s_ep%s_saddlenode.txt" %(eqNum,alpha[i,j],epsilon[j]) ,'a+')
        eTarget = eSaddle + deltaE 
        print('Loading the periodic orbit family from data file',po_fam_file.name,'\n') 
        x0podata = np.loadtxt(po_fam_file.name)
        po_fam_file.close()
    
    
        #%
        po_brac_file = open("x0po_T_energyPO_eqPt%s_alpha%s_ep%s_brac%s_saddlenode.txt" %(eqNum,alpha[i,j],epsilon[j],deltaE),'a+')
        t = time.time()
        # [x0poTarget,TTarget] = bracket_POEnergy_bp(eTarget, x0podata, po_brac_file)
        x0poTarget,TTarget = diffcorr_UPOsHam2dof.poBracketEnergy(eTarget, x0podata, po_brac_file, \
                                                                  diffcorr_setup_saddlenode, \
                                                                  conv_coord_saddlenode, \
                                                                  diffcorr_acc_corr_saddlenode, \
                                                                  ham2dof_saddlenode, \
                                                                  half_period_saddlenode, \
                                                                  pot_energy_saddlenode, varEqns_saddlenode, \
                                                                  plot_iter_orbit_saddlenode, \
                                                                  parameters)
        poTarE_runtime = time.time()-t
        model_parameters_file = open("model_parameters_eqPt%s_alpha%s_ep%s_DelE%s_saddlenode.txt" %(eqNum,alpha[i,j],epsilon[j],deltaE),'a+')
        np.savetxt(model_parameters_file.name, parameters,fmt='%1.16e')
        model_parameters_file.close()
        po_brac_file.close()
    
    
        #%
        # target specific periodic orbit
        # Target PO of specific energy with high precision does not work for the
        # model 
    
        po_target_file = open("x0_diffcorr_alpha%s_ep%s_deltaE%s_saddlenode.txt" %(alpha[i,j],epsilon[j],deltaE),'a+')
                    
        [x0po, T,energyPO] = diffcorr_UPOsHam2dof.poTargetEnergy(x0poTarget,eTarget, \
                                                                po_target_file, \
                                                                diffcorr_setup_saddlenode, \
                                                                conv_coord_saddlenode, \
                                                                diffcorr_acc_corr_saddlenode, \
                                                                ham2dof_saddlenode, half_period_saddlenode, \
                                                                pot_energy_saddlenode, varEqns_saddlenode, \
                                                                plot_iter_orbit_saddlenode, \
                                                                parameters)
        
        po_target_file.close()


#%% Load periodic orbit data from ascii files
    
x0po = np.zeros((4,10))

for i in range(10):

    po_fam_file = open("x0_diffcorr_alpha%s_ep%s_deltaE%s_saddlenode.txt" %(alpha[i,3],epsilon[3],deltaE),'a+')
    print('Loading the periodic orbit family from data file',po_fam_file.name,'\n') 
    x0podata = np.loadtxt(po_fam_file.name)
    po_fam_file.close()
    x0po[:,i] = x0podata[0:4]


#%% Plotting the unstable periodic orbits at the specified energies

TSPAN = [0,30]
RelTol = 3.e-10
AbsTol = 1.e-10

plt.close('all')
axis_fs = 15


for i in range(10):
    parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha[i,3], mu, epsilon[3], omega])
    f= lambda t,x: ham2dof_saddlenode(t,x,parameters)
    soln = solve_ivp(f, TSPAN, x0po[:,i],method='RK45',dense_output=True, \
                     events = lambda t,x : half_period_saddlenode(t,x,parameters), \
                     rtol=RelTol, atol=AbsTol)
    
    te = soln.t_events[0]
    tt = [0,te[1]] 
    t,x,phi_t1,PHI = diffcorr_UPOsHam2dof.stateTransitMat(tt,x0po[:,i],parameters,varEqns_saddlenode)
    ax = plt.gca(projection='3d')
#    ax.plot(x[:,0],x[:,1],x[:,3],'-',color = linecolor[i], label = '$\Delta E$ = %.2f'%(deltaE))
#    ax.plot(x[:,0],x[:,1],-x[:,3],'-',color = linecolor[i])
    ax.plot(x[:,0],x[:,1],x[:,3],'-', label = '$\Delta E$ = %.2f'%(deltaE))
    ax.plot(x[:,0],x[:,1],-x[:,3],'-')
    ax.scatter(x[0,0],x[0,1],x[0,3],s=20,marker='*')
    ax.scatter(x[0,0],x[0,1],-x[0,3],s=20,marker='o')
    ax.plot(x[:,0], x[:,1], zs=0, zdir='z')


ax = plt.gca(projection='3d')
resX = 100
xVec = np.linspace(-4,4,resX)
yVec = np.linspace(-4,4,resX)
xMat, yMat = np.meshgrid(xVec, yVec)
cset1 = ax.contour(xMat, yMat, \
                   diffcorr_UPOsHam2dof.get_pot_surf_proj(xVec, yVec, pot_energy_saddlenode, \
                                                          parameters), [0.01,0.1,1,2,4], \
                                                          zdir='z', offset=0, 
                                                          linewidths = 1.0, cmap=cm.viridis, \
                                                          alpha = 0.8)

ax.scatter(eqPt[0], eqPt[1], s = 200, c = 'r', marker = 'X')
ax.set_xlabel('$x$', fontsize=axis_fs)
ax.set_ylabel('$y$', fontsize=axis_fs)
ax.set_zlabel('$p_y$', fontsize=axis_fs)
#ax.set_title('$\Delta E$ = %1.e,%1.e,%1.e,%1.e,%1.e' \
#                %(energyPO_1,energyPO_2,energyPO_3,energyPO_4,energyPO_5) ,fontsize=axis_fs)

legend = ax.legend(loc='upper left')
ax.set_xlim(-0.4, 0.25)
ax.set_ylim(-0.5, 0.5)
ax.set_zlim(-1, 1)

plt.grid()
plt.show()

plt.savefig('diffcorr_POfam_saddlenode.pdf',format='pdf',bbox_inches='tight')














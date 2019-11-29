"""
Script to visualize the potential energy surface in the 2 degrees of freedom 
saddle-node Hamiltonian
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pylab import rcParams
    
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

plt.style.use('seaborn') # use sans-serif fonts

mpl.rcParams['axes.spines.left'] = True   ## display axis spines
mpl.rcParams['axes.spines.bottom'] = True
mpl.rcParams['axes.spines.top'] = True
mpl.rcParams['axes.spines.right'] = True
# mpl.rcParams['xtick.top'] = True
# mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0

rcParams['figure.figsize'] = 8, 8 


# axis_fs=15

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal
#m = cm.ScalarMappable(cmap=cm.jet)
mpl.rcParams['font.weight'] = 'normal'

lw = 2.0 # linewidth for line plots
axis_fs=15

# Definition of the 2 DoF PES
def V_SN2dof(x,y,par):
    pot = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2 + \
        0.5*par[6]*y**2+0.5*par[5]*(x-y)**2
    return pot
    

def get_pot_surf_proj(xVec, yVec,par):            

    resX = np.size(xVec)
    resY = np.size(xVec)
    surfProj = np.zeros([resX, resY])
    for i in range(len(xVec)):
        for j in range(len(yVec)):
            surfProj[i,j] = V_SN2dof(xVec[j], yVec[i],par)

    return surfProj



#%% plot of 2 DoF PES surface

resX=100
MASS_A=1
MASS_B=1
EPSILON_S=0
alpha=1
#mu=alpha**(3/4)
mu=0.1
omega=1
epsilon=0.1
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega])
epsilon_c=2*np.sqrt(mu)*omega**2/(omega**2-2*np.sqrt(mu)) # epsilon needs to be smaller that this critical value
print("Critical value of epsilon is %s" %(epsilon_c))


xVec = np.linspace(-2.5,2.5,resX)
yVec = np.linspace(-2.5,2.5,resX)
xMat, yMat = np.meshgrid(xVec, yVec)

fig = plt.figure()
ax = plt.gca(projection = '3d')

surf =ax.plot_surface(xMat, yMat, get_pot_surf_proj(xVec, yVec,parameters), \
                    cmap = 'viridis')

#m.set_array(get_pot_surf_proj(xVec, yVec,parameters))
#ax.set_zlim(0, 1.0)

#legend = ax.legend(loc='best')
# ax.set_title(r'$\alpha=%s, \mu=%s$' %(alpha,mu))
ax.set_xlabel('$x$', fontsize=fal)
ax.set_ylabel('$y$', fontsize=fal)
ax.set_zlabel('$V(x,y)$', fontsize=fal)
ax.set_zlim(0, 10)

# plt.savefig('PES_2DoF_alp=%s_mu=%s_ome=%s_ep=%s.pdf' %(alpha,mu,omega,epsilon),format='pdf',bbox_inches='tight')
plt.savefig('PES_2dof_alp%.3f'%(alpha) + '_mu%.3f'%(mu) + \
    '_ome%.3f'%(omega) + '_eps%.3f'%(epsilon) + '.pdf', bbox_inches = 'tight')

plt.draw()
plt.pause(0.01)
plt.show() # block = False





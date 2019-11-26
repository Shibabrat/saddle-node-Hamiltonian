# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 16:32:46 2019

@author: wl16298
"""

#%% plot of 1 DoF PES with alpha and mu, keep alpha as a constant

from functools import partial
%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
x = np.linspace(-5, 10, 1000)
axis_fs=15
# Definition of the PES
def V(alpha,mu,x):
    pot = (1/3)*alpha*x**3 - math.sqrt(mu)*x**2
    return pot

ax = plt.gca()


plot1 = ax.plot(x,V(1,0.1,x),label=r'$\alpha=1,\mu=0.1$')

plot2 = ax.plot(x,V(1,0.5,x),label=r'$\alpha=1,\mu=0.5$')
 
plot3 = ax.plot(x,V(1,1.0,x),label=r'$\alpha=1,\mu=1.0$')

plot4 = ax.plot(x,V(1,2.0,x),label=r'$\alpha=1,\mu=2.0$')
 
ax.scatter(0, 0, s = 100, c = 'r', marker = 'X')
#ax.scatter(0, 1, s = 100, c = 'r', marker = 'X')
ax.set_xlabel('$x$', fontsize=axis_fs)
ax.set_ylabel('$V(x)$', fontsize=axis_fs)

legend = ax.legend(loc='best')

ax.set_xlim(-5, 15)
ax.set_ylim(-5, 5)
plt.grid()
plt.show()


#%% plot of 1 DoF PES with alpha and mu, keep mu as a constant

from functools import partial
%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
x = np.linspace(-5, 10, 1000)
axis_fs=15

# Definition of the PES
def V(alpha,mu,x):
    pot = (1/3)*alpha*x**3 - math.sqrt(mu)*x**2
    return pot

ax = plt.gca()


plot1 = ax.plot(x,V(0.5,1,x),label=r'$\alpha=0.5,\mu=1$')

plot2 = ax.plot(x,V(1.0,1,x),label=r'$\alpha=1.0,\mu=1$')
 
plot3 = ax.plot(x,V(1.5,1.0,x),label=r'$\alpha=1.5,\mu=1$')

plot4 = ax.plot(x,V(2,1.0,x),label=r'$\alpha=2.0,\mu=1$')
 
ax.scatter(0, 0, s = 100, c = 'r', marker = 'X')
#ax.scatter(0, 1, s = 100, c = 'r', marker = 'X')
ax.set_xlabel('$x$', fontsize=axis_fs)
ax.set_ylabel('$V(x)$', fontsize=axis_fs)

legend = ax.legend(loc='best')

ax.set_xlim(-5, 7)
ax.set_ylim(-10, 5)
plt.grid()
plt.show()

#%% plot of 2 DoF PES surface
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import math
from matplotlib import cm

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal
#m = cm.ScalarMappable(cmap=cm.jet)
axis_fs=15
resX=100
MASS_A=1
MASS_B=1
EPSILON_S=0
alpha=1
#mu=alpha**(3/4)
mu=0.1
omega=1
epsilon=0.1
parameters = np.array([MASS_A, MASS_B, EPSILON_S, alpha, mu, epsilon, omega]);
epsilon_c=2*np.sqrt(mu)*omega**2/(omega**2-2*np.sqrt(mu)) # epsilon needs to be smaller that this critical value
print("Critical value of epsilon is %s" %(epsilon_c))

# Definition of the 2 DoF PES
def V_SN2dof(x,y,par):
    pot = (1/3)*par[3]*x**3 - math.sqrt(par[4])*x**2 + 0.5*par[6]*y**2+0.5*par[5]*(x-y)**2
    return pot
    

def get_pot_surf_proj(xVec, yVec,par):            

    resX = np.size(xVec)
    resY = np.size(xVec)
    surfProj = np.zeros([resX, resY])
    for i in range(len(xVec)):
        for j in range(len(yVec)):
            surfProj[i,j] = V_SN2dof(xVec[j], yVec[i],par)

    return surfProj

xVec = np.linspace(-2.5,2.5,resX)
yVec = np.linspace(-2.5,2.5,resX)
xMat, yMat = np.meshgrid(xVec, yVec)
fig = plt.figure()
ax = plt.gca(projection = '3d')


surf =ax.plot_surface(xMat, yMat, get_pot_surf_proj(xVec, yVec,parameters))

#m.set_array(get_pot_surf_proj(xVec, yVec,parameters))
#ax.set_zlim(0, 1.0)

#legend = ax.legend(loc='best')
ax.set_title(r'$\alpha=%s, \mu=%s$' %(alpha,mu))
ax.set_xlabel('$x$', fontsize=axis_fs)
ax.set_ylabel('$y$', fontsize=axis_fs)
ax.set_zlabel('$V(x,y)$', fontsize=axis_fs)
ax.set_zlim(0, 10)

plt.savefig('PES_2DoF_alp=%s_mu=%s_ome=%s_ep=%s.pdf' %(alpha,mu,omega,epsilon),format='pdf',bbox_inches='tight')
plt.show()


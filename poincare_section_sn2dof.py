# Script to set-up Poincare section computation for the 2 dof saddle-node 
# Hamiltonian

import numpy as np
from scipy.integrate import solve_ivp
from functools import partial
import matplotlib.pyplot as plt
import matplotlib as mpl

from pylab import rcParams
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal


# Import system-bath module
import saddlenode2dof
import importlib
importlib.reload(saddlenode2dof)
import saddlenode2dof as sn2dof

# low tolerance
# rtol_val = 1e-8
# atol_val = 1e-10

# high tolerance
rtol_val = 1e-12
atol_val = 1e-14

# Model parameters
N = 2
MASS_A = 1.0
MASS_B = 1.0
MU = 4.00
ALPHA = 1.00
OMEGA = 3.00
EPSILON = 0.00

deltaEnergy = 0.05
saddleEnergy = 0
totalEnergy = saddleEnergy + deltaEnergy
printFlag = True
timespan = [0, 20]
if EPSILON > 0:
    sos_crossings_filename = 'coupled_sn2dof_alpha' + str(int(ALPHA)) \
                                + '_tau' + str(timespan[1]) + '.txt'
else:
    sos_crossings_filename = 'uncoupled_sn2dof_alpha' + str(int(ALPHA)) \
                                + '_tau' + str(timespan[1]) + '.txt'


params = [MASS_A, MASS_B, MU, ALPHA, OMEGA, EPSILON]


# Generate initial conditions on the constant energy surface
# py = sn2dof.momentum_fixed_energy(4, 2, 0.1, params, deltaEnergy)
# print(py)

energy_boundary = sn2dof.energysurface_intersect_sos(params, deltaEnergy, 1.0)


# Generate initial conditions on the surface of section 
numpts = 50
yMax = np.max(energy_boundary[:,0])
yMin = np.min(energy_boundary[:,0])

pyMax = np.max(energy_boundary[:,1])
pyMin = np.min(energy_boundary[:,1])

yGrid = np.linspace(yMin, yMax, numpts)
pyGrid = np.linspace(pyMin, pyMax, numpts)

[yMesh, pyMesh] = np.meshgrid(yGrid, pyGrid)

# Save points that on the energy surface by calculating the momentum coordinate
x_e = (1/ALPHA)*(2.0*np.sqrt(MU) \
        - (OMEGA**2*EPSILON)/(OMEGA**2 + EPSILON))
# xGrid = x_e*np.ones((numpts,))

initconds = []
for y in yGrid:
    for py in pyGrid:
        px = sn2dof.momentum_fixed_energy(x_e, y, py, params, totalEnergy)

        if ~np.isnan(px):
            initconds.append([x_e, y, px, py])
        else:
            continue

initconds = np.array(initconds)


intersections_sos = []
for i in range(np.size(initconds,0)):

    odes_func = partial(sn2dof.vector_field, params)
    cross_sos_func = partial(sn2dof.cross_sos_center_eqpt, params)
    cross_xneg5_func = partial(sn2dof.cross_xneg5, params)

    cross_sos_func.terminal = False
    cross_sos_func.direction = 0

    cross_xneg5_func.terminal = True
    cross_xneg5_func.direction = 0

    sol = solve_ivp(odes_func, timespan, initconds[i,:], method='RK45', \
                    dense_output = True, \
                    events = (cross_sos_func, cross_xneg5_func) , \
                    rtol=rtol_val, atol=atol_val, vectorized = True)

    # Now store the intersections of the trajectory that has the correct direction of crossing: p_x < 0
    intersections_sos.append(sol.sol(sol.t_events[0][1:]).T) # builds a nested list


    print(str(i) + ' initial conditions done.')    


## Unpack the list into an array
num_intersections = 0
intersections = np.array([])
for i in range(len(intersections_sos)):
    intersections = np.append(intersections, intersections_sos[i])
    num_intersections += np.size(intersections_sos[i],0)
        
intersections = np.reshape(intersections, \
    (int(len(intersections)/(2*N)), int(2*N)))


# Save the intersection to file
np.savetxt(sos_crossings_filename, intersections)

index_pos_px = np.where(intersections[:,2] > 0)
if printFlag:
    plt.close('all')
    fig = plt.figure(figsize=(6, 6))
    ax = fig.gca()

    ax.plot(energy_boundary[:,0], energy_boundary[:,1], linewidth = 2, color = 'm')
    # ax.scatter(initconds[:,1], initconds[:,3], s = 5)
    ax.scatter(intersections[index_pos_px,1], intersections[index_pos_px,3], s = 0.5, c = 'k')
    # for i in range(len(intersections_sos)):
    #     ax.scatter(intersections_sos[i][1:][:,1], \
    #         intersections_sos[i][1:][:,3], s = 1, c = 'k')

    ax.set_xlabel(r'$y$', fontsize = 25)
    ax.set_ylabel(r'$p_y$', fontsize = 25)
    # plt.tick_params(axis = 'both', which = 'major', labelsize = 15)
    plt.show()



    
    

    
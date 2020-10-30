#!/usr/bin/python

import sys
import os
import numpy as np
import math
import h5py
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import rc,rcParams
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.cm as cm
plt.rc('text', usetex = True )
rc('font', family='sanserif')
plt.rc('font', size = 20)
plt.rcParams['text.latex.preamble'] = [ r'\usepackage{amsmath}', r'\usepackage{amssymb}', r'\usepackage{color}']


### Open data file
data_file = h5py.File('2d_riemann_problem_383.h5','r')
# Read data
x_data     = data_file['x'][0,:,:];     x_data     = np.asarray( x_data.flatten() )
y_data     = data_file['y'][0,:,:];     y_data     = np.asarray( y_data.flatten() )
z_data     = data_file['z'][0,:,:];     z_data     = np.asarray( z_data.flatten() )
rho_data   = data_file['rho'][0,:,:];   rho_data   = np.asarray( rho_data.flatten() )
u_data     = data_file['u'][0,:,:];     u_data     = np.asarray( u_data.flatten() )
v_data     = data_file['v'][0,:,:];     v_data     = np.asarray( v_data.flatten() )
w_data     = data_file['w'][0,:,:];     w_data     = np.asarray( w_data.flatten() )
E_data     = data_file['E'][0,:,:];     E_data     = np.asarray( E_data.flatten() )
P_data     = data_file['P'][0,:,:];     P_data     = np.asarray( P_data.flatten() )
T_data     = data_file['T'][0,:,:];     T_data     = np.asarray( T_data.flatten() )
sos_data   = data_file['sos'][0,:,:];   sos_data   = np.asarray( sos_data.flatten() )
mu_data    = data_file['mu'][0,:,:];    mu_data    = np.asarray( mu_data.flatten() )
kappa_data = data_file['kappa'][0,:,:]; kappa_data = np.asarray( kappa_data.flatten() )
ke_data    = 0.5*( u_data*u_data + v_data*v_data + w_data*w_data ) 
e_data     = E_data - ke_data


########## PARULA COLORMAP ##########
from matplotlib.colors import LinearSegmentedColormap

cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
 [0.0589714286, 0.6837571429, 0.7253857143], 
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
 [0.7184095238, 0.7411333333, 0.3904761905], 
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
 [0.9763, 0.9831, 0.0538]]

parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
########## PARULA COLORMAP ##########


### Plot density vs. x- & y-direction

# Clear plot
plt.clf()

# Plot data
my_cmap = parula_map
my_norm = colors.Normalize( vmin = rho_data.min(), vmax = rho_data.max() )
plt.tricontour( x_data, y_data, rho_data, norm = my_norm, levels = np.arange( 0.0, 2.01, 7.5e-2 ), linewidths = 0.35, colors = 'black' )
cs = plt.tricontourf( x_data, y_data, rho_data, cmap = my_cmap, norm = my_norm, levels = np.arange( 0.0, 2.01, 1.0e-3 ) )

# Colorbar
cbar = plt.colorbar( cs, shrink = 0.95, pad = 0.02, ticks = [ 0.0, 0.4, 0.8, 1.2, 1.6, 2.0 ] )
cbar.ax.tick_params( labelsize = 16 ) 
plt.text( 0.95, 1.04, r'$\rho \thinspace \textrm{[kg/m}^\textrm{3}\textrm{]}$', fontsize = 16 )
plt.clim( rho_data.min(), rho_data.max() )

## Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.01, 0.2 ) )
plt.tick_params( axis = 'x', left = True, right = True, top = True, bottom = True, direction = 'inout', labelsize = 16 )
plt.ylim( 0.0, 1.0 )
plt.yticks( np.arange( 0.0, 1.01, 0.2 ) )
plt.tick_params( axis = 'y', left = True, right = True, top = True, bottom = True, direction = 'inout', labelsize = 16 )
plt.gca().set_aspect( 'equal', adjustable = 'box' )
ax = plt.gca()
ax.tick_params( axis = 'both', pad = 5.0 )
plt.xlabel( r'$x \thinspace \textrm{[m]}$' )
plt.ylabel( r'$y \thinspace \textrm{[m]}$' )
plt.savefig('density_vs_x_y_direction.eps', format = 'eps', bbox_inches = 'tight' )
plt.savefig('density_vs_x_y_direction.png', format = 'png', bbox_inches = 'tight', dpi = 600 )


### Plot pressure vs. x- & y-direction

# Clear plot
plt.clf()

# Plot data
my_cmap = parula_map
my_norm = colors.Normalize( vmin = P_data.min(), vmax = P_data.max() )
plt.tricontour( x_data, y_data, P_data, norm = my_norm, levels = np.arange( 0.0, 2.01, 7.5e-2 ), linewidths = 0.35, colors = 'black' )
cs = plt.tricontourf( x_data, y_data, P_data, cmap = my_cmap, norm = my_norm, levels = np.arange( 0.0, 2.01, 1.0e-3 ) )

# Colorbar
cbar = plt.colorbar( cs, shrink = 0.95, pad = 0.02, ticks = [ 0.0, 0.4, 0.8, 1.2, 1.6, 2.0 ] )
cbar.ax.tick_params( labelsize = 16 ) 
plt.text( 0.975, 1.04, r'$P \thinspace \textrm{[Pa]}$', fontsize = 16 )
plt.clim( P_data.min(), P_data.max() )

## Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.01, 0.2 ) )
plt.tick_params( axis = 'x', left = True, right = True, top = True, bottom = True, direction = 'inout', labelsize = 16 )
plt.ylim( 0.0, 1.0 )
plt.yticks( np.arange( 0.0, 1.01, 0.2 ) )
plt.tick_params( axis = 'y', left = True, right = True, top = True, bottom = True, direction = 'inout', labelsize = 16 )
plt.gca().set_aspect( 'equal', adjustable = 'box' )
ax = plt.gca()
ax.tick_params( axis = 'both', pad = 5.0 )
plt.xlabel( r'$x \thinspace \textrm{[m]}$' )
plt.ylabel( r'$y \thinspace \textrm{[m]}$' )
plt.savefig('pressure_vs_x_y_direction.eps', format = 'eps', bbox_inches = 'tight' )
plt.savefig('pressure_vs_x_y_direction.png', format = 'png', bbox_inches = 'tight', dpi = 600 )

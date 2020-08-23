#!/usr/bin/python

import sys
import os
import numpy as np
import h5py    
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib import rc,rcParams
from scipy.optimize import fsolve
plt.rc( 'text', usetex = True )
rc('font', family='sanserif')
plt.rc( 'font', size = 20 )
plt.rcParams['text.latex.preamble'] = [ r'\usepackage{amsmath}', r'\usepackage{amssymb}', r'\usepackage{color}' ]


### Open data file
data_file = h5py.File( '1d_sod_shock_tube_60.h5', 'r' )
#list( data_file.keys() )
x_data     = data_file['x'][0,0,:];     x_data     = np.asarray( x_data.flatten() )
y_data     = data_file['y'][0,0,:];     y_data     = np.asarray( y_data.flatten() )
z_data     = data_file['z'][0,0,:];     z_data     = np.asarray( z_data.flatten() )
rho_data   = data_file['rho'][0,0,:];   rho_data   = np.asarray( rho_data.flatten() )
u_data     = data_file['u'][0,0,:];     u_data     = np.asarray( u_data.flatten() )
v_data     = data_file['v'][0,0,:];     v_data     = np.asarray( v_data.flatten() )
w_data     = data_file['w'][0,0,:];     w_data     = np.asarray( w_data.flatten() )
E_data     = data_file['E'][0,0,:];     E_data     = np.asarray( E_data.flatten() )
P_data     = data_file['P'][0,0,:];     P_data     = np.asarray( P_data.flatten() )
T_data     = data_file['T'][0,0,:];     T_data     = np.asarray( T_data.flatten() )
sos_data   = data_file['sos'][0,0,:];   sos_data   = np.asarray( sos_data.flatten() )
mu_data    = data_file['mu'][0,0,:];    mu_data    = np.asarray( mu_data.flatten() )
kappa_data = data_file['kappa'][0,0,:]; kappa_data = np.asarray( kappa_data.flatten() )
ke_data    = 0.5*( u_data*u_data + v_data*v_data + w_data*w_data ) 
e_data     = E_data - ke_data


### Plot density vs. x-direction

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( x_data, rho_data, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.1, 0.2 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$x \thinspace \textrm{[m]}$' )
plt.ylim( 0.0, 1.1 )
plt.yticks( np.arange( 0.0, 1.2, 0.2 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$\rho \thinspace \textrm{[kg/m}^\textrm{3}\textrm{]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'density_vs_x_direction.eps', format = 'eps', bbox_inches = 'tight' )


### Plot u-velocity vs. x-direction

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( x_data, u_data, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.1, 0.2 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$x \thinspace \textrm{[m]}$' )
plt.ylim( -0.1, 1.6 )
plt.yticks( np.arange( 0.0, 1.7, 0.4 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$u \thinspace \textrm{[m/s]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='lower left' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'u_velocity_vs_x_direction.eps', format = 'eps', bbox_inches = 'tight' )


### Plot pressure vs. x-direction

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( x_data, P_data, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.1, 0.2 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$x \thinspace \textrm{[m]}$' )
plt.ylim( 0.0, 1.1 )
plt.yticks( np.arange( 0.0, 1.2, 0.2 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$P \thinspace \textrm{[Pa]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'pressure_vs_x_direction.eps', format = 'eps', bbox_inches = 'tight' )


### Plot internal energy vs. x-direction

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( x_data, e_data, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.1, 0.2 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$x \thinspace \textrm{[m]}$' )
plt.ylim( 1.8, 3.8 )
plt.yticks( np.arange( 1.8, 4.0, 0.4 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$e \thinspace \textrm{[J/kg]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper left' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'internal_energy_vs_x_direction.eps', format = 'eps', bbox_inches = 'tight' )

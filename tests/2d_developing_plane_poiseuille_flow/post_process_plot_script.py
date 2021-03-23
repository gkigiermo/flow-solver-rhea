#!/usr/bin/python

import sys
import os
import numpy as np
import h5py    
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib import rc,rcParams
#np.set_printoptions(threshold=sys.maxsize)
plt.rc( 'text', usetex = True )
rc('font', family='sanserif')
plt.rc( 'font', size = 20 )
plt.rcParams['text.latex.preamble'] = [ r'\usepackage{amsmath}', r'\usepackage{amssymb}', r'\usepackage{color}' ]


### Open data file
data_file = h5py.File( '2d_developing_plane_poiseuille_flow_427447.h5', 'r' )
# Read data
x_index    = 195
x_data     = data_file['x'][1,:,x_index];     x_data     = np.asarray( x_data.flatten() )
y_data     = data_file['y'][1,:,x_index];     y_data     = np.asarray( y_data.flatten() )
z_data     = data_file['z'][1,:,x_index];     z_data     = np.asarray( z_data.flatten() )
rho_data   = data_file['rho'][1,:,x_index];   rho_data   = np.asarray( rho_data.flatten() )
u_data     = data_file['u'][1,:,x_index];     u_data     = np.asarray( u_data.flatten() )
v_data     = data_file['v'][1,:,x_index];     v_data     = np.asarray( v_data.flatten() )
w_data     = data_file['w'][1,:,x_index];     w_data     = np.asarray( w_data.flatten() )
E_data     = data_file['E'][1,:,x_index];     E_data     = np.asarray( E_data.flatten() )
P_data     = data_file['P'][1,:,x_index];     P_data     = np.asarray( P_data.flatten() )
T_data     = data_file['T'][1,:,x_index];     T_data     = np.asarray( T_data.flatten() )
sos_data   = data_file['sos'][1,:,x_index];   sos_data   = np.asarray( sos_data.flatten() )
mu_data    = data_file['mu'][1,:,x_index];    mu_data    = np.asarray( mu_data.flatten() )
kappa_data = data_file['kappa'][1,:,x_index]; kappa_data = np.asarray( kappa_data.flatten() )
ke_data    = 0.5*( u_data*u_data + v_data*v_data + w_data*w_data )
e_data     = E_data - ke_data


### Theoretical plane Poiseuille velocity profile
delta        = 1.0          # Channel half-height [m]
u_tau        = 1.0          # Friction velocity [m/s]
rho_0        = 1.0          # Reference density [kg/m3]
mu           = 0.033333     # Dynamic viscosity [Pa s]
nu           = mu/rho_0     # Kinematic viscosity [m2/s]
G            = 1.0          # Constant pressure gradient [Pa/m]
h            = 2.0*delta    # Channel height [m]
num_y_points = 100          # Number of points in wall-normal direction

y_points  = np.linspace( 0.0, h, num = num_y_points )
u_profile = ( G/( 2.0*mu ) )*y_points*( h - y_points )
u_b       = np.sum( u_profile )/num_y_points
#print( u_b )


### Plot u+ vs. y+

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( y_data*u_tau/nu, u_data/u_tau, marker = 'p', s = 50, color = 'firebrick', zorder = 1, label = r'$\textrm{RHEA}$' )
plt.plot( y_points*u_tau/nu, u_profile/u_tau, linestyle = '-', linewidth = 1, color = 'black', zorder = 0, label = r'$\textrm{Poiseuille flow (periodic)}$' )

# Configure plot
plt.xlim( 0.0, 60.0 )
plt.xticks( np.arange( 0.0, 60.01, 15.0 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$y^{+}$' )
plt.ylim( 0.0, 20.0 )
plt.yticks( np.arange( 0.0, 20.1, 5.0 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$u^{+}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper left' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'u_plus_vs_y_plus.eps', format = 'eps', bbox_inches = 'tight' )

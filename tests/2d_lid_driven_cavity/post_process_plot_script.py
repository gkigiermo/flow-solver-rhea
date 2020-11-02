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
data_file = h5py.File( '2d_lid_driven_cavity_10000.h5', 'r' )
#list( data_file.keys() )
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


### Open reference solution files
y_ghia, u_ghia = np.loadtxt( 'ghia_u_velocity_solution.csv', delimiter=',', unpack = 'True' )
x_ghia, v_ghia = np.loadtxt( 'ghia_v_velocity_solution.csv', delimiter=',', unpack = 'True' )


### Reference parameters
L       = 0.1							# Size of cavity [m]
U_lid   = 100.0   						# Lid velocity [m/s]
epsilon = 0.5*L/( int( np.sqrt( len( y_data ) ) - 2.0 ) )	# Geometric epsilon based on grid resolution


### Plot u-velocity vs. y-direction (vertical line through geometric center of cavity)

# Prepare data along vertical line through geometric center of cavity
y_center_data = np.zeros( int( np.sqrt( len( y_data ) ) ) )
u_center_data = np.zeros( int( np.sqrt( len( u_data ) ) ) )
counter = 0
y_old = -1.0
for p in range( 0, len( x_data ) ):
    if( np.abs( x_data[p] - ( 0.5*L - 1.0e-1*epsilon ) ) < epsilon ):
        if( y_data[p] > y_old ):
            #print( x_data[p], y_data[p] )
            y_center_data[counter] = y_data[p]
            u_center_data[counter] = u_data[p]
            counter += 1
            y_old = y_data[p] 

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( y_ghia, u_ghia, marker = 'p', s = 50, color = 'black', label = r'$\textrm{Ghia et al.}$' )
plt.plot( y_center_data/L, u_center_data/U_lid, linestyle = '--', color = 'firebrick', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.1, 0.2 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$y/L$' )
plt.ylim( -0.5, 1.0 )
plt.yticks( np.arange( -0.5, 1.1, 0.3 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$u/U_\textrm{lid}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper left' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'u_velocity_vs_y_direction.eps', format = 'eps', bbox_inches = 'tight' )


### Plot v-velocity vs. x-direction (horizontal line through geometric center of cavity)

# Prepare data along vertical line through geometric center of cavity
x_center_data = np.zeros( int( np.sqrt( len( x_data ) ) ) )
v_center_data = np.zeros( int( np.sqrt( len( v_data ) ) ) )
counter = 0
x_old = -1.0
for p in range( 0, len( y_data ) ):
    if( np.abs( y_data[p] - ( 0.5*L - 1.0e-1*epsilon ) ) < epsilon ):
        if( x_data[p] > x_old ):
            #print( x_data[p], y_data[p] )
            x_center_data[counter] = x_data[p]
            v_center_data[counter] = v_data[p]
            counter += 1
            x_old = x_data[p] 

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( x_ghia, v_ghia, marker = 'p', s = 50, color = 'black', label = r'$\textrm{Ghia et al.}$' )
plt.plot( x_center_data/L, v_center_data/U_lid, linestyle = '--', color = 'firebrick', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 0.0, 1.0 )
plt.xticks( np.arange( 0.0, 1.1, 0.2 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
#plt.xscale( 'log' )
plt.xlabel( r'$x/L$' )
plt.ylim( -0.6, 0.6 )
plt.yticks( np.arange( -0.6, 0.61, 0.3 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$v/U_\textrm{lid}$' )
#legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'v_velocity_vs_x_direction.eps', format = 'eps', bbox_inches = 'tight' )

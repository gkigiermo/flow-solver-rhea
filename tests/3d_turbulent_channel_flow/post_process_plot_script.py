#!/usr/bin/python

import sys
import os
import numpy as np
import h5py    
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib import rc,rcParams
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.cm as cm
#np.set_printoptions(threshold=sys.maxsize)
plt.rc( 'text', usetex = True )
rc('font', family='sanserif')
plt.rc( 'font', size = 20 )
plt.rcParams['text.latex.preamble'] = [ r'\usepackage{amsmath}', r'\usepackage{amssymb}', r'\usepackage{color}' ]


### Open data file
data_file = h5py.File( '3d_turbulent_channel_flow_2100000.h5', 'r' )
#list( data_file.keys() )
y_data         = data_file['y'][:,:,:]
avg_u_data  = data_file['avg_u'][:,:,:]
avg_v_data  = data_file['avg_v'][:,:,:]
avg_w_data  = data_file['avg_w'][:,:,:]
rmsf_u_data = data_file['rmsf_u'][:,:,:]
rmsf_v_data = data_file['rmsf_v'][:,:,:]
rmsf_w_data = data_file['rmsf_w'][:,:,:]
num_points_x  = avg_u_data[0,0,:].size
num_points_y  = avg_u_data[0,:,0].size
num_points_z  = avg_u_data[:,0,0].size
num_points_xz = num_points_x*num_points_z


### Open reference solution file
y_plus_ref, u_plus_ref, rmsf_uu_plus_ref, rmsf_vv_plus_ref, rmsf_ww_plus_ref = np.loadtxt( 'reference_solution.csv', delimiter=',', unpack = 'True' )


### Reference parameters
rho_0  = 1.0				# Reference density [kg/m3]
u_tau  = 1.0				# Friction velocity [m/s]
delta  = 1.0				# Channel half-height [m]
Re_tau = 180.0				# Friction Reynolds number [-]
mu_ref = rho_0*u_tau*delta/Re_tau	# Dynamic viscosity [Pa s]
nu_ref = mu_ref/rho_0			# Kinematic viscosity [m2/s]


### Allocate averaged variables
avg_y_plus  = np.zeros( int( 0.5*num_points_y ) )
avg_u_plus  = np.zeros( int( 0.5*num_points_y ) )
rmsf_u_plus = np.zeros( int( 0.5*num_points_y ) )
rmsf_v_plus = np.zeros( int( 0.5*num_points_y ) )
rmsf_w_plus = np.zeros( int( 0.5*num_points_y ) )


### Average variables in space
for j in range( 0, num_points_y ):
    for i in range( 0, num_points_x ):
        for k in range( 0, num_points_z ):
            aux_j = j
            if( j > ( int( 0.5*num_points_y ) - 1 ) ):
                aux_j = num_points_y - j - 1
            avg_y_plus[aux_j]  += ( 0.5/num_points_xz )*y_data[k,aux_j,i]*( u_tau/nu_ref )
            avg_u_plus[aux_j]  += ( 0.5/num_points_xz )*avg_u_data[k,j,i]*( 1.0/u_tau )
            rmsf_u_plus[aux_j] += ( 0.5/num_points_xz )*rmsf_u_data[k,j,i]*( 1.0/u_tau )
            rmsf_v_plus[aux_j] += ( 0.5/num_points_xz )*rmsf_v_data[k,j,i]*( 1.0/u_tau )
            rmsf_w_plus[aux_j] += ( 0.5/num_points_xz )*rmsf_w_data[k,j,i]*( 1.0/u_tau )
#print( avg_y_plus )
#print( avg_u_plus )
#print( rmsf_u_plus )
#print( rmsf_v_plus )
#print( rmsf_w_plus )


### Plot u+ vs. y+

# Clear plot
plt.clf()

# Read & Plot data
plt.plot( y_plus_ref, u_plus_ref, linestyle = '-', linewidth = 1, color = 'black', zorder = 0, label = r'$\textrm{Moser et al., }Re_\tau = 180$' )
plt.scatter( avg_y_plus, avg_u_plus, marker = 'p', s = 50, color = 'firebrick', zorder = 1, label = r'$\textrm{RHEA}$' )

# Configure plot
plt.xlim( 1.0e-1, 2.0e2 )
plt.xticks( np.arange( 1.0e-1, 2.01e2, 1.0 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
plt.xscale( 'log' )
plt.xlabel( r'$y^{+}$' )
plt.ylim( 0.0, 20.0 )
plt.yticks( np.arange( 0.0, 20.1, 5.0 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$u^{+}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper left' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'u_plus_vs_y_plus.eps', format = 'eps', bbox_inches = 'tight' )


### Plot u_rms+, v_rms+, w-rms+ vs. y+

# Clear plot
plt.clf()

# Read & Plot data
plt.plot( y_plus_ref, rmsf_uu_plus_ref, linestyle = '-', linewidth = 1, color = 'black', zorder = 0 )
plt.plot( y_plus_ref, rmsf_vv_plus_ref, linestyle = '-', linewidth = 1, color = 'black', zorder = 0 )
plt.plot( y_plus_ref, rmsf_ww_plus_ref, linestyle = '-', linewidth = 1, color = 'black', zorder = 0 )
plt.scatter( avg_y_plus, rmsf_u_plus, marker = 'p', s = 50, color = 'firebrick', zorder = 1 )
plt.scatter( avg_y_plus, rmsf_v_plus, marker = 'p', s = 50, color = 'firebrick', zorder = 1 )
plt.scatter( avg_y_plus, rmsf_w_plus, marker = 'p', s = 50, color = 'firebrick', zorder = 1 )

# Configure plot
plt.xlim( 1.0e-1, 2.0e2 )
plt.xticks( np.arange( 1.0e-1, 2.01e2, 1.0 ) )
plt.tick_params( axis = 'x', bottom = True, top = True, labelbottom = 'True', labeltop = 'False', direction = 'in' )
plt.xscale( 'log' )
plt.xlabel( r'$y^{+}$' )
plt.ylim( 0.0, 3.0 )
plt.yticks( np.arange( 0.0, 3.1, 0.5 ) )
plt.tick_params( axis = 'y', left = True, right = True, labelleft = 'True', labelright = 'False', direction = 'in' )
#plt.yscale( 'log' )
plt.ylabel( r'$u_{\textrm{rms}}^{+}, v_{\textrm{rms}}^{+}, w_{\textrm{rms}}^{+}$' )
#legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper left' )
plt.text( 1.05, 1.0, r'$u_{\textrm{rms}}^{+}$' )
plt.text( 17.5, 0.2, r'$v_{\textrm{rms}}^{+}$' )
plt.text( 4.00, 0.9, r'$w_{\textrm{rms}}^{+}$' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'uvw_rms_plus_vs_y_plus.eps', format = 'eps', bbox_inches = 'tight' )

#!/usr/bin/python

import sys
import os
import numpy as np
import math
import h5py    


### Set fixed parameters
L       = 2.0*np.pi			                # Reference length [m]	
rho_0   = 1.0				                # Reference density [kg/m3]	
U_0     = 1.0				                # Reference velocity [m/s]
nu_ref  = 1.0				                # Kinematic viscosity [m2/s]
gamma_0 = 1.4				                # Ratio of heat capacities [-]
Ma      = 1.0e-1/np.sqrt( gamma_0 )	                # Mach number [-]
P_0     = rho_0*U_0*U_0/( gamma_0*Ma*Ma )               # Reference pressure [Pa]
time    = rho_0*L*L/( 8.0*nu_ref*np.pi*np.pi ) 		# Final time [s]


### Open data file
data_file = h5py.File( '2d_taylor_green_vortex_223.h5', 'r' )
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

### Calculate L2-norm relative errors between discrete and analytical solution
u_L2_norm_error = 0.0
u_L2_norm_den   = 0.0
v_L2_norm_error = 0.0
v_L2_norm_den   = 0.0
for i in range( 0, len( x_data ) ):
    if( ( x_data[i] > 0.0 ) and ( x_data[i] < ( 2.0*np.pi ) ) ):
        if( ( y_data[i] > 0.0 ) and ( y_data[i] < ( 2.0*np.pi ) ) ):
            ### Exact values
            u_exact =          U_0*np.sin( x_data[i] )*np.cos( y_data[i] )*math.exp( ( -2.0 )*nu_ref*time )
            v_exact = ( -1.0 )*U_0*np.cos( x_data[i] )*np.sin( y_data[i] )*math.exp( ( -2.0 )*nu_ref*time )
            ### L2-norm errors
            u_L2_norm_error += ( u_exact - u_data[i] )**2.0
            u_L2_norm_den   += u_exact**2.0
            v_L2_norm_error += ( v_exact - v_data[i] )**2.0
            v_L2_norm_den   += v_exact**2.0
u_L2_norm_error = np.sqrt( u_L2_norm_error/u_L2_norm_den )
v_L2_norm_error = np.sqrt( v_L2_norm_error/v_L2_norm_den )


### Print results
print( 'L2-norm relative errors N = 16:' + ' u_L2_norm = ' + str( u_L2_norm_error ) + ', v_L2_norm = ' + str( v_L2_norm_error ) )

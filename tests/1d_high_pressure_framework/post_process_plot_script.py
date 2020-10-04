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
data_file = h5py.File( '1d_high_pressure_framework_0.h5', 'r' )
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

### Open NIST file
T_nist, P_nist, rho_nist, V_nist, e_nist, h_nist, s_nist, cv_nist, cp_nist, sos_nist, JT_nist, mu_nist, kappa_nist = np.loadtxt( 'nist_nitrogen_4MPa_0-1000K.csv', delimiter=',', unpack = 'True' )


### Plot density vs. temperature

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( T_nist, rho_nist, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{NIST}$' )
plt.plot( T_data, rho_data, linestyle = '--', color = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.tick_params(axis = 'both', top = True, bottom = True, right = 'True', left = 'True', direction = 'inout')
plt.xlim( 100, 200 )
plt.xticks( np.arange( 100, 210, 25 ) )
#plt.xscale( 'log' )
plt.xlabel( r'$T \thinspace \textrm{[K]}$' )
plt.ylim( 50, 800 )
plt.yticks( np.arange( 50, 810, 150 ) )
#plt.yscale( 'log' )
plt.ylabel( r'$\rho \thinspace \textrm{[kg/m}^\textrm{3}\textrm{]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'density_vs_temperature.eps', format = 'eps', bbox_inches = 'tight' )


### Plot internal energy vs. temperature

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( T_nist, e_nist, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{NIST}$' )
plt.plot( T_data, e_data/1.0e3 + ( 387.676 - 76.809 ), linestyle = '--', color = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.tick_params(axis = 'both', top = True, bottom = True, right = 'True', left = True, direction = 'inout')
plt.xlim( 100, 200 )
plt.xticks( np.arange( 100, 210, 25 ) )
#plt.xscale( 'log' )
plt.xlabel( r'$T \thinspace \textrm{[K]}$' )
plt.ylim( -100, 150 )
plt.yticks( np.arange( -100, 160, 50 ) )
#plt.yscale( 'log' )
plt.ylabel( r'$e \thinspace \textrm{[kJ/kg]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper left' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'internal_energy_vs_temperature.eps', format = 'eps', bbox_inches = 'tight' )


### Plot speed of sound vs. temperature

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( T_nist, sos_nist, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{NIST}$' )
plt.plot( T_data, sos_data, linestyle = '--', color = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.tick_params(axis = 'both', top = True, bottom = True, right = 'True', left = True, direction = 'inout')
plt.xlim( 100, 200 )
plt.xticks( np.arange( 100, 210, 25 ) )
#plt.xscale( 'log' )
plt.xlabel( r'$T \thinspace \textrm{[K]}$' )
plt.ylim( 100, 700 )
plt.yticks( np.arange( 100, 710, 100 ) )
#plt.yscale( 'log' )
plt.ylabel( r'$c \thinspace \textrm{[m/s]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'sound_speed_vs_temperature.eps', format = 'eps', bbox_inches = 'tight' )


### Plot viscosity vs. temperature

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( T_nist, mu_nist, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{NIST}$' )
plt.plot( T_data, mu_data, linestyle = '--', color = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.tick_params(axis = 'both', top = True, bottom = True, right = 'True', left = True, direction = 'inout')
plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
plt.xlim( 100, 200 )
plt.xticks( np.arange( 100, 210, 25 ) )
#plt.xscale( 'log' )
plt.xlabel( r'$T \thinspace \textrm{[K]}$' )
plt.ylim( 0.0, 0.00015 )
plt.yticks( np.arange( 0.0, 0.00016, 0.00005 ) )
#plt.yscale( 'log' )
plt.ylabel( r'$\mu \thinspace \textrm{[Pa}\cdot\textrm{s]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'viscosity_vs_temperature.eps', format = 'eps', bbox_inches = 'tight' )


### Plot thermal conductivity vs. temperature

# Clear plot
plt.clf()

# Read & Plot data
plt.scatter( T_nist, kappa_nist, s = 25, facecolors = 'none', edgecolors = 'black', label = r'$\textrm{NIST}$' )
plt.plot( T_data, kappa_data, linestyle = '--', color = 'black', label = r'$\textrm{RHEA}$' )

# Configure plot
plt.tick_params(axis = 'both', top = True, bottom = True, right = 'True', left = True, direction = 'inout')
plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
plt.xlim( 100, 200 )
plt.xticks( np.arange( 100, 210, 25 ) )
#plt.xscale( 'log' )
plt.xlabel( r'$T \thinspace \textrm{[K]}$' )
plt.ylim( 0.0, 0.15 )
plt.yticks( np.arange( 0.0, 0.16, 0.05 ) )
#plt.yscale( 'log' )
plt.ylabel( r'$\kappa \thinspace \textrm{[W/(m}\cdot\textrm{K)]}$' )
legend = plt.legend( shadow = False, fancybox = False, frameon = False, loc='upper right' )
plt.tick_params( axis = 'both', pad = 7.5 )
plt.savefig( 'thermal_conductivity_vs_temperature.eps', format = 'eps', bbox_inches = 'tight' )

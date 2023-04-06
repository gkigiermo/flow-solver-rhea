#!/usr/bin/python

import sys
import os
import numpy as np
import h5py 

   
### Open data file
data_file = h5py.File( '2D_differentially_heated_cavity_27010.h5', 'r' )
#print( list( data_file.keys() ) )
x_data       = data_file['x'][1,:,:]
y_data       = data_file['y'][1,:,:]
T_data       = data_file['T'][1,:,:]
num_points_x = x_data[0,:].size
num_points_y = y_data[:,0].size

### Reference parameters
L       = 1.0               # Cavity size [m]
T_cw    = 295.0             # Cold-wall temperature [K]
T_hw    = 305.0             # Hot-wall temperature [K]
Delta_T = T_hw - T_cw		# Temperature difference [K]


### Calculate Nusselt number at west boundary
Nu_avg   = 0.0          # Average Nusselt number
Nu_max   = -1.0         # Maximum Nusselt number
y_Nu_max = -1.0         # Position of maximum Nusselt number
Nu_min   = 1.0e6        # Minimum Nusselt number
y_Nu_min = -1.0         # Position of minimum Nusselt number
for j in range( 1, num_points_y - 1 ):
    i = 0		                                                                # index of boundary point
    dT_dx = ( T_data[j,i+1] - T_data[j,i] )/( x_data[j,i+1] - x_data[j,i] )     # Temperature gradient at the wall
    Nu = ( (-1.0)*L/Delta_T )*dT_dx                                             # Nusselt number
    Nu_avg += ( 1.0/( num_points_y - 2.0 ) )*Nu                                 # Update average Nusselt number
    if( Nu > Nu_max ):                                                          # Update maximum Nusselt number
        Nu_max   = Nu
        y_Nu_max = y_data[j,i]/L
    if( Nu < Nu_min ):                                                          # Update minimum Nusselt number
        Nu_min   = Nu
        y_Nu_min = y_data[j,i]/L

### Print results
print( "Nusselt numbers at west (hot) wall:" )
print( "Average Nu:", Nu_avg, "... Vahl Davis (1983) value is Nu_avg = 2.242" )
print( "Maximum Nu:", Nu_max, "at y/L =", y_Nu_max, "... Vahl Davis (1983) value is Nu_max = 3.545 at y/L = 0.149" )
print( "Minimum Nu:", Nu_min, "at y/L =", y_Nu_min, "... Vahl Davis (1983) value is Nu_min = 0.592 at y/L = 1.0" )

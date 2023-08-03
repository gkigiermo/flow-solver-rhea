#!/usr/bin/python

import sys
import os
import numpy as np
import h5py


### Set name of h5 files and direction & modification
original_name_data_file = 'restart_data_file.h5'
modified_name_data_file = 'modified_restart_data_file.h5'
direction_to_modify     = 'z-direction'     # x-direction (x), y-direction (y), z-direction (z)
modified_delta          = 0.1               # Width of modified plane
format_dataset          = np.float64

### Open original & modified data files
original_data_file = h5py.File( original_name_data_file, 'r' )
modified_data_file = h5py.File( modified_name_data_file, 'w' )

#### Copy attributes from original to extruded file
attribute_time = original_data_file.attrs['Time']
modified_data_file.attrs.create( 'Time', attribute_time ) 
attribute_iteration = original_data_file.attrs['Iteration']
modified_data_file.attrs.create( 'Iteration', attribute_iteration ) 
attribute_averaging_time = original_data_file.attrs['AveragingTime']
modified_data_file.attrs.create( 'AveragingTime', attribute_averaging_time ) 

#### Modify dataset
# Works with only 1 inner cell (3 total cells) ... cannot be modified!
datasets = list( original_data_file.keys() )
for ds in datasets:
    ### Create modified dataset with zeros
    array_with_values = np.zeros( [ num_points_z, num_points_y, num_points_x ], dtype = format_dataset )
    ### Copy data from original to modified file
    original_data = original_data_file[ds][:,:,:]
    for i in range( 0, num_points_x ):
        for j in range( 0, num_points_y ):
            for k in range( 0, num_points_z ):
                array_with_values[k,j,i] = original_data[k,j,i]
                if( ( direction_to_modify == 'x-direction' ) and ( ds == 'x' ) ):
                    array_with_values[k,j,i] = ( 1.0*i - 0.5 )*modified_delta 
                if( ( direction_to_modify == 'y-direction' ) and ( ds == 'y' ) ):
                    array_with_values[k,j,i] = ( 1.0*j - 0.5 )*modified_delta 
                if( ( direction_to_modify == 'z-direction' ) and ( ds == 'z' ) ):
                    array_with_values[k,j,i] = ( 1.0*k - 0.5 )*modified_delta 
    ### Add modified dataset to file
    modified_data_file.create_dataset( ds, data = array_with_values )

### Close original & extruded data files
original_data_file.close()
modified_data_file.close()

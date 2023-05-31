#!/usr/bin/python

import sys
import os
import numpy as np
import h5py


### Set name of h5 files and direction of extrusion
original_name_data_file = 'restart_data_file.h5'
extruded_name_data_file = 'extruded_restart_data_file.h5'
extrusion_direction     = 'z-direction'    # x-direction (x), y-direction (y), z-direction(z)
format_dataset          = np.float64

### Open original & extruded data files
original_data_file = h5py.File( original_name_data_file, 'r' )
extruded_data_file = h5py.File( extruded_name_data_file, 'w' )

### Print initial size of datasets
x_data       = original_data_file['x'][:,:,:]
num_points_x = x_data[0,0,:].size
num_points_y = x_data[0,:,0].size
num_points_z = x_data[:,0,0].size
print( '\nInitial size of datasets:' )
print( 'x-direction:', num_points_x, ', y-direction:', num_points_y, ', z-direction:', num_points_z )

#### Copy attributes from original to extruded file
attribute_time = original_data_file.attrs['Time']
extruded_data_file.attrs.create( 'Time', attribute_time ) 
attribute_iteration = original_data_file.attrs['Iteration']
extruded_data_file.attrs.create( 'Iteration', attribute_iteration ) 
attribute_averaging_time = original_data_file.attrs['AveragingTime']
extruded_data_file.attrs.create( 'AveragingTime', attribute_averaging_time ) 

#### Create extruded datasets
num_extrusion_planes = 1                # Number of additional extrusion planes ... cannot be modified!
datasets = list( original_data_file.keys() )
for ds in datasets:
    if( extrusion_direction == 'x-direction' ):
        ### Create extruded dataset with zeros
        array_with_values = np.zeros( [ num_points_z, num_points_y, num_points_x + num_extrusion_planes ], dtype = format_dataset )
        ### Copy data from original to extruded file
        original_data = original_data_file[ds][:,:,:]
        for i in range( 0, num_points_x + num_extrusion_planes ):
            for j in range( 0, num_points_y ):
                for k in range( 0, num_points_z ):
                    if( i < num_points_x ):
                        array_with_values[k,j,i] = original_data[k,j,i]
                    else:
                        if( ds == 'x' ):
                            array_with_values[k,j,i] = 2.0*original_data[k,j,i-1] - original_data[k,j,i-2] 
                        else:
                            array_with_values[k,j,i] = original_data[k,j,i-1]
        ### Add extruded dataset to file
        extruded_data_file.create_dataset( ds, data = array_with_values )
    if( extrusion_direction == 'y-direction' ):
        ### Create extruded dataset with zeros
        array_with_values = np.zeros( [ num_points_z, num_points_y + num_extrusion_planes, num_points_x ], dtype = format_dataset )
        ### Copy data from original to extruded file
        original_data = original_data_file[ds][:,:,:]
        for i in range( 0, num_points_x ):
            for j in range( 0, num_points_y + num_extrusion_planes ):
                for k in range( 0, num_points_z ):
                    if( j < num_points_y ):
                        array_with_values[k,j,i] = original_data[k,j,i]
                    else:
                        if( ds == 'y' ):
                            array_with_values[k,j,i] = 2.0*original_data[k,j-1,i] - original_data[k,j-2,i] 
                        else:
                            array_with_values[k,j,i] = original_data[k,j-1,i]
        ### Add extruded dataset to file
        extruded_data_file.create_dataset( ds, data = array_with_values )
    if( extrusion_direction == 'z-direction' ):
        ### Create extruded dataset with zeros
        array_with_values = np.zeros( [ num_points_z + num_extrusion_planes, num_points_y, num_points_x ], dtype = format_dataset )
        ### Copy data from original to extruded file
        original_data = original_data_file[ds][:,:,:]
        for i in range( 0, num_points_x ):
            for j in range( 0, num_points_y ):
                for k in range( 0, num_points_z + num_extrusion_planes ):
                    if( k < num_points_z ):
                        array_with_values[k,j,i] = original_data[k,j,i]
                    else:
                        if( ds == 'z' ):
                            array_with_values[k,j,i] = 2.0*original_data[k-1,j,i] - original_data[k-2,j,i]
                        else:
                            array_with_values[k,j,i] = original_data[k-1,j,i]
        ### Add extruded dataset to file
        extruded_data_file.create_dataset( ds, data = array_with_values )

### Print final size of datasets
x_data       = extruded_data_file['x'][:,:,:]
num_points_x = x_data[0,0,:].size
num_points_y = x_data[0,:,0].size
num_points_z = x_data[:,0,0].size
print( '\nFinal size of datasets:' )
print( 'x-direction:', num_points_x, ', y-direction:', num_points_y, ', z-direction:', num_points_z )

### Close original & extruded data files
original_data_file.close()
extruded_data_file.close()

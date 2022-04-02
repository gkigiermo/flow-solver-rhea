#!/usr/bin/python

import sys
import os
import numpy as np
import h5py


### Set name of h5 file and name & format of new dataset
name_data_file = 'restart_data_file.h5'
name_dataset   = 'new_dataset'
format_dataset = np.float64

### Open data file
data_file = h5py.File( name_data_file, 'a' )

### Print initial datasets and obtain size of data
print( '\nInitial list of datasets:' )
print( list( data_file.keys() ) )
x_data       = data_file['x'][:,:,:]
num_points_x = x_data[0,0,:].size
num_points_y = x_data[0,:,0].size
num_points_z = x_data[:,0,0].size

### Create new dataset with zeros
array_with_zeros = np.zeros( [ num_points_z, num_points_y, num_points_x ], dtype = format_dataset )

### Add new dataset to file
data_file.create_dataset( name_dataset, data = array_with_zeros )

### Print final datasets
print( '\nFinal list of datasets:' )
print( list( data_file.keys() ) )

### Close data file
data_file.close()

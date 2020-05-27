#!/usr/bin/python

################################################################################################################
#                                                                                                              #
# RHEA - the open-source Reproducible Hybrid-architecture flow solver Engineered for Academia                  #
#                                                                                                              #
# Rhea was the Titaness great Mother of the Gods, and goddess of female fertility, motherhood, and generation. #
# Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.           #
#                                                                                                              #
#                                                                                                              #
# REHA is distributed under the MIT License:                                                                   #
#                                                                                                              #
# Copyright 2020 Lluís Jofre Cruanyes & Guillermo Oyarzún Altamirano                                           #
#                                                                                                              #
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated #
# documentation files (the "Software"), to deal in the Software without restriction, including without         #
# limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of    #
# the Software, and to permit persons to whom the Software is furnished to do so, subject to the following     #
# conditions:                                                                                                  #
#                                                                                                              #
# The above copyright notice and this permission notice shall be included in all copies or substantial         #
# portions of the Software.                                                                                    #
#                                                                                                              #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT        #
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.          #
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      #
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE          #
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                                       #
#                                                                                                              #
################################################################################################################


# Primitive variables:
#   - Density rho
#   - Velocities u, v, w
#   - Specific total energy E
#     ... E = e + ke
#     ... is the sum of internal energy e
#     ... and kinetic energy ke = (u*u + v*v + w*w)/2

# Conserved variables:
#   - Mass rho
#   - Momentum rho*u, rho*v, rho*w
#   - Total energy rho*E

# Thermodynamic state:
#   - Pressure P
#   - Temperature T
#   - Speed of sound sos

# Thermophysical properties:
#   - Specific gas constant R_specific
#   - Ratio of heat capacities gamma
#   - Dynamic viscosity mu
#   - Thermal conductivity kappa



########## PHYTON MODULES ##########
import os, sys
import numpy as np


########## SET PARAMETERS ##########

### Fluid properties 
R_specific = 287.058     # Specific gas constant of air [J/(kg K)]
gamma      = 1.4         # Heat capacity ratio of air [-]
mu         = 0.0         # Dynamic viscosity of air [Pa s]
kappa      = 0.0         # Thermal conductivity of air [W/(m k)]

### Problem parameters
initial_rho_L = 1.0      # Initial density left [kg/m3]
initial_u_L   = 0.75     # Initial velocity left [m/s]
initial_P_L   = 1.0      # Initial pressure left [Pa]
initial_rho_R = 0.125    # Initial density right [kg/m3]
initial_u_R   = 0.0      # Initial velocity right [m/s]
initial_P_R   = 0.1      # Initial pressure right [Pa]
disc_x_0      = 0.3      # Position (x-direction) of discontinuity [m]
L_x           = 1.0      # Size of domain in x-direction
L_y           = 1.0e-2   # Size of domain in y-direction
L_z           = 1.0e-2   # Size of domain in z-direction
initial_time  = 0.0      # Initial time [s]
final_time    = 0.2      # Final time [s]
name_file_out = 'output_data.csv'   # Name of output data [-]

### Computational parameters
num_grid_x        = 100  # Number of internal grid points in the x-direction
num_grid_y        = 1    # Number of internal grid points in the y-direction
num_grid_z        = 1    # Number of internal grid points in the z-direction
CFL               = 0.9  # CFL coefficient
max_num_time_iter = 1e6  # Maximum number of time iterations

### Fixed parameters
num_sptl_dim = 3         # Number of spatial dimensions (fixed value)
rk_order     = 3         # Time-integration Runge-Kutta order (fixed value)
epsilon      = 1.0e-15   # Small epsilon number (fixed value)


########## ALLOCATE MEMORY ##########

### Mesh coordinates ... two positions added for boundary points
mesh = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2, num_sptl_dim ] )             # 3-D positions of mesh

### Primitive, conserved and thermodynamic variables ... two positions added for boundary points
rho_field  = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of rho
u_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of u
v_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of v
w_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of w
E_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of E
rhou_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of rhou
rhov_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of rhov
rhow_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of rhow
rhoE_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of rhoE
P_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of P
T_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of T
sos_field  = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                     # 3-D field of sos

### Time integration variables ... two positions added for boundary points
rho_0_field  = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D old field of rho
rhou_0_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D old field of rhou
rhov_0_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D old field of rhov
rhow_0_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D old field of rhow
rhoE_0_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D old field of rhoE

### Time integration fluxes ... two positions added for boundary points
rho_rk_fluxes  = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2, rk_order ] )       # 3-D Runge-Kutta fluxes of rho
rhou_rk_fluxes = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2, rk_order ] )       # 3-D Runge-Kutta fluxes of rhou
rhov_rk_fluxes = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2, rk_order ] )       # 3-D Runge-Kutta fluxes of rhov
rhow_rk_fluxes = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2, rk_order ] )       # 3-D Runge-Kutta fluxes of rhow
rhoE_rk_fluxes = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2, rk_order ] )       # 3-D Runge-Kutta fluxes of rhoE

### Inviscid fluxes ... two positions added for boundary points
rho_inv_flux  = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D inviscid fluxes of rho
rhou_inv_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D inviscid fluxes of rhou
rhov_inv_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D inviscid fluxes of rhov
rhow_inv_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D inviscid fluxes of rhow
rhoE_inv_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D inviscid fluxes of rhoE

### Viscous fluxes ... two positions added for boundary points
rhou_vis_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D viscous fluxes of rhou
rhov_vis_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D viscous fluxes of rhov
rhow_vis_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D viscous fluxes of rhow
rhoE_vis_flux = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                  # 3-D viscous fluxes of rhoE

### Body forces ... two positions added for boundary points
f_rhou_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D field of rhou
f_rhov_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D field of rhov
f_rhow_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D field of rhow
f_rhoE_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                   # 3-D field of rhoE



########## DEFINITIONS ##########

### Define centroids of spatial discretization
def spatial_discretization( grid ):

    # Cartesian uniform mesh
    delta_x = L_x/num_grid_x
    delta_y = L_y/num_grid_y
    delta_z = L_z/num_grid_z

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                grid[i][j][k][0] = ( i - 0.5 )*delta_x
                grid[i][j][k][1] = ( j - 0.5 )*delta_y
                grid[i][j][k][2] = ( k - 0.5 )*delta_z
    #print( grid )


### Initialize u, v, w, P and T variables
def initialize_uvwPT( u, v, w, P, T, grid ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                if( grid[i][j][k][0] < disc_x_0 ):
                    u[i][j][k] = initial_u_L
                    v[i][j][k] = 0.0
                    w[i][j][k] = 0.0
                    P[i][j][k] = initial_P_L
                    T[i][j][k] = ( 1.0/( initial_rho_L*R_specific ) )*P[i][j][k]
                else:
                    u[i][j][k] = initial_u_R
                    v[i][j][k] = 0.0
                    w[i][j][k] = 0.0
                    P[i][j][k] = initial_P_R
                    T[i][j][k] = ( 1.0/( initial_rho_R*R_specific ) )*P[i][j][k]
    #print( u )
    #print( v )
    #print( w )
    #print( P )
    #print( T )


### Initialize thermodynamic variables
def initialize_thermodynamics( rho, E, sos, u, v, w, P, T ):

    # Ideal-gas model:
    # rho = P/(R_specific*T) is density
    # e = P/(rho*(gamma - 1)) is specific internal energy
    # ke = (u*u + v*v + w*w)/2 is specific kinetic energy
    # E = e + ke is total energy
    # sos = sqrt(gamma*P/rho) is speed of sound

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                rho[i][j][k] = ( 1.0/( R_specific*T[i][j][k] ) )*P[i][j][k]        
                e            = ( 1.0/( rho[i][j][k]*( gamma - 1.0 ) ) )*P[i][j][k]
                ke           = 0.5*( u[i][j][k]**2.0 + v[i][j][k]**2.0 + w[i][j][k]**2.0 )
                E[i][j][k]   = e + ke
                sos[i][j][k] = np.sqrt( gamma*( ( 1.0/rho[i][j][k] )*P[i][j][k] ) )
    #print( rho )
    #print( E )
    #print( sos )


### Update conserved variables from primitive variables
def update_conserved( conserved, primitive, rho ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                conserved[i][j][k] = rho[i][j][k]*primitive[i][j][k]
    #print( conserved )


### Update field
def update_field( field_a, field_b ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                field_a[i][j][k] = field_b[i][j][k]
    #print( field_a )


### Calculate time step
def time_step( rho, u, v, w, sos, grid ):

    # Inviscid time step size for explicit schemes:
    # E. F. Toro.
    # Riemann solvers and numerical methods for fluid dynamics.
    # Springer, 2009.

    # Viscous time step size for explicit schemes:
    # E. Turkel, R.C. Swanson, V. N. Vatsa, J.A. White.
    # Multigrid for hypersonic viscous two- and three-dimensional flows.
    # NASA Contractor Report 187603, 1991.

    # Calculate specific heat capacities
    c_v = R_specific/( gamma - 1.0 )
    c_p = c_v*gamma

    # Calculate Prandtl (Pr) number
    Pr = gamma
    if( kappa > epsilon ):
        Pr = ( 1.0/( kappa ) )*c_p*mu

    # Initialize to largest float value
    delta_t = float( 'inf' )

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                ## Geometric stuff
                delta_x = 0.5*( grid[i+1][j][k][0] - grid[i-1][j][k][0] ) 
                delta_y = 0.5*( grid[i][j+1][k][1] - grid[i][j-1][k][1] ) 
                delta_z = 0.5*( grid[i][j][k+1][2] - grid[i][j][k-1][2] )                
                ## x-direction inviscid & viscous terms
                S_x     = abs( u[i][j][k] ) + sos[i][j][k]
                delta_t = min( delta_t, ( 1.0/S_x )*CFL*delta_x )
                delta_t = min( delta_t, ( 1.0/( mu*gamma + epsilon ) )*CFL*Pr*rho[i][j][k]*( delta_x**2.0 ) )
                ## y-direction inviscid & viscous terms
                S_y     = abs( v[i][j][k] ) + sos[i][j][k]
                delta_t = min( delta_t, ( 1.0/S_y )*CFL*delta_y )
                delta_t = min( delta_t, ( 1.0/( mu*gamma + epsilon ) )*CFL*Pr*rho[i][j][k]*( delta_y**2.0 ) )
                ## z-direction inviscid & viscous terms
                S_z     = abs( w[i][j][k] ) + sos[i][j][k]
                delta_t = min( delta_t, ( 1.0/S_z )*CFL*delta_z )
                delta_t = min( delta_t, ( 1.0/( mu*gamma + epsilon ) )*CFL*Pr*rho[i][j][k]*( delta_z**2.0 ) )
    #print( delta_t )

    # Return minimum time step
    return( delta_t )


### Calculate source terms
def source_terms( f_rhou, f_rhov, f_rhow, f_rhoE, rho, u, v, w, mesh ):

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                f_rhou[i][j][k] = 0.0
                f_rhov[i][j][k] = 0.0
                f_rhow[i][j][k] = 0.0
                f_rhoE[i][j][k] = 0.0
                #f_rhoE[i][j][k] = ( -1.0 )*( f_rhou[i][j][k]*u[i][j][k] + f_rhov[i][j][k]*v[i][j][k] + f_rhow[i][j][k]*w[i][j][k] )
    #print( f_rhou )
    #print( f_rhov )
    #print( f_rhow )
    #print( f_rhoE )


### calculate wave speeds
def waves_speed( rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R ):

    # HLLC approximate Riemann solver:
    # E. F. Toro.
    # Riemann solvers and numerical methods for fluid dynamics.
    # Springer, 2009.    

    rho_bar = 0.5*( rho_L + rho_R )
    a_bar   = 0.5*( a_L + a_R )
    P_pvrs  = 0.5*( P_L + P_R ) - 0.5*( u_R - u_L )*rho_bar*a_bar
    P_star  = max( 0.0, P_pvrs )
    q_L     = 1.0
    if( P_star > P_L ):
        q_L = np.sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_L ) - 1.0 ) )
    q_R     = 1.0
    if( P_star > P_R ):
        q_R = np.sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_R ) - 1.0 ) )                    
    S_L     = u_L - a_L*q_L 
    S_R     = u_R + a_R*q_R
    #print( S_L, S_R )

    # Return wave speed estimates
    return( S_L, S_R )


### Calculate HLLC flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE 
def HLLC_flux( F_L, F_R, U_L, U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type ):

    # HLLC approximate Riemann solver:
    # E. F. Toro.
    # Riemann solvers and numerical methods for fluid dynamics.
    # Springer, 2009.    

    S_L, S_R = waves_speed( rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R )
    S_star   = ( P_R - P_L + rho_L*u_L*( S_L - u_L ) - rho_R*u_R*( S_R - u_R ) )/( rho_L*( S_L - u_L ) - rho_R*( S_R - u_R ) )
    U_star_L = rho_L*( ( S_L - u_L )/( S_L - S_star ) )
    U_star_R = rho_R*( ( S_R - u_R )/( S_R - S_star ) )
    if( var_type == 0 ):
        U_star_L *= 1.0
        U_star_R *= 1.0        
    elif( var_type == 1 ):
        U_star_L *= S_star
        U_star_R *= S_star
    elif( var_type == 2 ):
        U_star_L *= v_L
        U_star_R *= v_R
    elif( var_type == 3 ):
        U_star_L *= w_L
        U_star_R *= w_R
    elif( var_type == 4 ):
        U_star_L *= ( E_L + ( S_star - u_L )*( S_star + P_L/( rho_L*( S_L - u_L ) ) ) )
        U_star_R *= ( E_R + ( S_star - u_R )*( S_star + P_R/( rho_R*( S_R - u_R ) ) ) )
    F_star_L = F_L + S_L*( U_star_L - U_L )
    F_star_R = F_R + S_R*( U_star_R - U_R )
    F = 0.0
    if( 0.0 <= S_L ):
        F = F_L
    elif( ( S_L <= 0.0 ) and ( 0.0 <= S_star ) ):
        F = F_star_L
    elif( ( S_star <= 0.0 ) and ( 0.0 <= S_R ) ):
        F = F_star_R
    elif( 0.0 >= S_R ):
        F = F_R
    #print( F )

    # Return F value
    return( F )


### Calculate inviscid fluxes
def inviscid_fluxes( rho_inv, rhou_inv, rhov_inv, rhow_inv, rhoE_inv, rho, u, v, w, E, P, sos, grid ):

    # First-order Godunov-type unsplit method for Euler equations:
    # E. F. Toro.
    # Riemann solvers and numerical methods for fluid dynamics.
    # Springer, 2009.

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                ## Geometric stuff
                delta_x = 0.5*( grid[i+1][j][k][0] - grid[i-1][j][k][0] )
                delta_y = 0.5*( grid[i][j+1][k][1] - grid[i][j-1][k][1] ) 
                delta_z = 0.5*( grid[i][j][k+1][2] - grid[i][j][k-1][2] )                
                ## x-direction i+1/2
                index_L = i;                  index_R = i + 1
                rho_L   = rho[index_L][j][k]; rho_R   = rho[index_R][j][k] 
                u_L     = u[index_L][j][k];   u_R     = u[index_R][j][k]
                v_L     = v[index_L][j][k];   v_R     = v[index_R][j][k]
                w_L     = w[index_L][j][k];   w_R     = w[index_R][j][k]
                E_L     = E[index_L][j][k];   E_R     = E[index_R][j][k]
                P_L     = P[index_L][j][k];   P_R     = P[index_R][j][k]
                a_L     = sos[index_L][j][k]; a_R     = sos[index_R][j][k]
                # rho
                var_type = 0
                rho_F_L  = rho[index_L][j][k]*u[index_L][j][k]
                rho_F_R  = rho[index_R][j][k]*u[index_R][j][k]
                rho_U_L  = rho[index_L][j][k]
                rho_U_R  = rho[index_R][j][k]
                rho_F_p  = HLLC_flux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 1
                rhou_F_L = rho[index_L][j][k]*u[index_L][j][k]*u[index_L][j][k] + P[index_L][j][k]
                rhou_F_R = rho[index_R][j][k]*u[index_R][j][k]*u[index_R][j][k] + P[index_R][j][k]
                rhou_U_L = rho[index_L][j][k]*u[index_L][j][k]
                rhou_U_R = rho[index_R][j][k]*u[index_R][j][k]
                rhou_F_p = HLLC_flux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_L = rho[index_L][j][k]*u[index_L][j][k]*v[index_L][j][k]
                rhov_F_R = rho[index_R][j][k]*u[index_R][j][k]*v[index_R][j][k]
                rhov_U_L = rho[index_L][j][k]*v[index_L][j][k]
                rhov_U_R = rho[index_R][j][k]*v[index_R][j][k]
                rhov_F_p = HLLC_flux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_L = rho[index_L][j][k]*u[index_L][j][k]*w[index_L][j][k]
                rhow_F_R = rho[index_R][j][k]*u[index_R][j][k]*w[index_R][j][k]
                rhow_U_L = rho[index_L][j][k]*w[index_L][j][k]
                rhow_U_R = rho[index_R][j][k]*w[index_R][j][k]
                rhow_F_p = HLLC_flux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_L = rho[index_L][j][k]*u[index_L][j][k]*E[index_L][j][k] + u[index_L][j][k]*P[index_L][j][k] 
                rhoE_F_R = rho[index_R][j][k]*u[index_R][j][k]*E[index_R][j][k] + u[index_R][j][k]*P[index_R][j][k] 
                rhoE_U_L = rho[index_L][j][k]*E[index_L][j][k]
                rhoE_U_R = rho[index_R][j][k]*E[index_R][j][k]
                rhoE_F_p = HLLC_flux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )                
                ## x-direction i-1/2
                index_L = i - 1;              index_R = i
                rho_L   = rho[index_L][j][k]; rho_R   = rho[index_R][j][k] 
                u_L     = u[index_L][j][k];   u_R     = u[index_R][j][k]
                v_L     = v[index_L][j][k];   v_R     = v[index_R][j][k]
                w_L     = w[index_L][j][k];   w_R     = w[index_R][j][k]
                E_L     = E[index_L][j][k];   E_R     = E[index_R][j][k]
                P_L     = P[index_L][j][k];   P_R     = P[index_R][j][k]
                a_L     = sos[index_L][j][k]; a_R     = sos[index_R][j][k]
                # rho
                var_type = 0
                rho_F_L  = rho[index_L][j][k]*u[index_L][j][k]
                rho_F_R  = rho[index_R][j][k]*u[index_R][j][k]
                rho_U_L  = rho[index_L][j][k]
                rho_U_R  = rho[index_R][j][k]
                rho_F_m  = HLLC_flux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 1
                rhou_F_L = rho[index_L][j][k]*u[index_L][j][k]*u[index_L][j][k] + P[index_L][j][k]
                rhou_F_R = rho[index_R][j][k]*u[index_R][j][k]*u[index_R][j][k] + P[index_R][j][k]
                rhou_U_L = rho[index_L][j][k]*u[index_L][j][k]
                rhou_U_R = rho[index_R][j][k]*u[index_R][j][k]
                rhou_F_m = HLLC_flux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_L = rho[index_L][j][k]*u[index_L][j][k]*v[index_L][j][k]
                rhov_F_R = rho[index_R][j][k]*u[index_R][j][k]*v[index_R][j][k]
                rhov_U_L = rho[index_L][j][k]*v[index_L][j][k]
                rhov_U_R = rho[index_R][j][k]*v[index_R][j][k]
                rhov_F_m = HLLC_flux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_L = rho[index_L][j][k]*u[index_L][j][k]*w[index_L][j][k]
                rhow_F_R = rho[index_R][j][k]*u[index_R][j][k]*w[index_R][j][k]
                rhow_U_L = rho[index_L][j][k]*w[index_L][j][k]
                rhow_U_R = rho[index_R][j][k]*w[index_R][j][k]
                rhow_F_m = HLLC_flux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_L = rho[index_L][j][k]*u[index_L][j][k]*E[index_L][j][k] + u[index_L][j][k]*P[index_L][j][k] 
                rhoE_F_R = rho[index_R][j][k]*u[index_R][j][k]*E[index_R][j][k] + u[index_R][j][k]*P[index_R][j][k]
                rhoE_U_L = rho[index_L][j][k]*E[index_L][j][k]
                rhoE_U_R = rho[index_R][j][k]*E[index_R][j][k]
                rhoE_F_m = HLLC_flux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                ## Fluxes x-direction
                rho_inv[i][j][k]  = ( 1.0/delta_x )*( rho_F_p - rho_F_m )
                rhou_inv[i][j][k] = ( 1.0/delta_x )*( rhou_F_p - rhou_F_m )
                rhov_inv[i][j][k] = ( 1.0/delta_x )*( rhov_F_p - rhov_F_m )
                rhow_inv[i][j][k] = ( 1.0/delta_x )*( rhow_F_p - rhow_F_m )
                rhoE_inv[i][j][k] = ( 1.0/delta_x )*( rhoE_F_p - rhoE_F_m )
                ## y-direction j+1/2
                index_L = j;                  index_R = j + 1
                rho_L   = rho[i][index_L][k]; rho_R   = rho[i][index_R][k] 
                u_L     = v[i][index_L][k];   u_R     = v[i][index_R][k]
                v_L     = u[i][index_L][k];   v_R     = u[i][index_R][k]
                w_L     = w[i][index_L][k];   w_R     = w[i][index_R][k]
                E_L     = E[i][index_L][k];   E_R     = E[i][index_R][k]
                P_L     = P[i][index_L][k];   P_R     = P[i][index_R][k]
                a_L     = sos[i][index_L][k]; a_R     = sos[i][index_R][k]
                # rho
                var_type = 0
                rho_F_L  = rho[i][index_L][k]*v[i][index_L][k]
                rho_F_R  = rho[i][index_R][k]*v[i][index_R][k]
                rho_U_L  = rho[i][index_L][k]
                rho_U_R  = rho[i][index_R][k]
                rho_F_p  = HLLC_flux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 2
                rhou_F_L = rho[i][index_L][k]*v[i][index_L][k]*u[i][index_L][k]
                rhou_F_R = rho[i][index_R][k]*v[i][index_R][k]*u[i][index_R][k] 
                rhou_U_L = rho[i][index_L][k]*u[i][index_L][k]
                rhou_U_R = rho[i][index_R][k]*u[i][index_R][k]
                rhou_F_p = HLLC_flux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 1
                rhov_F_L = rho[i][index_L][k]*v[i][index_L][k]*v[i][index_L][k] + P[i][index_L][k]
                rhov_F_R = rho[i][index_R][k]*v[i][index_R][k]*v[i][index_R][k] + P[i][index_R][k]
                rhov_U_L = rho[i][index_L][k]*v[i][index_L][k]
                rhov_U_R = rho[i][index_R][k]*v[i][index_R][k]
                rhov_F_p = HLLC_flux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_L = rho[i][index_L][k]*v[i][index_L][k]*w[i][index_L][k]
                rhow_F_R = rho[i][index_R][k]*v[i][index_R][k]*w[i][index_R][k]
                rhow_U_L = rho[i][index_L][k]*w[i][index_L][k]
                rhow_U_R = rho[i][index_R][k]*w[i][index_R][k]
                rhow_F_p = HLLC_flux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_L = rho[i][index_L][k]*v[i][index_L][k]*E[i][index_L][k] + v[i][index_L][k]*P[i][index_L][k] 
                rhoE_F_R = rho[i][index_R][k]*v[i][index_R][k]*E[i][index_R][k] + v[i][index_R][k]*P[i][index_R][k] 
                rhoE_U_L = rho[i][index_L][k]*E[i][index_L][k]
                rhoE_U_R = rho[i][index_R][k]*E[i][index_R][k]
                rhoE_F_p = HLLC_flux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )                
                ## y-direction j-1/2
                index_L = j - 1;              index_R = j
                rho_L   = rho[i][index_L][k]; rho_R   = rho[i][index_R][k] 
                u_L     = v[i][index_L][k];   u_R     = v[i][index_R][k]
                v_L     = u[i][index_L][k];   v_R     = u[i][index_R][k]
                w_L     = w[i][index_L][k];   w_R     = w[i][index_R][k]
                E_L     = E[i][index_L][k];   E_R     = E[i][index_R][k]
                P_L     = P[i][index_L][k];   P_R     = P[i][index_R][k]
                a_L     = sos[i][index_L][k]; a_R     = sos[i][index_R][k]
                # rho
                var_type = 0
                rho_F_L  = rho[i][index_L][k]*v[i][index_L][k]
                rho_F_R  = rho[i][index_R][k]*v[i][index_R][k]
                rho_U_L  = rho[i][index_L][k]
                rho_U_R  = rho[i][index_R][k]
                rho_F_m  = HLLC_flux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 2
                rhou_F_L = rho[i][index_L][k]*v[i][index_L][k]*u[i][index_L][k]
                rhou_F_R = rho[i][index_R][k]*v[i][index_R][k]*u[i][index_R][k] 
                rhou_U_L = rho[i][index_L][k]*u[i][index_L][k]
                rhou_U_R = rho[i][index_R][k]*u[i][index_R][k]
                rhou_F_m = HLLC_flux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 1
                rhov_F_L = rho[i][index_L][k]*v[i][index_L][k]*v[i][index_L][k] + P[i][index_L][k]
                rhov_F_R = rho[i][index_R][k]*v[i][index_R][k]*v[i][index_R][k] + P[i][index_R][k]
                rhov_U_L = rho[i][index_L][k]*v[i][index_L][k]
                rhov_U_R = rho[i][index_R][k]*v[i][index_R][k]
                rhov_F_m = HLLC_flux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_L = rho[i][index_L][k]*v[i][index_L][k]*w[i][index_L][k]
                rhow_F_R = rho[i][index_R][k]*v[i][index_R][k]*w[i][index_R][k]
                rhow_U_L = rho[i][index_L][k]*w[i][index_L][k]
                rhow_U_R = rho[i][index_R][k]*w[i][index_R][k]
                rhow_F_m = HLLC_flux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_L = rho[i][index_L][k]*v[i][index_L][k]*E[i][index_L][k] + v[i][index_L][k]*P[i][index_L][k] 
                rhoE_F_R = rho[i][index_R][k]*v[i][index_R][k]*E[i][index_R][k] + v[i][index_R][k]*P[i][index_R][k]
                rhoE_U_L = rho[i][index_L][k]*E[i][index_L][k]
                rhoE_U_R = rho[i][index_R][k]*E[i][index_R][k]
                rhoE_F_m = HLLC_flux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                ## Fluxes y-direction
                rho_inv[i][j][k]  += ( 1.0/delta_y )*( rho_F_p - rho_F_m )
                rhou_inv[i][j][k] += ( 1.0/delta_y )*( rhou_F_p - rhou_F_m )
                rhov_inv[i][j][k] += ( 1.0/delta_y )*( rhov_F_p - rhov_F_m )
                rhow_inv[i][j][k] += ( 1.0/delta_y )*( rhow_F_p - rhow_F_m )
                rhoE_inv[i][j][k] += ( 1.0/delta_y )*( rhoE_F_p - rhoE_F_m )
                ## z-direction k+1/2
                index_L = k;                  index_R = k + 1
                rho_L   = rho[i][j][index_L]; rho_R   = rho[i][j][index_R] 
                u_L     = w[i][j][index_L];   u_R     = w[i][j][index_R]
                v_L     = v[i][j][index_L];   v_R     = v[i][j][index_R]
                w_L     = u[i][j][index_L];   w_R     = u[i][j][index_R]
                E_L     = E[i][j][index_L];   E_R     = E[i][j][index_R]
                P_L     = P[i][j][index_L];   P_R     = P[i][j][index_R]
                a_L     = sos[i][j][index_L]; a_R     = sos[i][j][index_R]
                # rho
                var_type = 0
                rho_F_L  = rho[i][j][index_L]*w[i][j][index_L]
                rho_F_R  = rho[i][j][index_R]*w[i][j][index_R]
                rho_U_L  = rho[i][j][index_L]
                rho_U_R  = rho[i][j][index_R]
                rho_F_p  = HLLC_flux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 3
                rhou_F_L = rho[i][j][index_L]*w[i][j][index_L]*u[i][j][index_L]
                rhou_F_R = rho[i][j][index_R]*w[i][j][index_R]*u[i][j][index_R]
                rhou_U_L = rho[i][j][index_L]*u[i][j][index_L]
                rhou_U_R = rho[i][j][index_R]*u[i][j][index_R]
                rhou_F_p = HLLC_flux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_L = rho[i][j][index_L]*w[i][j][index_L]*v[i][j][index_L]
                rhov_F_R = rho[i][j][index_R]*w[i][j][index_R]*v[i][j][index_R]
                rhov_U_L = rho[i][j][index_L]*v[i][j][index_L]
                rhov_U_R = rho[i][j][index_R]*v[i][j][index_R]
                rhov_F_p = HLLC_flux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 1
                rhow_F_L = rho[i][j][index_L]*w[i][j][index_L]*w[i][j][index_L] + P[i][j][index_L]
                rhow_F_R = rho[i][j][index_R]*w[i][j][index_R]*w[i][j][index_R] + P[i][j][index_R]
                rhow_U_L = rho[i][j][index_L]*w[i][j][index_L]
                rhow_U_R = rho[i][j][index_R]*w[i][j][index_R]
                rhow_F_p = HLLC_flux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_L = rho[i][j][index_L]*w[i][j][index_L]*E[i][j][index_L] + w[i][j][index_L]*P[i][j][index_L] 
                rhoE_F_R = rho[i][j][index_R]*w[i][j][index_R]*E[i][j][index_R] + w[i][j][index_R]*P[i][j][index_R] 
                rhoE_U_L = rho[i][j][index_L]*E[i][j][index_L]
                rhoE_U_R = rho[i][j][index_R]*E[i][j][index_R]
                rhoE_F_p = HLLC_flux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )                
                ## z-direction k-1/2
                index_L = k - 1;              index_R = k
                rho_L   = rho[i][j][index_L]; rho_R   = rho[i][j][index_R] 
                u_L     = w[i][j][index_L];   u_R     = w[i][j][index_R]
                v_L     = v[i][j][index_L];   v_R     = v[i][j][index_R]
                w_L     = u[i][j][index_L];   w_R     = u[i][j][index_R]
                E_L     = E[i][j][index_L];   E_R     = E[i][j][index_R]
                P_L     = P[i][j][index_L];   P_R     = P[i][j][index_R]
                a_L     = sos[i][j][index_L]; a_R     = sos[i][j][index_R]
                # rho
                var_type = 0
                rho_F_L  = rho[i][j][index_L]*w[i][j][index_L]
                rho_F_R  = rho[i][j][index_R]*w[i][j][index_R]
                rho_U_L  = rho[i][j][index_L]
                rho_U_R  = rho[i][j][index_R]
                rho_F_m  = HLLC_flux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 3
                rhou_F_L = rho[i][j][index_L]*w[i][j][index_L]*u[i][j][index_L]
                rhou_F_R = rho[i][j][index_R]*w[i][j][index_R]*u[i][j][index_R]
                rhou_U_L = rho[i][j][index_L]*u[i][j][index_L]
                rhou_U_R = rho[i][j][index_R]*u[i][j][index_R]
                rhou_F_m = HLLC_flux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_L = rho[i][j][index_L]*w[i][j][index_L]*v[i][j][index_L]
                rhov_F_R = rho[i][j][index_R]*w[i][j][index_R]*v[i][j][index_R]
                rhov_U_L = rho[i][j][index_L]*v[i][j][index_L]
                rhov_U_R = rho[i][j][index_R]*v[i][j][index_R]
                rhov_F_m = HLLC_flux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 1
                rhow_F_L = rho[i][j][index_L]*w[i][j][index_L]*w[i][j][index_L] + P[i][j][index_L]
                rhow_F_R = rho[i][j][index_R]*w[i][j][index_R]*w[i][j][index_R] + P[i][j][index_R]
                rhow_U_L = rho[i][j][index_L]*w[i][j][index_L]
                rhow_U_R = rho[i][j][index_R]*w[i][j][index_R]
                rhow_F_m = HLLC_flux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_L = rho[i][j][index_L]*w[i][j][index_L]*E[i][j][index_L] + w[i][j][index_L]*P[i][j][index_L] 
                rhoE_F_R = rho[i][j][index_R]*w[i][j][index_R]*E[i][j][index_R] + w[i][j][index_R]*P[i][j][index_R]
                rhoE_U_L = rho[i][j][index_L]*E[i][j][index_L]
                rhoE_U_R = rho[i][j][index_R]*E[i][j][index_R]
                rhoE_F_m = HLLC_flux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                ## Fluxes z-direction
                rho_inv[i][j][k]  += ( 1.0/delta_z )*( rho_F_p - rho_F_m )
                rhou_inv[i][j][k] += ( 1.0/delta_z )*( rhou_F_p - rhou_F_m )
                rhov_inv[i][j][k] += ( 1.0/delta_z )*( rhov_F_p - rhov_F_m )
                rhow_inv[i][j][k] += ( 1.0/delta_z )*( rhow_F_p - rhow_F_m )
                rhoE_inv[i][j][k] += ( 1.0/delta_z )*( rhoE_F_p - rhoE_F_m )
    #print( rho_inv )
    #print( rhou_inv )
    #print( rhov_inv )
    #print( rhow_inv )
    #print( rhoE_inv )


### Calculate viscous fluxes
def viscous_fluxes( rhou_vis, rhov_vis, rhow_vis, rhoE_vis, u, v, w, T, f_rhou, f_rhov, f_rhow, f_rhoE, grid ):

    # Second-order central finite differences for derivatives:
    # P. Moin.
    # Fundamentals of engineering numerical analysis.
    # Cambridge University Press, 2010.

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                # Geometric stuff
                delta_x = grid[i+1][j][k][0] - grid[i-1][j][k][0] 
                delta_y = grid[i][j+1][k][1] - grid[i][j-1][k][1] 
                delta_z = grid[i][j][k+1][2] - grid[i][j][k-1][2] 
                # Divergence of tau tensor terms
                div_tau_xx  = ( 2.0/delta_x )*( 2.0*mu*( ( ( u[i+1][j][k] - u[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                       - ( ( u[i][j][k] - u[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                div_tau_xx -= ( 2.0/delta_x )*( ( 2.0/3.0 )*mu*( ( ( u[i+1][j][k] - u[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                               - ( ( u[i][j][k] - u[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                div_tau_xx -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( ( ( v[i+1][j+1][k] - v[i+1][j-1][k] )/delta_y ) 
                                                               - ( ( v[i-1][j+1][k] - v[i-1][j-1][k] )/delta_y ) ) )                
                div_tau_xx -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( ( ( w[i+1][j][k+1] - w[i+1][j][k-1] )/delta_z ) 
                                                               - ( ( w[i-1][j][k+1] - w[i-1][j][k-1] )/delta_z ) ) )                
                div_tau_xy  = ( 2.0/delta_x )*( mu*( ( ( v[i+1][j][k] - v[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                   - ( ( v[i][j][k] - v[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                div_tau_xy += ( 1.0/delta_x )*( mu*( ( ( u[i+1][j+1][k] - u[i+1][j-1][k] )/delta_y ) 
                                                   - ( ( u[i-1][j+1][k] - u[i-1][j-1][k] )/delta_y ) ) )
                div_tau_xz  = ( 2.0/delta_x )*( mu*( ( ( w[i+1][j][k] - w[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                   - ( ( w[i][j][k] - w[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                div_tau_xz += ( 1.0/delta_x )*( mu*( ( ( u[i+1][j][k+1] - u[i+1][j][k-1] )/delta_z ) 
                                                   - ( ( u[i-1][j][k+1] - u[i-1][j][k-1] )/delta_z ) ) )
                div_tau_yx  = ( 2.0/delta_y )*( mu*( ( ( u[i][j+1][k] - u[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                   - ( ( u[i][j][k] - u[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                div_tau_yx += ( 1.0/delta_y )*( mu*( ( ( v[i+1][j+1][k] - v[i-1][j+1][k] )/delta_x ) 
                                                   - ( ( v[i+1][j-1][k] - v[i-1][j-1][k] )/delta_x ) ) )
                div_tau_yy  = ( 2.0/delta_y )*( 2.0*mu*( ( ( v[i][j+1][k] - v[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                       - ( ( v[i][j][k] - v[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                div_tau_yy -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( ( ( u[i+1][j+1][k] - u[i-1][j+1][k] )/delta_x ) 
                                                               - ( ( u[i+1][j-1][k] - u[i-1][j-1][k] )/delta_x ) ) )
                div_tau_yy -= ( 2.0/delta_y )*( ( 2.0/3.0 )*mu*( ( ( v[i][j+1][k] - v[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                               - ( ( v[i][j][k] - v[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                div_tau_yy -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( ( ( w[i][j+1][k+1] - w[i][j+1][k-1] )/delta_z ) 
                                                               - ( ( w[i][j-1][k+1] - w[i][j-1][k-1] )/delta_z ) ) )                
                div_tau_yz  = ( 2.0/delta_y )*( mu*( ( ( w[i][j+1][k] - w[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                   - ( ( w[i][j][k] - w[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                div_tau_yz += ( 1.0/delta_y )*( mu*( ( ( v[i][j+1][k+1] - v[i][j+1][k-1] )/delta_z ) 
                                                   - ( ( v[i][j-1][k+1] - v[i][j-1][k-1] )/delta_z ) ) )                
                div_tau_zx  = ( 2.0/delta_z )*( mu*( ( ( u[i][j][k+1] - u[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                   - ( ( u[i][j][k] - u[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                div_tau_zx += ( 1.0/delta_z )*( mu*( ( ( w[i+1][j][k+1] - w[i-1][j][k+1] )/delta_x ) 
                                                   - ( ( w[i+1][j][k-1] - w[i-1][j][k-1] )/delta_x ) ) )
                div_tau_zy  = ( 2.0/delta_z )*( mu*( ( ( v[i][j][k+1] - v[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                   - ( ( v[i][j][k] - v[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                div_tau_zy += ( 1.0/delta_z )*( mu*( ( ( w[i][j+1][k+1] - w[i][j-1][k+1] )/delta_y ) 
                                                   - ( ( w[i][j+1][k-1] - w[i][j-1][k-1] )/delta_y ) ) )                
                div_tau_zz  = ( 2.0/delta_z )*( 2.0*mu*( ( ( w[i][j][k+1] - w[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                       - ( ( w[i][j][k] - w[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                div_tau_zz -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( ( ( u[i+1][j][k+1] - u[i-1][j][k+1] )/delta_x ) 
                                                               - ( ( u[i+1][j][k-1] - u[i-1][j][k-1] )/delta_x ) ) )
                div_tau_zz -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( ( ( v[i][j+1][k+1] - v[i][j-1][k+1] )/delta_y ) 
                                                               - ( ( v[i][j+1][k-1] - v[i][j-1][k-1] )/delta_y ) ) )                
                div_tau_zz -= ( 2.0/delta_z )*( ( 2.0/3.0 )*mu*( ( ( w[i][j][k+1] - w[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                               - ( ( w[i][j][k] - w[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                # Fourier term
                div_q = ( -1.0 )*( ( 2.0/delta_x )*( kappa*( ( ( T[i+1][j][k] - T[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                           - ( ( T[i][j][k] - T[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) ) 
                                 + ( 2.0/delta_y )*( kappa*( ( ( T[i][j+1][k] - T[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                           - ( ( T[i][j][k] - T[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                                 + ( 2.0/delta_z )*( kappa*( ( ( T[i][j][k+1] - T[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                           - ( ( T[i][j][k] - T[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) ) )
                # Divergence of tau*velocity terms
                vel_p = 0.5*( u[i+1][j][k] + u[i][j][k] ); vel_m = 0.5*( u[i][j][k] + u[i-1][j][k] )
                div_tauuvw_x  = ( 2.0/delta_x )*( 2.0*mu*( vel_p*( ( u[i+1][j][k] - u[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                         - vel_m*( ( u[i][j][k] - u[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                vel_p = 0.5*( u[i+1][j][k] + u[i][j][k] ); vel_m = 0.5*( u[i][j][k] + u[i-1][j][k] )
                div_tauuvw_x -= ( 2.0/delta_x )*( ( 2.0/3.0 )*mu*( vel_p*( ( u[i+1][j][k] - u[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                                 - vel_m*( ( u[i][j][k] - u[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                vel_p = 0.5*( u[i+1][j][k] + u[i+1][j][k] ); vel_m = 0.5*( u[i-1][j][k] + u[i-1][j][k] )
                div_tauuvw_x -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( vel_p*( ( v[i+1][j+1][k] - v[i+1][j-1][k] )/delta_y ) 
                                                                 - vel_m*( ( v[i-1][j+1][k] - v[i-1][j-1][k] )/delta_y ) ) )                
                vel_p = 0.5*( u[i+1][j][k] + u[i+1][j][k] ); vel_m = 0.5*( u[i-1][j][k] + u[i-1][j][k] )
                div_tauuvw_x -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( vel_p*( ( w[i+1][j][k+1] - w[i+1][j][k-1] )/delta_z ) 
                                                                 - vel_m*( ( w[i-1][j][k+1] - w[i-1][j][k-1] )/delta_z ) ) )                
                vel_p = 0.5*( v[i+1][j][k] + v[i][j][k] ); vel_m = 0.5*( v[i][j][k] + v[i-1][j][k] )
                div_tauuvw_x += ( 2.0/delta_x )*( mu*( vel_p*( ( v[i+1][j][k] - v[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                     - vel_m*( ( v[i][j][k] - v[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                vel_p = 0.5*( v[i+1][j][k] + v[i+1][j][k] ); vel_m = 0.5*( v[i-1][j][k] + v[i-1][j][k] )
                div_tauuvw_x += ( 1.0/delta_x )*( mu*( vel_p*( ( u[i+1][j+1][k] - u[i+1][j-1][k] )/delta_y ) 
                                                     - vel_m*( ( u[i-1][j+1][k] - u[i-1][j-1][k] )/delta_y ) ) )
                vel_p = 0.5*( w[i+1][j][k] + w[i][j][k] ); vel_m = 0.5*( w[i][j][k] + w[i-1][j][k] )
                div_tauuvw_x += ( 2.0/delta_x )*( mu*( vel_p*( ( w[i+1][j][k] - w[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] ) ) 
                                                     - vel_m*( ( w[i][j][k] - w[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) ) ) )
                vel_p = 0.5*( w[i+1][j][k] + w[i+1][j][k] ); vel_m = 0.5*( w[i-1][j][k] + w[i-1][j][k] )
                div_tauuvw_x += ( 1.0/delta_x )*( mu*( vel_p*( ( u[i+1][j][k+1] - u[i+1][j][k-1] )/delta_z ) 
                                                     - vel_m*( ( u[i-1][j][k+1] - u[i-1][j][k-1] )/delta_z ) ) )
                vel_p = 0.5*( u[i][j+1][k] + u[i][j][k] ); vel_m = 0.5*( u[i][j][k] + u[i][j-1][k] )
                div_tauuvw_y  = ( 2.0/delta_y )*( mu*( vel_p*( ( u[i][j+1][k] - u[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                     - vel_m*( ( u[i][j][k] - u[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                vel_p = 0.5*( u[i][j+1][k] + u[i][j+1][k] ); vel_m = 0.5*( u[i][j-1][k] + u[i][j-1][k] )
                div_tauuvw_y += ( 1.0/delta_y )*( mu*( vel_p*( ( v[i+1][j+1][k] - v[i-1][j+1][k] )/delta_x ) 
                                                     - vel_m*( ( v[i+1][j-1][k] - v[i-1][j-1][k] )/delta_x ) ) )
                vel_p = 0.5*( v[i][j+1][k] + v[i][j][k] ); vel_m = 0.5*( v[i][j][k] + v[i][j-1][k] )
                div_tauuvw_y += ( 2.0/delta_y )*( 2.0*mu*( vel_p*( ( v[i][j+1][k] - v[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                         - vel_m*( ( v[i][j][k] - v[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                vel_p = 0.5*( v[i][j+1][k] + v[i][j+1][k] ); vel_m = 0.5*( v[i][j-1][k] + v[i][j-1][k] )
                div_tauuvw_y -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( vel_p*( ( u[i+1][j+1][k] - u[i-1][j+1][k] )/delta_x ) 
                                                                 - vel_m*( ( u[i+1][j-1][k] - u[i-1][j-1][k] )/delta_x ) ) )
                vel_p = 0.5*( v[i][j+1][k] + v[i][j][k] ); vel_m = 0.5*( v[i][j][k] + v[i][j-1][k] )
                div_tauuvw_y -= ( 2.0/delta_y )*( ( 2.0/3.0 )*mu*( vel_p*( ( v[i][j+1][k] - v[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                                 - vel_m*( ( v[i][j][k] - v[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                vel_p = 0.5*( v[i][j+1][k] + v[i][j+1][k] ); vel_m = 0.5*( v[i][j-1][k] + v[i][j-1][k] )
                div_tauuvw_y -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( vel_p*( ( w[i][j+1][k+1] - w[i][j+1][k-1] )/delta_z ) 
                                                                 - vel_m*( ( w[i][j-1][k+1] - w[i][j-1][k-1] )/delta_z ) ) )                
                vel_p = 0.5*( w[i][j+1][k] + w[i][j][k] ); vel_m = 0.5*( w[i][j][k] + w[i][j-1][k] )
                div_tauuvw_y += ( 2.0/delta_y )*( mu*( vel_p*( ( w[i][j+1][k] - w[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] ) ) 
                                                     - vel_m*( ( w[i][j][k] - w[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) ) ) )
                vel_p = 0.5*( w[i][j+1][k] + w[i][j+1][k] ); vel_m = 0.5*( w[i][j-1][k] + w[i][j-1][k] )
                div_tauuvw_y += ( 1.0/delta_y )*( mu*( vel_p*( ( v[i][j+1][k+1] - v[i][j+1][k-1] )/delta_z ) 
                                                     - vel_m*( ( v[i][j-1][k+1] - v[i][j-1][k-1] )/delta_z ) ) )                
                vel_p = 0.5*( u[i][j][k+1] + u[i][j][k] ); vel_m = 0.5*( u[i][j][k] + u[i][j][k-1] )
                div_tauuvw_z  = ( 2.0/delta_z )*( mu*( vel_p*( ( u[i][j][k+1] - u[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                     - vel_m*( ( u[i][j][k] - u[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                vel_p = 0.5*( u[i][j][k+1] + u[i][j][k+1] ); vel_m = 0.5*( u[i][j][k-1] + u[i][j][k-1] )
                div_tauuvw_z += ( 1.0/delta_z )*( mu*( vel_p*( ( w[i+1][j][k+1] - w[i-1][j][k+1] )/delta_x ) 
                                                     - vel_m*( ( w[i+1][j][k-1] - w[i-1][j][k-1] )/delta_x ) ) )
                vel_p = 0.5*( v[i][j][k+1] + v[i][j][k] ); vel_m = 0.5*( v[i][j][k] + v[i][j][k-1] )
                div_tauuvw_z += ( 2.0/delta_z )*( mu*( vel_p*( ( v[i][j][k+1] - v[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                     - vel_m*( ( v[i][j][k] - v[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                vel_p = 0.5*( v[i][j][k+1] + v[i][j][k+1] ); vel_m = 0.5*( v[i][j][k-1] + v[i][j][k-1] )
                div_tauuvw_z += ( 1.0/delta_z )*( mu*( vel_p*( ( w[i][j+1][k+1] - w[i][j-1][k+1] )/delta_y ) 
                                                     - vel_m*( ( w[i][j+1][k-1] - w[i][j-1][k-1] )/delta_y ) ) )                
                vel_p = 0.5*( w[i][j][k+1] + w[i][j][k] ); vel_m = 0.5*( w[i][j][k] + w[i][j][k-1] )
                div_tauuvw_z += ( 2.0/delta_z )*( 2.0*mu*( vel_p*( ( w[i][j][k+1] - w[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                         - vel_m*( ( w[i][j][k] - w[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                vel_p = 0.5*( w[i][j][k+1] + w[i][j][k+1] ); vel_m = 0.5*( w[i][j][k-1] + w[i][j][k-1] )
                div_tauuvw_z -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( vel_p*( ( u[i+1][j][k+1] - u[i-1][j][k+1] )/delta_x ) 
                                                                 - vel_m*( ( u[i+1][j][k-1] - u[i-1][j][k-1] )/delta_x ) ) )
                vel_p = 0.5*( w[i][j][k+1] + w[i][j][k+1] ); vel_m = 0.5*( w[i][j][k-1] + w[i][j][k-1] )
                div_tauuvw_z -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( vel_p*( ( v[i][j+1][k+1] - v[i][j-1][k+1] )/delta_y ) 
                                                                 - vel_m*( ( v[i][j+1][k-1] - v[i][j-1][k-1] )/delta_y ) ) )                
                vel_p = 0.5*( w[i][j][k+1] + w[i][j][k] ); vel_m = 0.5*( w[i][j][k] + w[i][j][k-1] )
                div_tauuvw_z -= ( 2.0/delta_z )*( ( 2.0/3.0 )*mu*( vel_p*( ( w[i][j][k+1] - w[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] ) ) 
                                                                 - vel_m*( ( w[i][j][k] - w[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) )
                # Work of viscous stresses
                div_tauuvw = div_tauuvw_x + div_tauuvw_y + div_tauuvw_z 
                # Work of sources
                f_rhouvw = f_rhou[i][j][k]*u[i][j][k] + f_rhov[i][j][k]*v[i][j][k] + f_rhow[i][j][k]*w[i][j][k] 
                # Viscous fluxes
                rhou_vis[i][j][k] = div_tau_xx + div_tau_yx + div_tau_zx + f_rhou[i][j][k] 
                rhov_vis[i][j][k] = div_tau_xy + div_tau_yy + div_tau_zy + f_rhov[i][j][k]
                rhow_vis[i][j][k] = div_tau_xz + div_tau_yz + div_tau_zz + f_rhow[i][j][k]
                rhoE_vis[i][j][k] = ( -1.0 )*div_q + div_tauuvw + f_rhouvw + f_rhoE[i][j][k]
    #print( rhou_vis )
    #print( rhov_vis )
    #print( rhow_vis )
    #print( rhoE_vis )


### Sum inviscid and viscous fluxes
def sum_fluxes( rho_tot, rhou_tot, rhov_tot, rhow_tot, rhoE_tot, rho_inv, rhou_inv, rhov_inv, rhow_inv, rhoE_inv, rhou_vis, rhov_vis, rhow_vis, rhoE_vis, rk_iter ):

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                rho_tot[i][j][k][rk_iter]  = ( -1.0 )*rho_inv[i][j][k]
                rhou_tot[i][j][k][rk_iter] = ( -1.0 )*rhou_inv[i][j][k] + rhou_vis[i][j][k]
                rhov_tot[i][j][k][rk_iter] = ( -1.0 )*rhov_inv[i][j][k] + rhov_vis[i][j][k]
                rhow_tot[i][j][k][rk_iter] = ( -1.0 )*rhow_inv[i][j][k] + rhow_vis[i][j][k]
                rhoE_tot[i][j][k][rk_iter] = ( -1.0 )*rhoE_inv[i][j][k] + rhoE_vis[i][j][k]
    #print( rho_tot )
    #print( rhou_tot )
    #print( rhov_tot )
    #print( rhow_tot )
    #print( rhoE_tot )


### Time integration of conserved variables
def time_integration( y, y_0, h, k_s, rk_iter ):

    # Explicit third-order strong-stability-preserving Runge-Kutta (SSP-RK3) method:
    # S. Gottlieb, C.-W. Shu & E. Tadmor.
    # Strong stability-preserving high-order time discretization methods.
    # SIAM Review 43, 89-112, 2001.

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                if( rk_iter == 0 ):
                    y[i][j][k] = y_0[i][j][k] + h*k_s[i][j][k][0]
                elif( rk_iter == 1 ):
                    y[i][j][k] = y_0[i][j][k] + ( h/4.0 )*( k_s[i][j][k][0] + k_s[i][j][k][1] )
                elif( rk_iter == 2 ):
                    y[i][j][k] = y_0[i][j][k] + ( h/6.0 )*( k_s[i][j][k][0] + k_s[i][j][k][1] + 4.0*k_s[i][j][k][2] )
    #print( y )


### Update boundaries
def update_boundaries( rho, rhou, rhov, rhow, rhoE ):

    # Boundary conditions:
    # West:  transmissive
    # East:  transmissive
    # South: transmissive
    # North: transmissive
    # Back:  transmissive
    # Front: transmissive

    # West boundary points
    i = 0
    for j in range( 0, num_grid_y + 2 ):    
        for k in range( 0, num_grid_z + 2 ):    
            rho[i][j][k]  = rho[i+1][j][k]
            rhou[i][j][k] = rhou[i+1][j][k]
            rhov[i][j][k] = rhov[i+1][j][k]
            rhow[i][j][k] = rhow[i+1][j][k]
            rhoE[i][j][k] = rhoE[i+1][j][k]

    # East boundary points
    i = num_grid_x + 1    
    for j in range( 0, num_grid_y + 2 ):    
        for k in range( 0, num_grid_z + 2 ):    
            rho[i][j][k]  = rho[i-1][j][k]
            rhou[i][j][k] = rhou[i-1][j][k]
            rhov[i][j][k] = rhov[i-1][j][k]
            rhow[i][j][k] = rhow[i-1][j][k]
            rhoE[i][j][k] = rhoE[i-1][j][k]

    # South boundary points
    j = 0    
    for i in range( 0, num_grid_x + 2 ):    
        for k in range( 0, num_grid_z + 2 ):    
            rho[i][j][k]  = rho[i][j+1][k]
            rhou[i][j][k] = rhou[i][j+1][k]
            rhov[i][j][k] = rhov[i][j+1][k]
            rhow[i][j][k] = rhow[i][j+1][k]
            rhoE[i][j][k] = rhoE[i][j+1][k]

    # North boundary points
    j = num_grid_y + 1    
    for i in range( 0, num_grid_x + 2 ):    
        for k in range( 0, num_grid_z + 2 ):    
            rho[i][j][k]  = rho[i][j-1][k]            
            rhou[i][j][k] = rhou[i][j-1][k]            
            rhov[i][j][k] = rhov[i][j-1][k]            
            rhow[i][j][k] = rhow[i][j-1][k]            
            rhoE[i][j][k] = rhoE[i][j-1][k]            

    # Back boundary points
    k = 0    
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            rho[i][j][k]  = rho[i][j][k+1]
            rhou[i][j][k] = rhou[i][j][k+1]
            rhov[i][j][k] = rhov[i][j][k+1]
            rhow[i][j][k] = rhow[i][j][k+1]
            rhoE[i][j][k] = rhoE[i][j][k+1]

    # Front boundary points
    k = num_grid_z + 1    
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            rho[i][j][k]  = rho[i][j][k-1]
            rhou[i][j][k] = rhou[i][j][k-1]
            rhov[i][j][k] = rhov[i][j][k-1]
            rhow[i][j][k] = rhow[i][j][k-1]
            rhoE[i][j][k] = rhoE[i][j][k-1]
    #print( rho )    
    #print( rhou )    
    #print( rhov )    
    #print( rhow )    
    #print( rhoE )    


### Update primitive variables from conserved variables
def update_primitive( primitive, conserved, rho ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                primitive[i][j][k] = ( 1.0/rho[i][j][k] )*conserved[i][j][k]
    #print( primitive )


### Update thermodynamic variables from primitive variables
def thermodynamic_state( P, T, sos, rho, u, v, w, E ):

    # Ideal-gas model:
    # c_v = R_specific/(gamma - 1) is specific heat capcity at constant volume
    # ke = (u*u + v*v + w*w)/2 is specific kinetic energy
    # e = E - ke is specific internal energy
    # P = e*rho*(gamma - 1) is pressure
    # T = e/c_v is temperature
    # sos = sqrt(gamma*P/rho) is speed of sound

    # Specific heat capacity at constant volume
    c_v = R_specific/( gamma - 1.0 )

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                ke           = 0.5*( u[i][j][k]**2.0 + v[i][j][k]**2.0 + w[i][j][k]**2.0 )
                e            = E[i][j][k] - ke
                P[i][j][k]   = e*rho[i][j][k]*( gamma - 1.0 ) 
                T[i][j][k]   = ( 1.0/c_v )*e
                sos[i][j][k] = np.sqrt( gamma*( ( 1.0/rho[i][j][k] )*P[i][j][k] ) )
    #print( P )
    #print( T )
    #print( sos )


### Output data to file
def data_output( rho, u, v, w, E, P, T, sos, grid ):
    
    # Write 3-D data:
    # x, y, z, rho, u, v, w, E, P, T, sos

    # Open output file
    data_file_out = open( name_file_out, 'wt' )

    # Header string
    header_string = '# x [m], y [m], z [m], rho [kg/m3], u [m/s], v [m/s], w [m/s], E [J/kg], P [Pa], T [K], sos [m/s]\n'
    data_file_out.write( header_string )

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):  
                output_string  = str( grid[i][j][k][0] ) + ',' + str( grid[i][j][k][1] ) + ',' + str( grid[i][j][k][2] ) 
                output_string += ','
                output_string += str( rho[i][j][k] )
                output_string += ','
                output_string += str( u[i][j][k] ) + ',' + str( v[i][j][k] ) + ',' + str( w[i][j][k] )
                output_string += ','
                output_string += str( E[i][j][k] )
                output_string += ','
                output_string += str( P[i][j][k] ) + ',' + str( T[i][j][k] ) + ',' + str( sos[i][j][k] )
                output_string += '\n' 
                data_file_out.write( output_string )

    # Close output file
    data_file_out.close()



########## MAIN ##########

### START SIMULATION
print( 'RHEA: START SIMULATION' )

### Define spatial discretization
spatial_discretization( mesh )

### Initialize u, v, w, P and T variables
initialize_uvwPT( u_field, v_field, w_field, P_field, T_field, mesh )

### Initialize thermodynamic variables
initialize_thermodynamics( rho_field, E_field, sos_field, u_field, v_field, w_field, P_field, T_field )

### Update conserved variables from primitive variables
update_conserved( rhou_field, u_field, rho_field )
update_conserved( rhov_field, v_field, rho_field )
update_conserved( rhow_field, w_field, rho_field )
update_conserved( rhoE_field, E_field, rho_field )

### Update old fields of conserved variables
update_field( rho_0_field,  rho_field )  
update_field( rhou_0_field, rhou_field )  
update_field( rhov_0_field, rhov_field )  
update_field( rhow_0_field, rhow_field )  
update_field( rhoE_0_field, rhoE_field )

### Initialize time
time = initial_time

### Iterate solver in time
for t in range( 0, int( max_num_time_iter ) ):

    ### Calculate time step
    delta_t = time_step( rho_field, u_field, v_field, w_field, sos_field, mesh )
    if( ( time + delta_t ) > final_time ):
        delta_t = final_time - time
        #print( delta_t )

    ### Print time iteration iformation
    print( 'Time iteration ' + str( t ) + ': t = ' + str( time ) + ' [s], delta_t = ' + str( delta_t ) + ' [s]' )

    ### Runge-Kutta sub-steps
    for rk in range( 0, rk_order ):

        ### Calculate source terms
        source_terms( f_rhou_field, f_rhov_field, f_rhow_field, f_rhoE_field, rho_field, u_field, v_field, w_field, mesh )

        ### Calculate inviscid and viscous fluxes
        inviscid_fluxes( rho_inv_flux, rhou_inv_flux, rhov_inv_flux, rhow_inv_flux, rhoE_inv_flux, rho_field, u_field, v_field, w_field, E_field, P_field, sos_field, mesh )
        viscous_fluxes( rhou_vis_flux, rhov_vis_flux, rhow_vis_flux, rhoE_vis_flux, u_field, v_field, w_field, T_field, f_rhou_field, f_rhov_field, f_rhow_field, f_rhoE_field, mesh )
        sum_fluxes( rho_rk_fluxes, rhou_rk_fluxes, rhov_rk_fluxes, rhow_rk_fluxes, rhoE_rk_fluxes, rho_inv_flux, rhou_inv_flux, rhov_inv_flux, rhow_inv_flux, rhoE_inv_flux, rhou_vis_flux, rhov_vis_flux, rhow_vis_flux, rhoE_vis_flux, rk )

        ### Advance conserved variables in time
        time_integration( rho_field,  rho_0_field,  delta_t, rho_rk_fluxes,  rk )
        time_integration( rhou_field, rhou_0_field, delta_t, rhou_rk_fluxes, rk )
        time_integration( rhov_field, rhov_0_field, delta_t, rhov_rk_fluxes, rk )
        time_integration( rhow_field, rhow_0_field, delta_t, rhow_rk_fluxes, rk )
        time_integration( rhoE_field, rhoE_0_field, delta_t, rhoE_rk_fluxes, rk )
        update_boundaries( rho_field, rhou_field, rhov_field, rhow_field, rhoE_field )

        ### Update primitive variables from conserved variables
        update_primitive( u_field, rhou_field, rho_field )
        update_primitive( v_field, rhov_field, rho_field )
        update_primitive( w_field, rhow_field, rho_field )
        update_primitive( E_field, rhoE_field, rho_field )
       
        ### Update thermodynamic variables from primitive variables
        thermodynamic_state( P_field, T_field, sos_field, rho_field, u_field, v_field, w_field, E_field )

    ### Update old fields of conserved variables
    update_field( rho_0_field,  rho_field )  
    update_field( rhou_0_field, rhou_field )  
    update_field( rhov_0_field, rhov_field )  
    update_field( rhow_0_field, rhow_field )  
    update_field( rhoE_0_field, rhoE_field )  

    ### Update time
    time += delta_t

    ### Check if simulation is completed (time > final_time)
    if( time >= final_time ):
        break

### Output data to file
data_output( rho_field, u_field, v_field, w_field, E_field, P_field, T_field, sos_field, mesh )

### Print data output information
print( 'Data output at time' + ': t = ' + str( time ) + ' [s]' )

### END SIMULATION
print( 'RHEA: END SIMULATION' )

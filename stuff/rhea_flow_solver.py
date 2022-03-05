#!/usr/bin/python

################################################################################################################
#                                                                                                              #
# RHEA - an open-source Reproducible Hybrid-architecture flow solver Engineered for Academia                  #
#                                                                                                              #
# Rhea was the Titaness great Mother of the Gods, and goddess of female fertility, motherhood, and generation. #
# Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.           #
#                                                                                                              #
#                                                                                                              #
# REHA is released under the MIT License:                                                                      #
#                                                                                                              #
# Copyright (c) 2020 Lluis Jofre Cruanyes & Guillermo Oyarzun Altamirano.                                      #
#                                                                                                              #
# Permission is hereby granted, free of charge, to any person obtaining a copy                                 #
# of this software and associated documentation files (the "Software"), to deal                                #
# in the Software without restriction, including without limitation the rights                                 #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell                                    #
# copies of the Software, and to permit persons to whom the Software is                                        #
# furnished to do so, subject to the following conditions:                                                     #
#                                                                                                              #
# The above copyright notice and this permission notice shall be included in all                               #
# copies or substantial portions of the Software.                                                              #
#                                                                                                              #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR                                   #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                                     #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE                                  #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                                       #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,                                #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE                                #
# SOFTWARE.                                                                                                    #
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
R_specific  = 287.058					# Specific gas constant of the fluid [J/(kg K)]
gamma       = 1.4					# Heat capacity ratio of the fluid [-]
c_p         = gamma*R_specific/( gamma - 1.0 )		# Isobaric heat capacity
Re_L        = 100.0					# Reynolds number based on L [-]
Ma          = 1.0e-1/np.sqrt( gamma )			# Mach number [-]
Pr          = 0.71					# Prandtl number [-]
rho_0       = 1.0					# Reference density [kg/m^3]
L           = 1.0					# Cavity size [m]
U_lid       = 1.0 					# Lid velocity [m/s]
P_0         = rho_0*U_lid*U_lid/( gamma*Ma*Ma )		# Reference pressure [Pa] 
T_0         = ( 1.0/( rho_0*R_specific ) )*P_0		# Reference temperature [K] 
mu_value    = rho_0*U_lid*L/Re_L       			# Dynamic viscosity of the fluid [Pa s]
kappa_value = c_p*mu_value/Pr				# Thermal conductivity of the fluid [W/(m k)]

### Problem parameters
x_0           = 0.0                     # Domain origin in x-direction [m]
y_0           = 0.0                     # Domain origin in y-direction [m]
z_0           = 0.0                     # Domain origin in z-direction [m]
L_x           = L                  	# Size of domain in x-direction
L_y           = L	      		# Size of domain in y-direction
L_z           = 0.1*L_x         	# Size of domain in z-direction
initial_time  = 0.0   			# Initial time [s]
final_time    = 50.0      		# Final time [s]
name_file_out = 'output_data.csv'	# Name of output data [-]

### Computational parameters
num_grid_x        = 32 			# Number of internal grid points in the x-direction
num_grid_y        = 32			# Number of internal grid points in the y-direction
num_grid_z        = 1			# Number of internal grid points in the z-direction
# Stretching factors: x = L*eta + A*( 0.5*L - L*eta )*( 1.0 - eta )*eta, with eta = ( l - 0.5 )/num_grid 
# A < 0: stretching at ends; A = 0: uniform; A > 0: stretching at center
A_x               = -1.0                # Stretching factor in x-direction
A_y               = -1.0                # Stretching factor in y-direction
A_z               = 0.0                 # Stretching factor in z-direction
CFL               = 0.9			# CFL coefficient
max_num_time_iter = 1e6			# Maximum number of time iterations
output_iter       = 1e2			# Output data every given number of iterations

### Fixed parameters
num_sptl_dim = 3         		# Number of spatial dimensions (fixed value)
rk_order     = 3         		# Time-integration Runge-Kutta order (fixed value)
epsilon      = 1.0e-15   		# Small epsilon number (fixed value)


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

### Transport coefficients ... two positions added for boundary points
mu_field    = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                    # 3-D field of mu
kappa_field = np.zeros( [ num_grid_x + 2, num_grid_y + 2, num_grid_z + 2 ] )                    # 3-D field of kappa

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

### Initialize u, v, w, P and T variables
def initialize_uvwPT( u, v, w, P, T, grid ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                u[i][j][k] = 0.0
                v[i][j][k] = 0.0
                w[i][j][k] = 0.0
                P[i][j][k] = P_0
                T[i][j][k] = T_0
    #print( u )
    #print( v )
    #print( w )
    #print( P )
    #print( T )


### Update boundaries
def update_boundaries( rho, rhou, rhov, rhow, rhoE, u, v, w, P, T, grid ):

    # General form: w_g*phi_g + w_in*phi_in = phi_b
    # phi_g is ghost cell value
    # phi_in is inner cell value
    # phi_b is boundary value/flux
    # w_g is ghost cell weight
    # w_in is inner cell weight

    # Boundary conditions:
    # West:  Dirichlet & Neumann
    # East:  Dirichlet & Neumann
    # South: Dirichlet & Neumann
    # North: Dirichlet & Neumann
    # Back:  Periodic
    # Front: Periodic

    # West boundary points
    i = 0
    for j in range( 1, num_grid_y + 1 ):    
        for k in range( 1, num_grid_z + 1 ):
            wg_g  = 1.0 - ( x_0 - grid[i][j][k][0] )/( grid[i+1][j][k][0] - grid[i][j][k][0] );
            wg_in = 1.0 - ( grid[i+1][j][k][0] - x_0 )/( grid[i+1][j][k][0] - grid[i][j][k][0] ); 
            P_in  = P[i+1][j][k]    
            T_in  = T[i+1][j][k]    
            u_in  = rhou[i+1][j][k]/rho[i+1][j][k]
            v_in  = rhov[i+1][j][k]/rho[i+1][j][k]
            w_in  = rhow[i+1][j][k]/rho[i+1][j][k]
            P_g   = P_in;					# Neumann
            T_g   = ( T_0 - wg_in*T_in )/wg_g;			# Dirichlet
            u_g   = ( 0.0 - wg_in*u_in )/wg_g;			# Dirichlet
            v_g   = ( 0.0 - wg_in*v_in )/wg_g;			# Dirichlet
            w_g   = ( 0.0 - wg_in*w_in )/wg_g;			# Dirichlet
            rho_g = ( 1.0/( R_specific*T_g ) )*P_g		# Ideal-gas law
            e_g   = ( 1.0/( rho_g*( gamma - 1.0 ) ) )*P_g	# Specific internal energy
            ke_g  = 0.5*( u_g**2.0 + v_g**2.0 + w_g**2.0 )	# Specific kinetic energy
            E_g   = e_g + ke_g					# Specific total energy
            rho[i][j][k]  = rho_g
            rhou[i][j][k] = rho_g*u_g
            rhov[i][j][k] = rho_g*v_g
            rhow[i][j][k] = rho_g*w_g
            rhoE[i][j][k] = rho_g*E_g

    # East boundary points
    i = num_grid_x + 1    
    for j in range( 1, num_grid_y + 1 ):    
        for k in range( 1, num_grid_z + 1 ):    
            wg_g  = 1.0 - ( grid[i][j][k][0] - ( x_0 + L ) )/( grid[i][j][k][0] - grid[i-1][j][k][0] );
            wg_in = 1.0 - ( ( x_0 + L ) - grid[i-1][j][k][0] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ); 
            P_in  = P[i-1][j][k]    
            T_in  = T[i-1][j][k]    
            u_in  = rhou[i-1][j][k]/rho[i-1][j][k]
            v_in  = rhov[i-1][j][k]/rho[i-1][j][k]
            w_in  = rhow[i-1][j][k]/rho[i-1][j][k]
            P_g   = P_in;					# Neumann
            T_g   = ( T_0 - wg_in*T_in )/wg_g;			# Dirichlet
            u_g   = ( 0.0 - wg_in*u_in )/wg_g;			# Dirichlet
            v_g   = ( 0.0 - wg_in*v_in )/wg_g;			# Dirichlet
            w_g   = ( 0.0 - wg_in*w_in )/wg_g;			# Dirichlet
            rho_g = ( 1.0/( R_specific*T_g ) )*P_g		# Ideal-gas law
            e_g   = ( 1.0/( rho_g*( gamma - 1.0 ) ) )*P_g	# Specific internal energy
            ke_g  = 0.5*( u_g**2.0 + v_g**2.0 + w_g**2.0 )	# Specific kinetic energy
            E_g   = e_g + ke_g					# Specific total energy
            rho[i][j][k]  = rho_g
            rhou[i][j][k] = rho_g*u_g
            rhov[i][j][k] = rho_g*v_g
            rhow[i][j][k] = rho_g*w_g
            rhoE[i][j][k] = rho_g*E_g

    # South boundary points
    j = 0    
    for i in range( 1, num_grid_x + 1 ):    
        for k in range( 1, num_grid_z + 1 ):    
            wg_g  = 1.0 - ( y_0 - grid[i][j][k][1] )/( grid[i][j+1][k][1] - grid[i][j][k][1] );
            wg_in = 1.0 - ( grid[i][j+1][k][1] - y_0 )/( grid[i][j+1][k][1] - grid[i][j][k][1] ); 
            P_in  = P[i][j+1][k]    
            T_in  = T[i][j+1][k]    
            u_in  = rhou[i][j+1][k]/rho[i][j+1][k]
            v_in  = rhov[i][j+1][k]/rho[i][j+1][k]
            w_in  = rhow[i][j+1][k]/rho[i][j+1][k]
            P_g   = P_in;					# Neumann
            T_g   = ( T_0 - wg_in*T_in )/wg_g;			# Dirichlet
            u_g   = ( 0.0 - wg_in*u_in )/wg_g;		        # Dirichlet
            v_g   = ( 0.0 - wg_in*v_in )/wg_g;			# Dirichlet
            w_g   = ( 0.0 - wg_in*w_in )/wg_g;			# Dirichlet
            rho_g = ( 1.0/( R_specific*T_g ) )*P_g		# Ideal-gas law
            e_g   = ( 1.0/( rho_g*( gamma - 1.0 ) ) )*P_g	# Specific internal energy
            ke_g  = 0.5*( u_g**2.0 + v_g**2.0 + w_g**2.0 )	# Specific kinetic energy
            E_g   = e_g + ke_g					# Specific total energy
            rho[i][j][k]  = rho_g
            rhou[i][j][k] = rho_g*u_g
            rhov[i][j][k] = rho_g*v_g
            rhow[i][j][k] = rho_g*w_g
            rhoE[i][j][k] = rho_g*E_g

    # North boundary points
    j = num_grid_y + 1    
    for i in range( 1, num_grid_x + 1 ):    
        for k in range( 1, num_grid_z + 1 ):    
            wg_g  = 1.0 - ( grid[i][j][k][1] - ( y_0 + L ) )/( grid[i][j][k][1] - grid[i][j-1][k][1] );
            wg_in = 1.0 - ( ( y_0 + L ) - grid[i][j-1][k][1] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ); 
            P_in  = P[i][j-1][k]    
            T_in  = T[i][j-1][k]    
            u_in  = rhou[i][j-1][k]/rho[i][j-1][k]
            v_in  = rhov[i][j-1][k]/rho[i][j-1][k]
            w_in  = rhow[i][j-1][k]/rho[i][j-1][k]
            P_g   = P_in;					# Neumann
            T_g   = ( T_0 - wg_in*T_in )/wg_g;			# Dirichlet
            u_g   = ( U_lid - wg_in*u_in )/wg_g;		# Dirichlet
            v_g   = ( 0.0 - wg_in*v_in )/wg_g;			# Dirichlet
            w_g   = ( 0.0 - wg_in*w_in )/wg_g;			# Dirichlet
            rho_g = ( 1.0/( R_specific*T_g ) )*P_g		# Ideal-gas law
            e_g   = ( 1.0/( rho_g*( gamma - 1.0 ) ) )*P_g	# Specific internal energy
            ke_g  = 0.5*( u_g**2.0 + v_g**2.0 + w_g**2.0 )	# Specific kinetic energy
            E_g   = e_g + ke_g					# Specific total energy
            rho[i][j][k]  = rho_g
            rhou[i][j][k] = rho_g*u_g
            rhov[i][j][k] = rho_g*v_g
            rhow[i][j][k] = rho_g*w_g
            rhoE[i][j][k] = rho_g*E_g            

    # Back boundary points
    k = 0    
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            rho[i][j][k]  = rho[i][j][k+1]
            rhou[i][j][k] = rhou[i][j][k+1]
            rhov[i][j][k] = rhov[i][j][k+1]
            rhow[i][j][k] = rhow[i][j][k+1]
            rhoE[i][j][k] = rhoE[i][j][k+1]

    # Front boundary points
    k = num_grid_z + 1    
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            rho[i][j][k]  = rho[i][j][k-1]
            rhou[i][j][k] = rhou[i][j][k-1]
            rhov[i][j][k] = rhov[i][j][k-1]
            rhow[i][j][k] = rhow[i][j][k-1]
            rhoE[i][j][k] = rhoE[i][j][k-1]

    # Fill x-direction edge boundary points
    for i in range( 1, num_grid_x + 1 ):    
        j = 0; k = 0
        rho[i][j][k]  = 0.5*( rho[i][j+1][k] + rho[i][j][k+1] )
        rhou[i][j][k] = 0.5*( rhou[i][j+1][k] + rhou[i][j][k+1] )
        rhov[i][j][k] = 0.5*( rhov[i][j+1][k] + rhov[i][j][k+1] )
        rhow[i][j][k] = 0.5*( rhow[i][j+1][k] + rhow[i][j][k+1] )
        rhoE[i][j][k] = 0.5*( rhoE[i][j+1][k] + rhoE[i][j][k+1] )
        j = 0; k = num_grid_z + 1
        rho[i][j][k]  = 0.5*( rho[i][j+1][k] + rho[i][j][k-1] )
        rhou[i][j][k] = 0.5*( rhou[i][j+1][k] + rhou[i][j][k-1] )
        rhov[i][j][k] = 0.5*( rhov[i][j+1][k] + rhov[i][j][k-1] )
        rhow[i][j][k] = 0.5*( rhow[i][j+1][k] + rhow[i][j][k-1] )
        rhoE[i][j][k] = 0.5*( rhoE[i][j+1][k] + rhoE[i][j][k-1] )
        j = num_grid_y + 1; k = 0
        rho[i][j][k]  = 0.5*( rho[i][j-1][k] + rho[i][j][k+1] )
        rhou[i][j][k] = 0.5*( rhou[i][j-1][k] + rhou[i][j][k+1] )
        rhov[i][j][k] = 0.5*( rhov[i][j-1][k] + rhov[i][j][k+1] )
        rhow[i][j][k] = 0.5*( rhow[i][j-1][k] + rhow[i][j][k+1] )
        rhoE[i][j][k] = 0.5*( rhoE[i][j-1][k] + rhoE[i][j][k+1] )
        j = num_grid_y + 1; k = num_grid_z + 1
        rho[i][j][k]  = 0.5*( rho[i][j-1][k] + rho[i][j][k-1] )
        rhou[i][j][k] = 0.5*( rhou[i][j-1][k] + rhou[i][j][k-1] )
        rhov[i][j][k] = 0.5*( rhov[i][j-1][k] + rhov[i][j][k-1] )
        rhow[i][j][k] = 0.5*( rhow[i][j-1][k] + rhow[i][j][k-1] )
        rhoE[i][j][k] = 0.5*( rhoE[i][j-1][k] + rhoE[i][j][k-1] )

    # Fill y-direction edge boundary points
    for j in range( 1, num_grid_y + 1 ):    
        i = 0; k = 0
        rho[i][j][k]  = 0.5*( rho[i+1][j][k] + rho[i][j][k+1] )
        rhou[i][j][k] = 0.5*( rhou[i+1][j][k] + rhou[i][j][k+1] )
        rhov[i][j][k] = 0.5*( rhov[i+1][j][k] + rhov[i][j][k+1] )
        rhow[i][j][k] = 0.5*( rhow[i+1][j][k] + rhow[i][j][k+1] )
        rhoE[i][j][k] = 0.5*( rhoE[i+1][j][k] + rhoE[i][j][k+1] )
        i = 0; k = num_grid_z + 1
        rho[i][j][k]  = 0.5*( rho[i+1][j][k] + rho[i][j][k-1] )
        rhou[i][j][k] = 0.5*( rhou[i+1][j][k] + rhou[i][j][k-1] )
        rhov[i][j][k] = 0.5*( rhov[i+1][j][k] + rhov[i][j][k-1] )
        rhow[i][j][k] = 0.5*( rhow[i+1][j][k] + rhow[i][j][k-1] )
        rhoE[i][j][k] = 0.5*( rhoE[i+1][j][k] + rhoE[i][j][k-1] )
        i = num_grid_x + 1; k = 0
        rho[i][j][k]  = 0.5*( rho[i-1][j][k] + rho[i][j][k+1] )
        rhou[i][j][k] = 0.5*( rhou[i-1][j][k] + rhou[i][j][k+1] )
        rhov[i][j][k] = 0.5*( rhov[i-1][j][k] + rhov[i][j][k+1] )
        rhow[i][j][k] = 0.5*( rhow[i-1][j][k] + rhow[i][j][k+1] )
        rhoE[i][j][k] = 0.5*( rhoE[i-1][j][k] + rhoE[i][j][k+1] )
        i = num_grid_x + 1; k = num_grid_z + 1
        rho[i][j][k]  = 0.5*( rho[i-1][j][k] + rho[i][j][k-1] )
        rhou[i][j][k] = 0.5*( rhou[i-1][j][k] + rhou[i][j][k-1] )
        rhov[i][j][k] = 0.5*( rhov[i-1][j][k] + rhov[i][j][k-1] )
        rhow[i][j][k] = 0.5*( rhow[i-1][j][k] + rhow[i][j][k-1] )
        rhoE[i][j][k] = 0.5*( rhoE[i-1][j][k] + rhoE[i][j][k-1] )

    # Fill z-direction edge boundary points
    for k in range( 1, num_grid_z + 1 ):    
        i = 0; j = 0
        rho[i][j][k]  = 0.5*( rho[i+1][j][k] + rho[i][j+1][k] )
        rhou[i][j][k] = 0.5*( rhou[i+1][j][k] + rhou[i][j+1][k] )
        rhov[i][j][k] = 0.5*( rhov[i+1][j][k] + rhov[i][j+1][k] )
        rhow[i][j][k] = 0.5*( rhow[i+1][j][k] + rhow[i][j+1][k] )
        rhoE[i][j][k] = 0.5*( rhoE[i+1][j][k] + rhoE[i][j+1][k] )
        i = 0; j = num_grid_y + 1
        rho[i][j][k]  = 0.5*( rho[i+1][j][k] + rho[i][j-1][k] )
        rhou[i][j][k] = 0.5*( rhou[i+1][j][k] + rhou[i][j-1][k] )
        rhov[i][j][k] = 0.5*( rhov[i+1][j][k] + rhov[i][j-1][k] )
        rhow[i][j][k] = 0.5*( rhow[i+1][j][k] + rhow[i][j-1][k] )
        rhoE[i][j][k] = 0.5*( rhoE[i+1][j][k] + rhoE[i][j-1][k] )
        i = num_grid_x + 1; j = 0
        rho[i][j][k]  = 0.5*( rho[i-1][j][k] + rho[i][j+1][k] )
        rhou[i][j][k] = 0.5*( rhou[i-1][j][k] + rhou[i][j+1][k] )
        rhov[i][j][k] = 0.5*( rhov[i-1][j][k] + rhov[i][j+1][k] )
        rhow[i][j][k] = 0.5*( rhow[i-1][j][k] + rhow[i][j+1][k] )
        rhoE[i][j][k] = 0.5*( rhoE[i-1][j][k] + rhoE[i][j+1][k] )
        i = num_grid_x + 1; j = num_grid_y + 1
        rho[i][j][k]  = 0.5*( rho[i-1][j][k] + rho[i][j-1][k] )
        rhou[i][j][k] = 0.5*( rhou[i-1][j][k] + rhou[i][j-1][k] )
        rhov[i][j][k] = 0.5*( rhov[i-1][j][k] + rhov[i][j-1][k] )
        rhow[i][j][k] = 0.5*( rhow[i-1][j][k] + rhow[i][j-1][k] )
        rhoE[i][j][k] = 0.5*( rhoE[i-1][j][k] + rhoE[i][j-1][k] )

    # Fill corner boundary points
    i = 0; j = 0; k = 0
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i+1][j][k] + rho[i][j+1][k] + rho[i][j][k+1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i+1][j][k] + rhou[i][j+1][k] + rhou[i][j][k+1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i+1][j][k] + rhov[i][j+1][k] + rhov[i][j][k+1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i+1][j][k] + rhow[i][j+1][k] + rhow[i][j][k+1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i+1][j][k] + rhoE[i][j+1][k] + rhoE[i][j][k+1] )
    i = num_grid_x + 1; j = 0; k = 0
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i-1][j][k] + rho[i][j+1][k] + rho[i][j][k+1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i-1][j][k] + rhou[i][j+1][k] + rhou[i][j][k+1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i-1][j][k] + rhov[i][j+1][k] + rhov[i][j][k+1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i-1][j][k] + rhow[i][j+1][k] + rhow[i][j][k+1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i-1][j][k] + rhoE[i][j+1][k] + rhoE[i][j][k+1] )
    i = 0; j = num_grid_y + 1; k = 0
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i+1][j][k] + rho[i][j-1][k] + rho[i][j][k+1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i+1][j][k] + rhou[i][j-1][k] + rhou[i][j][k+1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i+1][j][k] + rhov[i][j-1][k] + rhov[i][j][k+1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i+1][j][k] + rhow[i][j-1][k] + rhow[i][j][k+1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i+1][j][k] + rhoE[i][j-1][k] + rhoE[i][j][k+1] )
    i = num_grid_x + 1; j = num_grid_y + 1; k = 0
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i-1][j][k] + rho[i][j-1][k] + rho[i][j][k+1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i-1][j][k] + rhou[i][j-1][k] + rhou[i][j][k+1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i-1][j][k] + rhov[i][j-1][k] + rhov[i][j][k+1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i-1][j][k] + rhow[i][j-1][k] + rhow[i][j][k+1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i-1][j][k] + rhoE[i][j-1][k] + rhoE[i][j][k+1] )
    i = 0; j = 0; k = num_grid_z + 1
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i+1][j][k] + rho[i][j+1][k] + rho[i][j][k-1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i+1][j][k] + rhou[i][j+1][k] + rhou[i][j][k-1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i+1][j][k] + rhov[i][j+1][k] + rhov[i][j][k-1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i+1][j][k] + rhow[i][j+1][k] + rhow[i][j][k-1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i+1][j][k] + rhoE[i][j+1][k] + rhoE[i][j][k-1] )
    i = num_grid_x + 1; j = 0; k = num_grid_z + 1
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i-1][j][k] + rho[i][j+1][k] + rho[i][j][k-1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i-1][j][k] + rhou[i][j+1][k] + rhou[i][j][k-1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i-1][j][k] + rhov[i][j+1][k] + rhov[i][j][k-1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i-1][j][k] + rhow[i][j+1][k] + rhow[i][j][k-1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i-1][j][k] + rhoE[i][j+1][k] + rhoE[i][j][k-1] )
    i = 0; j = num_grid_y + 1; k = num_grid_z + 1
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i+1][j][k] + rho[i][j-1][k] + rho[i][j][k-1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i+1][j][k] + rhou[i][j-1][k] + rhou[i][j][k-1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i+1][j][k] + rhov[i][j-1][k] + rhov[i][j][k-1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i+1][j][k] + rhow[i][j-1][k] + rhow[i][j][k-1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i+1][j][k] + rhoE[i][j-1][k] + rhoE[i][j][k-1] )
    i = num_grid_x + 1; j = num_grid_y + 1; k = num_grid_z + 1
    rho[i][j][k]  = ( 1.0/3.0 )*( rho[i-1][j][k] + rho[i][j-1][k] + rho[i][j][k-1] )
    rhou[i][j][k] = ( 1.0/3.0 )*( rhou[i-1][j][k] + rhou[i][j-1][k] + rhou[i][j][k-1] )
    rhov[i][j][k] = ( 1.0/3.0 )*( rhov[i-1][j][k] + rhov[i][j-1][k] + rhov[i][j][k-1] )
    rhow[i][j][k] = ( 1.0/3.0 )*( rhow[i-1][j][k] + rhow[i][j-1][k] + rhow[i][j][k-1] )
    rhoE[i][j][k] = ( 1.0/3.0 )*( rhoE[i-1][j][k] + rhoE[i][j-1][k] + rhoE[i][j][k-1] )

    #print( rho )    
    #print( rhou )    
    #print( rhov )    
    #print( rhow )    
    #print( rhoE )


### Calculate transport coefficients
def transport_coefficients( mu, kappa ):

    # Constant-value model:

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                mu[i][j][k]    = mu_value        
                kappa[i][j][k] = kappa_value
    #print( mu )
    #print( kappa )


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
    #print( f_rhou )
    #print( f_rhov )
    #print( f_rhow )
    #print( f_rhoE )


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


### Calculate time step
def time_step( rho, u, v, w, sos, mu, kappa, grid ):

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
                ## Calculate Prandtl (Pr) number
                Pr = gamma
                if( kappa[i][j][k] > epsilon ):
                    Pr = ( 1.0/( kappa[i][j][k] ) )*c_p*mu[i][j][k]
                ## x-direction inviscid & viscous terms
                S_x     = abs( u[i][j][k] ) + sos[i][j][k]
                delta_t = min( delta_t, ( 1.0/S_x )*CFL*delta_x )
                delta_t = min( delta_t, ( 1.0/( mu[i][j][k]*gamma + epsilon ) )*CFL*Pr*rho[i][j][k]*( delta_x**2.0 ) )
                ## y-direction inviscid & viscous terms
                S_y     = abs( v[i][j][k] ) + sos[i][j][k]
                delta_t = min( delta_t, ( 1.0/S_y )*CFL*delta_y )
                delta_t = min( delta_t, ( 1.0/( mu[i][j][k]*gamma + epsilon ) )*CFL*Pr*rho[i][j][k]*( delta_y**2.0 ) )
                ## z-direction inviscid & viscous terms
                S_z     = abs( w[i][j][k] ) + sos[i][j][k]
                delta_t = min( delta_t, ( 1.0/S_z )*CFL*delta_z )
                delta_t = min( delta_t, ( 1.0/( mu[i][j][k]*gamma + epsilon ) )*CFL*Pr*rho[i][j][k]*( delta_z**2.0 ) )
    #print( delta_t )

    # Return minimum time step
    return( delta_t )


### Define centroids of spatial discretization
def spatial_discretization( grid ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                # Evenly-spaced (dummy variable) cell centroids
                eta_x = ( i - 0.5 )/num_grid_x
                eta_y = ( j - 0.5 )/num_grid_y
                eta_z = ( k - 0.5 )/num_grid_z
                # Unevenly-spaced (stretched) cell centroids
                grid[i][j][k][0] = x_0 + L_x*eta_x + A_x*( 0.5*L_x - L_x*eta_x )*( 1.0 - eta_x )*eta_x
                grid[i][j][k][1] = y_0 + L_y*eta_y + A_y*( 0.5*L_y - L_y*eta_y )*( 1.0 - eta_y )*eta_y
                grid[i][j][k][2] = z_0 + L_z*eta_z + A_z*( 0.5*L_z - L_z*eta_z )*( 1.0 - eta_z )*eta_z
                # Adjust (symmetric) boundary cell centroids
                if( grid[i][j][k][0] < x_0 ):
                    eta_x = ( 1.0 - 0.5 )/num_grid_x
                    grid[i][j][k][0] = x_0 - ( L_x*eta_x + A_x*( 0.5*L_x - L_x*eta_x )*( 1.0 - eta_x )*eta_x )
                if( grid[i][j][k][0] > ( x_0 + L_x ) ):
                    eta_x = ( num_grid_x - 0.5 )/num_grid_x
                    grid[i][j][k][0] = x_0 + 2.0*L_x - ( L_x*eta_x + A_x*( 0.5*L_x - L_x*eta_x )*( 1.0 - eta_x )*eta_x )
                if( grid[i][j][k][1] < y_0 ):
                    eta_y = ( 1.0 - 0.5 )/num_grid_y
                    grid[i][j][k][1] = y_0 - ( L_y*eta_y + A_y*( 0.5*L_y - L_y*eta_y )*( 1.0 - eta_y )*eta_y )
                if( grid[i][j][k][1] > ( y_0 + L_y ) ):
                    eta_y = ( num_grid_y - 0.5 )/num_grid_y
                    grid[i][j][k][1] = y_0 + 2.0*L_y - ( L_y*eta_y + A_y*( 0.5*L_y - L_y*eta_y )*( 1.0 - eta_y )*eta_y )
                if( grid[i][j][k][2] < z_0 ):
                    eta_z = ( 1.0 - 0.5 )/num_grid_z
                    grid[i][j][k][2] = z_0 - ( L_z*eta_z + A_z*( 0.5*L_z - L_z*eta_z )*( 1.0 - eta_z )*eta_z )
                if( grid[i][j][k][2] > ( z_0 + L_z ) ):
                    eta_z = ( num_grid_z - 0.5 )/num_grid_z
                    grid[i][j][k][2] = z_0 + 2.0*L_z - ( L_z*eta_z + A_z*( 0.5*L_z - L_z*eta_z )*( 1.0 - eta_z )*eta_z )
    #print( grid )


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
def HLLC_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type ):

    # HLLC approximate Riemann solver:
    # E. F. Toro.
    # Riemann solvers and numerical methods for fluid dynamics.
    # Springer, 2009.    

    F_L = rho_L*u_L
    F_R = rho_R*u_R
    U_L = rho_L
    U_R = rho_R
    if( var_type == 0 ):
        F_L *= 1.0
        F_R *= 1.0
        U_L *= 1.0
        U_R *= 1.0
    elif( var_type == 1 ):
        F_L *= u_L; F_L += P_L
        F_R *= u_R; F_R += P_R
        U_L *= u_L
        U_R *= u_R
    elif( var_type == 2 ):
        F_L *= v_L
        F_R *= v_R
        U_L *= v_L
        U_R *= v_R
    elif( var_type == 3 ):
        F_L *= w_L
        F_R *= w_R
        U_L *= w_L
        U_R *= w_R
    elif( var_type == 4 ):
        F_L *= E_L; F_L += u_L*P_L
        F_R *= E_R; F_R += u_R*P_R
        U_L *= E_L
        U_R *= E_R

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


### Calculate KGP flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE 
def KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type ):

    # Kennedy, Gruber & Pirozzoli (KGP) scheme:
    # G. Coppola , F. Capuano , S. Pirozzoli, L. de Luca.
    # Numerically stable formulations of convective terms for turbulent compressible flows.
    # Journal of Computational Physics, 382, 86-104, 2019.

    F = ( 1.0/8.0 )*( rho_L + rho_R )*( u_L + u_R )
    if( var_type == 0 ):
        F *= 1.0 + 1.0;
    elif ( var_type == 1 ):
        F *= u_L + u_R
        F += ( 1.0/2.0 )*( P_L + P_R )
    elif ( var_type == 2 ):
        F *= v_L + v_R
    elif ( var_type == 3 ):
        F *= w_L + w_R
    elif ( var_type == 4 ):
        F *= E_L + P_L/rho_L + E_R + P_R/rho_R

    return( F )


### Calculate inviscid fluxes
def inviscid_fluxes( rho_inv, rhou_inv, rhov_inv, rhow_inv, rhoE_inv, rho, u, v, w, E, P, sos, grid ):

    # Unsplit method for Euler equations:
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
                rho_F_p  = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 1
                rhou_F_p = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_p = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_p = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_p = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )                
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
                rho_F_m  = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 1
                rhou_F_m = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_m = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_m = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_m = KGP_flux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                ## Fluxes x-direction
                rho_inv[i][j][k]  = ( 1.0/delta_x )*( rho_F_p - rho_F_m )
                rhou_inv[i][j][k] = ( 1.0/delta_x )*( rhou_F_p - rhou_F_m )
                rhov_inv[i][j][k] = ( 1.0/delta_x )*( rhov_F_p - rhov_F_m )
                rhow_inv[i][j][k] = ( 1.0/delta_x )*( rhow_F_p - rhow_F_m )
                rhoE_inv[i][j][k] = ( 1.0/delta_x )*( rhoE_F_p - rhoE_F_m )
                ## y-direction j+1/2
                index_L = j;                  index_R = j + 1
                rho_L   = rho[i][index_L][k]; rho_R   = rho[i][index_R][k] 
                u_L     = u[i][index_L][k];   u_R     = u[i][index_R][k]
                v_L     = v[i][index_L][k];   v_R     = v[i][index_R][k]
                w_L     = w[i][index_L][k];   w_R     = w[i][index_R][k]
                E_L     = E[i][index_L][k];   E_R     = E[i][index_R][k]
                P_L     = P[i][index_L][k];   P_R     = P[i][index_R][k]
                a_L     = sos[i][index_L][k]; a_R     = sos[i][index_R][k]
                # rho
                var_type = 0
                rho_F_p  = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 2
                rhou_F_p = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 1
                rhov_F_p = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_p = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_p = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )                
                ## y-direction j-1/2
                index_L = j - 1;              index_R = j
                rho_L   = rho[i][index_L][k]; rho_R   = rho[i][index_R][k] 
                u_L     = u[i][index_L][k];   u_R     = u[i][index_R][k]
                v_L     = v[i][index_L][k];   v_R     = v[i][index_R][k]
                w_L     = w[i][index_L][k];   w_R     = w[i][index_R][k]
                E_L     = E[i][index_L][k];   E_R     = E[i][index_R][k]
                P_L     = P[i][index_L][k];   P_R     = P[i][index_R][k]
                a_L     = sos[i][index_L][k]; a_R     = sos[i][index_R][k]
                # rho
                var_type = 0
                rho_F_m  = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 2
                rhou_F_m = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 1
                rhov_F_m = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 3
                rhow_F_m = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_m = KGP_flux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                ## Fluxes y-direction
                rho_inv[i][j][k]  += ( 1.0/delta_y )*( rho_F_p - rho_F_m )
                rhou_inv[i][j][k] += ( 1.0/delta_y )*( rhou_F_p - rhou_F_m )
                rhov_inv[i][j][k] += ( 1.0/delta_y )*( rhov_F_p - rhov_F_m )
                rhow_inv[i][j][k] += ( 1.0/delta_y )*( rhow_F_p - rhow_F_m )
                rhoE_inv[i][j][k] += ( 1.0/delta_y )*( rhoE_F_p - rhoE_F_m )
                ## z-direction k+1/2
                index_L = k;                  index_R = k + 1
                rho_L   = rho[i][j][index_L]; rho_R   = rho[i][j][index_R] 
                u_L     = u[i][j][index_L];   u_R     = u[i][j][index_R]
                v_L     = v[i][j][index_L];   v_R     = v[i][j][index_R]
                w_L     = w[i][j][index_L];   w_R     = w[i][j][index_R]
                E_L     = E[i][j][index_L];   E_R     = E[i][j][index_R]
                P_L     = P[i][j][index_L];   P_R     = P[i][j][index_R]
                a_L     = sos[i][j][index_L]; a_R     = sos[i][j][index_R]
                # rho
                var_type = 0
                rho_F_p  = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 3
                rhou_F_p = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_p = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 1
                rhow_F_p = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_p = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )                
                ## z-direction k-1/2
                index_L = k - 1;              index_R = k
                rho_L   = rho[i][j][index_L]; rho_R   = rho[i][j][index_R] 
                u_L     = u[i][j][index_L];   u_R     = u[i][j][index_R]
                v_L     = v[i][j][index_L];   v_R     = v[i][j][index_R]
                w_L     = w[i][j][index_L];   w_R     = w[i][j][index_R]
                E_L     = E[i][j][index_L];   E_R     = E[i][j][index_R]
                P_L     = P[i][j][index_L];   P_R     = P[i][j][index_R]
                a_L     = sos[i][j][index_L]; a_R     = sos[i][j][index_R]
                # rho
                var_type = 0
                rho_F_m  = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
                # rhou
                var_type = 3
                rhou_F_m = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )               
                # rhov
                var_type = 2
                rhov_F_m = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )         
                # rhow
                var_type = 1
                rhow_F_m = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )     
                # rhoE
                var_type = 4
                rhoE_F_m = KGP_flux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type )
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
def viscous_fluxes( rhou_vis, rhov_vis, rhow_vis, rhoE_vis, u, v, w, T, mu, kappa, grid ):

    # Second-order central finite differences for derivatives:
    # P. Moin.
    # Fundamentals of engineering numerical analysis.
    # Cambridge University Press, 2010.

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                ## Geometric stuff
                delta_x = 0.5*( grid[i+1][j][k][0] - grid[i-1][j][k][0] ) 
                delta_y = 0.5*( grid[i][j+1][k][1] - grid[i][j-1][k][1] ) 
                delta_z = 0.5*( grid[i][j][k+1][2] - grid[i][j][k-1][2] )
                ## Velocity derivatives
                d_u_x = ( u[i+1][j][k] - u[i-1][j][k] )/delta_x
                d_u_y = ( u[i][j+1][k] - u[i][j-1][k] )/delta_y
                d_u_z = ( u[i][j][k+1] - u[i][j][k-1] )/delta_z
                d_v_x = ( v[i+1][j][k] - v[i-1][j][k] )/delta_x
                d_v_y = ( v[i][j+1][k] - v[i][j-1][k] )/delta_y
                d_v_z = ( v[i][j][k+1] - v[i][j][k-1] )/delta_z
                d_w_x = ( w[i+1][j][k] - w[i-1][j][k] )/delta_x
                d_w_y = ( w[i][j+1][k] - w[i][j-1][k] )/delta_y
                d_w_z = ( w[i][j][k+1] - w[i][j][k-1] )/delta_z
                ## Temperature derivatives
                d_T_x = ( T[i+1][j][k] - T[i-1][j][k] )/delta_x;
                d_T_y = ( T[i][j+1][k] - T[i][j-1][k] )/delta_y;
                d_T_z = ( T[i][j][k+1] - T[i][j][k-1] )/delta_z;
                ## Transport coefficients derivatives
                d_mu_x    = ( mu[i+1][j][k] - mu[i-1][j][k] )/delta_x;
                d_mu_y    = ( mu[i][j+1][k] - mu[i][j-1][k] )/delta_y;
                d_mu_z    = ( mu[i][j][k+1] - mu[i][j][k-1] )/delta_z;
                d_kappa_x = ( kappa[i+1][j][k] - kappa[i-1][j][k] )/delta_x;
                d_kappa_y = ( kappa[i][j+1][k] - kappa[i][j-1][k] )/delta_y;
                d_kappa_z = ( kappa[i][j][k+1] - kappa[i][j][k-1] )/delta_z;
                ## Divergence of velocity
                div_uvw = d_u_x + d_v_y + d_w_z
                ## Viscous stresses ( symmetric tensor )
                tau_xx = 2.0*mu[i][j][k]*( d_u_x - ( div_uvw/3.0 ) )
                tau_xy = mu[i][j][k]*( d_u_y + d_v_x )
                tau_xz = mu[i][j][k]*( d_u_z + d_w_x )
                tau_yy = 2.0*mu[i][j][k]*( d_v_y - ( div_uvw/3.0 ) )
                tau_yz = mu[i][j][k]*( d_v_z + d_w_y )
                tau_zz = 2.0*mu[i][j][k]*( d_w_z - ( div_uvw/3.0 ) )
                ## Divergence of viscous stresses
                div_tau_x = mu[i][j][k]*( ( 1.00/delta_x )*( ( u[i+1][j][k] - u[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] )
                                                           - ( u[i][j][k] - u[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) )
                                        + ( 1.00/delta_y )*( ( u[i][j+1][k] - u[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] )
                                                           - ( u[i][j][k] - u[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) )
                                        + ( 1.00/delta_z )*( ( u[i][j][k+1] - u[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] )
                                                           - ( u[i][j][k] - u[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) )			   + ( 1.0/3.0 )*mu[i][j][k]*( ( 1.00/delta_x )*( ( u[i+1][j][k] - u[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] )
                                                           - ( u[i][j][k] - u[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) )
                                        + ( 0.25/delta_x )*( ( v[i+1][j+1][k] - v[i+1][j-1][k] )/delta_y
                                                           - ( v[i-1][j+1][k] - v[i-1][j-1][k] )/delta_y )
                                        + ( 0.25/delta_x )*( ( w[i+1][j][k+1] - w[i+1][j][k-1] )/delta_z
                                                           - ( w[i-1][j][k+1] - w[i-1][j][k-1] )/delta_z ) )                                                                   + d_mu_x*tau_xx + d_mu_y*tau_xy + d_mu_z*tau_xz
                div_tau_y = mu[i][j][k]*( ( 1.00/delta_x )*( ( v[i+1][j][k] - v[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] )
                                                           - ( v[i][j][k] - v[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) )
                                        + ( 1.00/delta_y )*( ( v[i][j+1][k] - v[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] )
                                                           - ( v[i][j][k] - v[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) )
                                        + ( 1.00/delta_z )*( ( v[i][j][k+1] - v[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] )
                                                           - ( v[i][j][k] - v[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) ) 		           + ( 1.0/3.0 )*mu[i][j][k]*( ( 0.25/delta_y )*( ( u[i+1][j+1][k] - u[i-1][j+1][k] )/delta_x
                                                           - ( u[i+1][j-1][k] - u[i-1][j-1][k] )/delta_x )
                                        + ( 1.00/delta_y )*( ( v[i][j+1][k] - v[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] )
                                                           - ( v[i][j][k] - v[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) )
                                        + ( 0.25/delta_y )*( ( w[i][j+1][k+1] - w[i][j+1][k-1] )/delta_z
                                                           - ( w[i][j-1][k+1] - w[i][j-1][k-1] )/delta_z ) )                                                                   + d_mu_x*tau_xy + d_mu_y*tau_yy + d_mu_z*tau_yz
                div_tau_z = mu[i][j][k]*( ( 1.00/delta_x )*( ( w[i+1][j][k] - w[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] )
                                                           - ( w[i][j][k] - w[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) )
                                        + ( 1.00/delta_y )*( ( w[i][j+1][k] - w[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] )
                                                           - ( w[i][j][k] - w[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) )
                                        + ( 1.00/delta_z )*( ( w[i][j][k+1] - w[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] )
                                                           - ( w[i][j][k] - w[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) )			   + ( 1.0/3.0 )*mu[i][j][k]*( ( 0.25/delta_z )*( ( u[i+1][j][k+1] - u[i-1][j][k+1] )/delta_x
                                                           - ( u[i+1][j][k-1] - u[i-1][j][k-1] )/delta_x )
                                        + ( 0.25/delta_z )*( ( v[i][j+1][k+1] - v[i][j-1][k+1] )/delta_y
                                                           - ( v[i][j+1][k-1] - v[i][j-1][k-1] )/delta_y )
                                        + ( 1.00/delta_z )*( ( w[i][j][k+1] - w[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] )
                                                           - ( w[i][j][k] - w[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) )                                       + d_mu_x*tau_xz + d_mu_y*tau_yz + d_mu_z*tau_zz
                ## Fourier term
                div_q = ( -1.0 )*kappa[i][j][k]*( ( 1.0/delta_x )*( ( T[i+1][j][k] - T[i][j][k] )/( grid[i+1][j][k][0] - grid[i][j][k][0] )
                                                                  - ( T[i][j][k] - T[i-1][j][k] )/( grid[i][j][k][0] - grid[i-1][j][k][0] ) )
                                                + ( 1.0/delta_y )*( ( T[i][j+1][k] - T[i][j][k] )/( grid[i][j+1][k][1] - grid[i][j][k][1] )
                                                                  - ( T[i][j][k] - T[i][j-1][k] )/( grid[i][j][k][1] - grid[i][j-1][k][1] ) )
                                                + ( 1.0/delta_z )*( ( T[i][j][k+1] - T[i][j][k] )/( grid[i][j][k+1][2] - grid[i][j][k][2] )
                                                                  - ( T[i][j][k] - T[i][j][k-1] )/( grid[i][j][k][2] - grid[i][j][k-1][2] ) ) )                            - d_kappa_x*d_T_x - d_kappa_y*d_T_y - d_kappa_z*d_T_z
                ## Work of viscous stresses
                div_uvw_tau = u[i][j][k]*div_tau_x + v[i][j][k]*div_tau_y + w[i][j][k]*div_tau_z								                 + tau_xx*d_u_x + tau_xy*d_u_y + tau_xz*d_u_z													      + tau_xy*d_v_x + tau_yy*d_v_y + tau_yz*d_v_z													   + tau_xz*d_w_x + tau_yz*d_w_y + tau_zz*d_w_z
                ## Viscous fluxes
                rhou_vis[i][j][k] = div_tau_x
                rhov_vis[i][j][k] = div_tau_y
                rhow_vis[i][j][k] = div_tau_z
                rhoE_vis[i][j][k] = ( -1.0 )*div_q + div_uvw_tau
    #print( rhou_vis )
    #print( rhov_vis )
    #print( rhow_vis )
    #print( rhoE_vis )


### Sum inviscid & viscous fluxes and source terms
def sum_fluxes_source_terms( rho_tot, rhou_tot, rhov_tot, rhow_tot, rhoE_tot, rho_inv, rhou_inv, rhov_inv, rhow_inv, rhoE_inv, rhou_vis, rhov_vis, rhow_vis, rhoE_vis, f_rhou, f_rhov, f_rhow, f_rhoE, rk_iter ):

    # Internal points
    for i in range( 1, num_grid_x + 1 ):    
        for j in range( 1, num_grid_y + 1 ):    
            for k in range( 1, num_grid_z + 1 ):
                rho_tot[i][j][k][rk_iter]  = ( -1.0 )*rho_inv[i][j][k]
                rhou_tot[i][j][k][rk_iter] = ( -1.0 )*rhou_inv[i][j][k] + rhou_vis[i][j][k] + f_rhou[i][j][k]
                rhov_tot[i][j][k][rk_iter] = ( -1.0 )*rhov_inv[i][j][k] + rhov_vis[i][j][k] + f_rhov[i][j][k]
                rhow_tot[i][j][k][rk_iter] = ( -1.0 )*rhow_inv[i][j][k] + rhow_vis[i][j][k] + f_rhow[i][j][k]
                rhoE_tot[i][j][k][rk_iter] = ( -1.0 )*rhoE_inv[i][j][k] + rhoE_vis[i][j][k] + f_rhoE[i][j][k]
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



########## MAIN ##########

### START SIMULATION
print( 'RHEA: START SIMULATION' )

### Define spatial discretization
spatial_discretization( mesh )

### Initialize u, v, w, P and T variables
initialize_uvwPT( u_field, v_field, w_field, P_field, T_field, mesh )

### Initialize thermodynamic variables
initialize_thermodynamics( rho_field, E_field, sos_field, u_field, v_field, w_field, P_field, T_field )

### Calculate transport coefficients
transport_coefficients( mu_field, kappa_field )

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
    delta_t = time_step( rho_field, u_field, v_field, w_field, sos_field, mu_field, kappa_field, mesh )
    if( ( time + delta_t ) > final_time ):
        delta_t = final_time - time
        #print( delta_t )

    ### Print time iteration iformation
    print( 'Time iteration ' + str( t ) + ': t = ' + str( time ) + ' [s], delta_t = ' + str( delta_t ) + ' [s]' )

    ### Output data to file
    if( t % output_iter == 0 ):
        data_output( rho_field, u_field, v_field, w_field, E_field, P_field, T_field, sos_field, mesh )

    ### Runge-Kutta sub-steps
    for rk in range( 0, rk_order ):

        ### Calculate transport coefficients
        transport_coefficients( mu_field, kappa_field )

        ### Calculate inviscid fluxes
        inviscid_fluxes( rho_inv_flux, rhou_inv_flux, rhov_inv_flux, rhow_inv_flux, rhoE_inv_flux, rho_field, u_field, v_field, w_field, E_field, P_field, sos_field, mesh )

        ### Calculate viscous fluxes
        viscous_fluxes( rhou_vis_flux, rhov_vis_flux, rhow_vis_flux, rhoE_vis_flux, u_field, v_field, w_field, T_field, mu_field, kappa_field, mesh )

        ### Calculate source terms
        source_terms( f_rhou_field, f_rhov_field, f_rhow_field, f_rhoE_field, rho_field, u_field, v_field, w_field, mesh )

        ### Sum fluxes & source terms
        sum_fluxes_source_terms( rho_rk_fluxes, rhou_rk_fluxes, rhov_rk_fluxes, rhow_rk_fluxes, rhoE_rk_fluxes, rho_inv_flux, rhou_inv_flux, rhov_inv_flux, rhow_inv_flux, rhoE_inv_flux, rhou_vis_flux, rhov_vis_flux, rhow_vis_flux, rhoE_vis_flux, f_rhou_field, f_rhov_field, f_rhow_field, f_rhoE_field, rk )

        ### Advance conserved variables in time
        time_integration( rho_field,  rho_0_field,  delta_t, rho_rk_fluxes,  rk )
        time_integration( rhou_field, rhou_0_field, delta_t, rhou_rk_fluxes, rk )
        time_integration( rhov_field, rhov_0_field, delta_t, rhov_rk_fluxes, rk )
        time_integration( rhow_field, rhow_0_field, delta_t, rhow_rk_fluxes, rk )
        time_integration( rhoE_field, rhoE_0_field, delta_t, rhoE_rk_fluxes, rk )

        ### Update primitive variables from conserved variables
        update_primitive( u_field, rhou_field, rho_field )
        update_primitive( v_field, rhov_field, rho_field )
        update_primitive( w_field, rhow_field, rho_field )
        update_primitive( E_field, rhoE_field, rho_field )
       
        ### Update thermodynamic variables from primitive variables
        thermodynamic_state( P_field, T_field, sos_field, rho_field, u_field, v_field, w_field, E_field )

        ### Update boundaries
        update_boundaries( rho_field, rhou_field, rhov_field, rhow_field, rhoE_field, u_field, v_field, w_field, P_field, T_field, mesh )

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

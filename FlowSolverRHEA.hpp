//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                              //
// RHEA - an open-source Reproducible Hybrid-architecture flow solver Engineered for Academia                   //
//                                                                                                              //
// Rhea was the Titaness great Mother of the Gods, and goddess of female fertility, motherhood, and generation. //
// Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.           //
//                                                                                                              //
//                                                                                                              //
// REHA is released under the MIT License:                                                                      //
//                                                                                                              //
// Copyright (c) 2020 Lluis Jofre Cruanyes & Guillermo Oyarzun Altamirano.                                      //
//                                                                                                              //
// Permission is hereby granted, free of charge, to any person obtaining a copy                                 //
// of this software and associated documentation files (the "Software"), to deal                                //
// in the Software without restriction, including without limitation the rights                                 //
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell                                    //
// copies of the Software, and to permit persons to whom the Software is                                        //
// furnished to do so, subject to the following conditions:                                                     //
//                                                                                                              //
// The above copyright notice and this permission notice shall be included in all                               //
// copies or substantial portions of the Software.                                                              //
//                                                                                                              //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR                                   //
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                                     //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE                                  //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                                       //
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,                                //
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE                                //
// SOFTWARE.                                                                                                    //
//                                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _RHEA_FLOW_SOLVER_
#define _RHEA_FLOW_SOLVER_

////////// INCLUDES //////////
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <mpi.h>
#include "src/parameters.h"
#include "src/domain.h"
#include "src/comm_scheme.h"
#include "src/parvec.h"

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class FlowSolverRHEA;						/// Flow solver RHEA

////////// FUNCTION DECLARATION //////////


////////// FlowSolverRHEA CLASS //////////
class FlowSolverRHEA {
   
    ////////// VARIABLES & PARAMETERS DESCRIPTION //////////

    /// Primitive variables:
    ///   - Density rho
    ///   - Velocities u, v, w
    ///   - Specific total energy E = e + ke
    ///     ... sum of internal energy e and kinetic energy ke = (u*u + v*v + w*w)/2

    /// Conserved variables:
    ///   - Mass rho
    ///   - Momentum rho*u, rho*v, rho*w
    ///   - Total energy rho*E

    /// Thermodynamic state:
    ///   - Pressure P
    ///   - Temperature T
    ///   - Speed of sound sos

    /// Thermophysical properties:
    ///   - Specific gas constant R_specific
    ///   - Ratio of heat capacities gamma
    ///   - Dynamic viscosity mu
    ///   - Thermal conductivity kappa
 
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        FlowSolverRHEA();					/// Default constructor
        FlowSolverRHEA(const string configuration_file);	/// Parametrized constructor
        FlowSolverRHEA(const FlowSolverRHEA &in);		/// Copy constructor
        virtual ~FlowSolverRHEA();				/// Destructor

    protected:

        ////////// SOLVER PARAMETERS //////////
	
        /// Fluid properties 
        double R_specific;					/// Specific gas constant [J/(kg K)]
        double gamma;						/// Heat capacity ratio [-]
        double mu;						/// Dynamic viscosity [Pa s]
        double kappa;						/// Thermal conductivity [W/(m k)]       

        /// Problem parameters
        double x_0;           	        	 	 	/// Domain origin in x-direction [m]
        double y_0;             	        		/// Domain origin in y-direction [m]
        double z_0;                     			/// Domain origin in z-direction [m]
        double L_x;   						/// Domain size in x-direction
        double L_y;      					/// Domain size in y-direction
        double L_z;						/// Domain size in z-direction
        double initial_time;   					/// Initial time [s]
        double final_time;		      			/// Final time [s]
        string configuration_file;				/// Configuration file (YAML language)	

        /// Computational parameters
        double num_grid_x;					/// Number of internal grid points in the x-direction
        double num_grid_y;					/// Number of internal grid points in the y-direction
        double num_grid_z;					/// Number of internal grid points in the z-direction
        double CFL;						/// CFL coefficient
        double max_num_time_iter;				/// Maximum number of time iterations
        double output_iter;					/// Output data every given number of iterations

        /// Boundary conditions
        int bocos[6];						/// Array of boundary conditions

        /// Parallelization scheme
        int np_x;						/// Number of processes in x-direction
        int np_y;						/// Number of processes in y-direction
        int np_z;						/// Number of processes in z-direction

	////////// SOLVER (PARALLELIZED) VARIABLES //////////
	
        /// Primitive, conserved and thermodynamic variables
        parvec *rho_field;					/// 3-D field of rho
        parvec *u_field;					/// 3-D field of u
        parvec *v_field;					/// 3-D field of v
        parvec *w_field;					/// 3-D field of w
        parvec *E_field;					/// 3-D field of E
        parvec *rhou_field;					/// 3-D field of rhou
        parvec *rhov_field;					/// 3-D field of rhov
        parvec *rhow_field;					/// 3-D field of rhow
        parvec *rhoE_field;					/// 3-D field of rhoE
        parvec *P_field;					/// 3-D field of P
        parvec *T_field;					/// 3-D field of T
        parvec *sos_field;					/// 3-D field of sos

        /// Time-integration variables
        parvec *rho_0_field;					/// 3-D old field of rho
        parvec *rhou_0_field;					/// 3-D old field of rhou
        parvec *rhov_0_field;					/// 3-D old field of rhov
        parvec *rhow_0_field;					/// 3-D old field of rhow
        parvec *rhoE_0_field;					/// 3-D old field of rhoE

        /// Time-integration fluxes
        parvec *rho_rk1_flux;					/// 3-D Runge-Kutta flux 1 of rho
        parvec *rho_rk2_flux;					/// 3-D Runge-Kutta flux 2 of rho
        parvec *rho_rk3_flux;					/// 3-D Runge-Kutta flux 3 of rho
        parvec *rhou_rk1_flux;					/// 3-D Runge-Kutta flux 1 of rhou
        parvec *rhou_rk2_flux;					/// 3-D Runge-Kutta flux 2 of rhou
        parvec *rhou_rk3_flux;					/// 3-D Runge-Kutta flux 3 of rhou
        parvec *rhov_rk1_flux;					/// 3-D Runge-Kutta flux 1 of rhov
        parvec *rhov_rk2_flux;					/// 3-D Runge-Kutta flux 2 of rhov
        parvec *rhov_rk3_flux;					/// 3-D Runge-Kutta flux 3 of rhov
        parvec *rhow_rk1_flux;					/// 3-D Runge-Kutta flux 1 of rhow
        parvec *rhow_rk2_flux;					/// 3-D Runge-Kutta flux 2 of rhow
        parvec *rhow_rk3_flux;					/// 3-D Runge-Kutta flux 3 of rhow
        parvec *rhoE_rk1_flux;					/// 3-D Runge-Kutta flux 1 of rhoE
        parvec *rhoE_rk2_flux;					/// 3-D Runge-Kutta flux 2 of rhoE
        parvec *rhoE_rk3_flux;					/// 3-D Runge-Kutta flux 3 of rhoE

        /// Inviscid fluxes
        parvec *rho_inv_flux;					/// 3-D inviscid fluxes of rho
        parvec *rhou_inv_flux;					/// 3-D inviscid fluxes of rhou
        parvec *rhov_inv_flux;					/// 3-D inviscid fluxes of rhov
        parvec *rhow_inv_flux;					/// 3-D inviscid fluxes of rhow
        parvec *rhoE_inv_flux;					/// 3-D inviscid fluxes of rhoE

        /// Viscous fluxes
        parvec *rhou_vis_flux;					/// 3-D viscous fluxes of rhou
        parvec *rhov_vis_flux;					/// 3-D viscous fluxes of rhov
        parvec *rhow_vis_flux;					/// 3-D viscous fluxes of rhow
        parvec *rhoE_vis_flux;					/// 3-D viscous fluxes of rhoE

        /// Volumetric & external forces
        parvec *f_rhou_field;					/// 3-D field of rhou
        parvec *f_rhov_field;					/// 3-D field of rhov
        parvec *f_rhow_field;					/// 3-D field of rhow
        parvec *f_rhoE_field;					/// 3-D field of rhoE

	////////// COMPUTATIONAL DOMAIN & PARALLEL TOPOLOGY //////////
        domain *dom;						/// Computational domain
        comm_scheme *topo;					/// Communication scheme (parallel topology)

	////////// SOLVER METHODS //////////
	

    private:

};

#endif /*_RHEA_FLOW_SOLVER_*/

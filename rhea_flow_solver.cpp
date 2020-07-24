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

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <mpi.h>
#include "src/parameters.h"
#include "src/domain.h"
#include "src/comm_scheme.h"
#include "src/parvec.h"

using namespace std;


// Primitive variables:
//   - Density rho
//   - Velocities u, v, w
//   - Specific total energy E
//     ... E = e + ke
//     ... is the sum of internal energy e
//     ... and kinetic energy ke = (u*u + v*v + w*w)/2
//
// Conserved variables:
//   - Mass rho
//   - Momentum rho*u, rho*v, rho*w
//   - Total energy rho*E
//
// Thermodynamic state:
//   - Pressure P
//   - Temperature T
//   - Speed of sound sos
//
// Thermophysical properties:
//   - Specific gas constant R_specific
//   - Ratio of heat capacities gamma
//   - Dynamic viscosity mu
//   - Thermal conductivity kappa

// Fixed parameters
const double rk_order = 3;			// Time-integration Runge-Kutta order (fixed value)
const double epsilon  = 1.0e-15;		// Small epsilon number (fixed value)
const double Pi       = 2.0*asin(1.0);		// Pi number

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);

    // The lines below need to be introduced via configuration (input) file
    const int _DEBUG_ = 0;

    ////////// SET PARAMETERS //////////
    const double Re_tau  = 100.0;				// Friction Reynolds number [-]
    const double delta   = 1.0;					// Channel half-height [m]
    const double u_tau   = 1.0;					// Friction velocity [m/s]
    const double rho_ref = 1.0;					// Reference density [kg/m3]
    const double P_ref   = 101325.0;				// Reference pressure [Pa]

    /// Fluid properties 
    const double R_specific = 287.058;				// Specific gas constant of air [J/(kg K)]
    const double gamma      = 1.4;				// Heat capacity ratio of air [-]
    const double mu         = rho_ref*u_tau*delta/Re_tau;	// Dynamic viscosity of air [Pa s]
    const double kappa      = 0.0;				// Thermal conductivity of air [W/(m k)]

    /// Problem parameters
    const double x_0           = 0.0;                     	// Domain origin in x-direction [m]
    const double y_0           = 0.0;                     	// Domain origin in y-direction [m]
    const double z_0           = 0.0;                     	// Domain origin in z-direction [m]
    const double L_x           = 4.0*Pi*delta;   		// Size of domain in x-direction
    const double L_y           = 2.0*delta;      		// Size of domain in y-direction
    const double L_z           = 4.0*Pi*delta/3.0;		// Size of domain in z-direction
    const double initial_time  = 0.0;   			// Initial time [s]
    const double final_time    = 1.0e3;      			// Final time [s]
    const string name_file_out = "output_data.csv";		// Name of output data [-]

    /// Computational parameters
    const double num_grid_x        = 64;			// Number of internal grid points in the x-direction
    const double num_grid_y        = 64;			// Number of internal grid points in the y-direction
    const double num_grid_z        = 64;			// Number of internal grid points in the z-direction
    const double CFL               = 0.9;			// CFL coefficient
    const double max_num_time_iter = 1e6;			// Maximum number of time iterations
    const double output_iter       = 1e2;			// Output data every given number of iterations

    /// Boundary conditions
    int bocos[6];
    bocos[_WEST_]  = _PERIODIC_;
    bocos[_EAST_]  = _PERIODIC_;
    bocos[_SOUTH_] = _PERIODIC_;
    bocos[_NORTH_] = _PERIODIC_;
    bocos[_BACK_]  = _PERIODIC_;
    bocos[_FRONT_] = _PERIODIC_;

    /// Parallelization scheme
    const int npx = 2;
    const int npy = 1;
    const int npz = 1;
    // The lines above need to be introduced via configuration (input) file

    //This is just for testing, later would be read from a file
    domain dom(L_x, L_y, L_z, x_0, y_0, z_0, RHEA_NX, RHEA_NY, RHEA_NZ);

    //To add the boundary conditions,  later would be read from a file
    dom.updateBocos(bocos);

    comm_scheme topo(&dom, npx, npy, npz);

    if(topo.getRank() == 0)
        dom.printDomain();

   for(int p = 0; p < npx*npy*npz; p++)
       topo.printCommSchemeToFile(p);


    ////////// ALLOCATE MEMORY //////////

    /// Primitive, conserved and thermodynamic variables
    parvec rho_field(&topo);		// 3-D field of rho
    parvec u_field(&topo);		// 3-D field of u
    parvec v_field(&topo);		// 3-D field of v
    parvec w_field(&topo);		// 3-D field of w
    parvec E_field(&topo);		// 3-D field of E
    parvec rhou_field(&topo);		// 3-D field of rhou
    parvec rhov_field(&topo);		// 3-D field of rhov
    parvec rhow_field(&topo);		// 3-D field of rhow
    parvec rhoE_field(&topo);		// 3-D field of rhoE
    parvec P_field(&topo);		// 3-D field of P
    parvec T_field(&topo);		// 3-D field of T
    parvec sos_field(&topo);		// 3-D field of sos

    /// Time integration variables
    parvec rho_0_field(&topo);		// 3-D old field of rho
    parvec rhou_0_field(&topo);		// 3-D old field of rhou
    parvec rhov_0_field(&topo);		// 3-D old field of rhov
    parvec rhow_0_field(&topo);		// 3-D old field of rhow
    parvec rhoE_0_field(&topo);		// 3-D old field of rhoE

    /// Time integration fluxes
    parvec rho_rk1_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rho
    parvec rho_rk2_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rho
    parvec rho_rk3_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rho
    parvec rhou_rk1_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhou
    parvec rhou_rk2_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhou
    parvec rhou_rk3_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhou
    parvec rhov_rk1_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhov
    parvec rhov_rk2_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhov
    parvec rhov_rk3_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhov
    parvec rhow_rk1_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhow
    parvec rhow_rk2_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhow
    parvec rhow_rk3_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhow
    parvec rhoE_rk1_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhoE
    parvec rhoE_rk2_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhoE
    parvec rhoE_rk3_fluxes(&topo);	// 3-D Runge-Kutta fluxes of rhoE

    /// Inviscid fluxes
    parvec rho_inv_flux(&topo);		// 3-D inviscid fluxes of rho
    parvec rhou_inv_flux(&topo);	// 3-D inviscid fluxes of rhou
    parvec rhov_inv_flux(&topo);	// 3-D inviscid fluxes of rhov
    parvec rhow_inv_flux(&topo);	// 3-D inviscid fluxes of rhow
    parvec rhoE_inv_flux(&topo);	// 3-D inviscid fluxes of rhoE

    /// Viscous fluxes
    parvec rhou_vis_flux(&topo);	// 3-D viscous fluxes of rhou
    parvec rhov_vis_flux(&topo);	// 3-D viscous fluxes of rhov
    parvec rhow_vis_flux(&topo);	// 3-D viscous fluxes of rhow
    parvec rhoE_vis_flux(&topo);	// 3-D viscous fluxes of rhoE

    /// Body forces
    parvec f_rhou_field(&topo);		// 3-D field of rhou
    parvec f_rhov_field(&topo);		// 3-D field of rhov
    parvec f_rhow_field(&topo);		// 3-D field of rhow
    parvec f_rhoE_field(&topo);		// 3-D field of rhoE


/*    parvec T(&topo);
    parvec Tnew(&topo);

    T = 0.0;
    Tnew = 0.0;

    int maxiter= 20000;


    for(int it=0; it < maxiter; it++)
    {

        for(int i = Tnew.ini_x; i <=  Tnew.fin_x; i++)
            for(int j = Tnew.ini_y; j <=  Tnew.fin_y; j++)
                for(int k = Tnew.ini_z; k <=  Tnew.fin_z; k++){
                    Tnew[I1D(i,j,k)] =  T[I1D(i,j,k)] + dt*kappa*( (2.0/(dom.x[i+1]-dom.x[i-1]))*( (T[I1D(i+1,j,k)] -T[I1D(i,j,k)])/(dom.x[i+1]-dom.x[i]) -  (T[I1D(i,j,k)] -T[I1D(i-1,j,k)])/(dom.x[i]-dom.x[i-1])) + 
                                                                   (2.0/(dom.y[j+1]-dom.y[j-1]))*( (T[I1D(i,j+1,k)] -T[I1D(i,j,k)])/(dom.y[j+1]-dom.y[j]) -  (T[I1D(i,j,k)] -T[I1D(i,j-1,k)])/(dom.y[j]-dom.y[j-1])) +        
                                                                   (2.0/(dom.z[k+1]-dom.z[k-1]))*( (T[I1D(i,j,k+1)] -T[I1D(i,j,k)])/(dom.z[k+1]-dom.z[k]) -  (T[I1D(i,j,k)] -T[I1D(i,j,k-1)])/(dom.z[k]-dom.z[k-1])) );    
                }
    
        Tnew.update();
        //WEST 
        for(int i = topo.iter_bound[_WEST_][_INIX_]; i <= topo.iter_bound[_WEST_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_WEST_][_INIY_]; j <= topo.iter_bound[_WEST_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_WEST_][_INIZ_]; k <=  topo.iter_bound[_WEST_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = 2*0.0 - Tnew[I1D(i+1,j,k)];
                }

        //EAST 
        for(int i = topo.iter_bound[_EAST_][_INIX_]; i <= topo.iter_bound[_EAST_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_EAST_][_INIY_]; j <= topo.iter_bound[_EAST_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_EAST_][_INIZ_]; k <=  topo.iter_bound[_EAST_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = 2*1.0 - Tnew[I1D(i-1,j,k)];
                }

        //SOUTH 
        for(int i = topo.iter_bound[_SOUTH_][_INIX_]; i <= topo.iter_bound[_SOUTH_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_SOUTH_][_INIY_]; j <= topo.iter_bound[_SOUTH_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_SOUTH_][_INIZ_]; k <=  topo.iter_bound[_SOUTH_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,j+1,k)];
                }

        //NORTH 
        for(int i = topo.iter_bound[_NORTH_][_INIX_]; i <= topo.iter_bound[_NORTH_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_NORTH_][_INIY_]; j <= topo.iter_bound[_NORTH_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_NORTH_][_INIZ_]; k <=  topo.iter_bound[_NORTH_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,(j-1),k)];
                }

        //BACK 
        for(int i = topo.iter_bound[_BACK_][_INIX_]; i <= topo.iter_bound[_BACK_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_BACK_][_INIY_]; j <= topo.iter_bound[_BACK_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_BACK_][_INIZ_]; k <=  topo.iter_bound[_BACK_][_ENDZ_]; k++)
                {
//                    Tnew[I1D(i,j,0)] = Tnew[I1D(i,j,_lNz_)];
                     Tnew[I1D(i,j,k)] = Tnew[I1D(i,j,k+1)];
                }

        //FRONT 
        for(int i = topo.iter_bound[_FRONT_][_INIX_]; i <= topo.iter_bound[_FRONT_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_FRONT_][_INIY_]; j <= topo.iter_bound[_FRONT_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_FRONT_][_INIZ_]; k <=  topo.iter_bound[_FRONT_][_ENDZ_]; k++)
                {
//                   Tnew[I1D(i,j,_lNz_+1)] = Tnew[I1D(i,j,1)];
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,j,k-1)];
               }



//        for(int i = topo.iter_common[_ALL_][_INIX_]; i <= topo.iter_common[_ALL_][_INIX_]; i++)
//            for(int j =  topo.iter_common[_ALL_][_INIY_]; j <=  topo.iter_common[_ALL_][_INIY_]; j++)
//                for(int k =  topo.iter_common[_ALL_][_INIZ_]; k <=  topo.iter_common[_ALL_][_INIZ_]; k++){
//                    T[I1D(i,j,k)] = Tnew[I1D(i,j,k)];
//                }
 
        for(int l=0; l < topo.getSize(); l++)
        {
            T[l] = Tnew[l];
        }
 
    }



    if(topo.getRank() == 0 && _DEBUG_ ==0 )
        for(int l=0; l < topo.getSize(); l++)
            cout<<l<<":  "<<Tnew[l]<<endl;
*/


    MPI_Finalize();
}

#include "myRHEA.hpp"

using namespace std;

/// pi number
//const double pi = 2.0*asin( 1.0 );

/// PROBLEM PARAMETERS ///
const double gamma_0 = 1.4;					/// Reference ratio of heat capacities
//const double Re_0    = pi;					/// Reynolds number
const double Ma_0    = 1.0e-2/sqrt( gamma_0 );			/// Mach number
const double rho_0   = 1.0;					/// Reference density	
const double U_0     = 1.0;					/// Reference velocity
//const double mu_0    = rho_0*U_0*pi/Re_0;			/// Dynamic viscosity	
const double P_0     = rho_0*U_0*U_0/( gamma_0*Ma_0*Ma_0 );	/// Reference pressure
//const double L       = 2.0*pi;					/// Domain size			


////////// myRHEA CLASS //////////

void myRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// All (inner, boundary & halo) points: u, v, w, P and T
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                u_field[I1D(i,j,k)] = ( -1.0 )*U_0*cos( mesh->x[i] )*sin( mesh->y[j] );
                v_field[I1D(i,j,k)] = U_0*sin( mesh->x[i] )*cos( mesh->y[j] );
                w_field[I1D(i,j,k)] = 0.0;
                P_field[I1D(i,j,k)] = P_0 - ( rho_0*U_0*U_0/4.0 )*( cos( 2.0*mesh->x[i] ) + cos( 2.0*mesh->y[j] ) );
                T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( rho_0*thermodynamics->getSpecificGasConstant() );
            }
        }
    }

    /// Update halo values
    u_field.update();
    v_field.update();
    w_field.update();
    P_field.update();
    T_field.update();

};

void myRHEA::calculateSourceTerms() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                f_rhou_field[I1D(i,j,k)] = 0.0;
                f_rhov_field[I1D(i,j,k)] = 0.0;
                f_rhow_field[I1D(i,j,k)] = 0.0;
                f_rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

};

void myRHEA::calculateInviscidFluxes() {

    /// Inner points: rho, rhou, rhov, rhow and rhoE
    int index_L, index_R;
    double delta_x, delta_y, delta_z;
    double rho_L, u_L, v_L, w_L, E_L, P_L;
    double rho_R, u_R, v_R, w_R, E_R, P_R;
    double rho_F_p, rho_F_m;
    double rhou_F_p, rhou_F_m;
    double rhov_F_p, rhov_F_m;
    double rhow_F_p, rhow_F_m;
    double rhoE_F_p, rhoE_F_m;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = 0.5*( mesh->x[i+1] - mesh->x[i-1] ); 
                delta_y = 0.5*( mesh->y[j+1] - mesh->y[j-1] ); 
                delta_z = 0.5*( mesh->z[k+1] - mesh->z[k-1] );
                /// x-direction i+1/2
                index_L = i;                           index_R = i + 1;
                rho_L   = rho_field[I1D(index_L,j,k)]; rho_R   = rho_field[I1D(index_R,j,k)]; 
                u_L     = u_field[I1D(index_L,j,k)];   u_R     = u_field[I1D(index_R,j,k)];
                v_L     = v_field[I1D(index_L,j,k)];   v_R     = v_field[I1D(index_R,j,k)];
                w_L     = w_field[I1D(index_L,j,k)];   w_R     = w_field[I1D(index_R,j,k)];
                E_L     = E_field[I1D(index_L,j,k)];   E_R     = E_field[I1D(index_R,j,k)];
                P_L     = P_field[I1D(index_L,j,k)];   P_R     = P_field[I1D(index_R,j,k)];		
                /// rho, rhou, rhov, rhow, rhoE
                rho_F_p  = 0.5*( rho_L*u_L + rho_R*u_R );
                rhou_F_p = 0.5*( rho_L*u_L*u_L + P_L + rho_R*u_R*u_R + P_R );
                rhov_F_p = 0.5*( rho_L*u_L*v_L + rho_R*u_R*v_R );
                rhow_F_p = 0.5*( rho_L*u_L*w_L + rho_R*u_R*w_R );
                rhoE_F_p = 0.5*( rho_L*u_L*E_L + u_L*P_L + rho_R*u_R*E_R + u_R*P_R );
                /// x-direction i-1/2
                index_L = i - 1;                       index_R = i;
                rho_L   = rho_field[I1D(index_L,j,k)]; rho_R   = rho_field[I1D(index_R,j,k)];
                u_L     = u_field[I1D(index_L,j,k)];   u_R     = u_field[I1D(index_R,j,k)];
                v_L     = v_field[I1D(index_L,j,k)];   v_R     = v_field[I1D(index_R,j,k)];
                w_L     = w_field[I1D(index_L,j,k)];   w_R     = w_field[I1D(index_R,j,k)];
                E_L     = E_field[I1D(index_L,j,k)];   E_R     = E_field[I1D(index_R,j,k)];
                P_L     = P_field[I1D(index_L,j,k)];   P_R     = P_field[I1D(index_R,j,k)];
                /// rho, rhou, rhov, rhow, rhoE
                rho_F_m  = 0.5*( rho_L*u_L + rho_R*u_R );
                rhou_F_m = 0.5*( rho_L*u_L*u_L + P_L + rho_R*u_R*u_R + P_R );
                rhov_F_m = 0.5*( rho_L*u_L*v_L + rho_R*u_R*v_R );
                rhow_F_m = 0.5*( rho_L*u_L*w_L + rho_R*u_R*w_R );
                rhoE_F_m = 0.5*( rho_L*u_L*E_L + u_L*P_L + rho_R*u_R*E_R + u_R*P_R );
		/// Fluxes x-direction
                rho_inv_flux[I1D(i,j,k)]  = ( rho_F_p - rho_F_m )/delta_x;
                rhou_inv_flux[I1D(i,j,k)] = ( rhou_F_p - rhou_F_m )/delta_x;
                rhov_inv_flux[I1D(i,j,k)] = ( rhov_F_p - rhov_F_m )/delta_x;
                rhow_inv_flux[I1D(i,j,k)] = ( rhow_F_p - rhow_F_m )/delta_x;
                rhoE_inv_flux[I1D(i,j,k)] = ( rhoE_F_p - rhoE_F_m )/delta_x;
                /// y-direction j+1/2
                index_L = j;                           index_R = j + 1;
                rho_L   = rho_field[I1D(i,index_L,k)]; rho_R   = rho_field[I1D(i,index_R,k)];
                u_L     = v_field[I1D(i,index_L,k)];   u_R     = v_field[I1D(i,index_R,k)];
                v_L     = u_field[I1D(i,index_L,k)];   v_R     = u_field[I1D(i,index_R,k)];
                w_L     = w_field[I1D(i,index_L,k)];   w_R     = w_field[I1D(i,index_R,k)];
                E_L     = E_field[I1D(i,index_L,k)];   E_R     = E_field[I1D(i,index_R,k)];
                P_L     = P_field[I1D(i,index_L,k)];   P_R     = P_field[I1D(i,index_R,k)];
                /// rho, rhou, rhov, rhow, rhoE
                rho_F_p  = 0.5*( rho_L*u_L + rho_R*u_R );
                rhou_F_p = 0.5*( rho_L*u_L*v_L + rho_R*u_R*v_R );
                rhov_F_p = 0.5*( rho_L*u_L*u_L + P_L + rho_R*u_R*u_R + P_R );
                rhow_F_p = 0.5*( rho_L*u_L*w_L + rho_R*u_R*w_R );
                rhoE_F_p = 0.5*( rho_L*u_L*E_L + u_L*P_L + rho_R*u_R*E_R + u_R*P_R );
		/// y-direction j-1/2
                index_L = j - 1;                       index_R = j;
                rho_L   = rho_field[I1D(i,index_L,k)]; rho_R   = rho_field[I1D(i,index_R,k)];
                u_L     = v_field[I1D(i,index_L,k)];   u_R     = v_field[I1D(i,index_R,k)];
                v_L     = u_field[I1D(i,index_L,k)];   v_R     = u_field[I1D(i,index_R,k)];
                w_L     = w_field[I1D(i,index_L,k)];   w_R     = w_field[I1D(i,index_R,k)];
                E_L     = E_field[I1D(i,index_L,k)];   E_R     = E_field[I1D(i,index_R,k)];
                P_L     = P_field[I1D(i,index_L,k)];   P_R     = P_field[I1D(i,index_R,k)];
                /// rho, rhou, rhov, rhow, rhoE
                rho_F_m  = 0.5*( rho_L*u_L + rho_R*u_R );
                rhou_F_m = 0.5*( rho_L*u_L*v_L + rho_R*u_R*v_R );
                rhov_F_m = 0.5*( rho_L*u_L*u_L + P_L + rho_R*u_R*u_R + P_R );
                rhow_F_m = 0.5*( rho_L*u_L*w_L + rho_R*u_R*w_R );
                rhoE_F_m = 0.5*( rho_L*u_L*E_L + u_L*P_L + rho_R*u_R*E_R + u_R*P_R );
		/// Fluxes y-direction
                rho_inv_flux[I1D(i,j,k)]  += ( rho_F_p - rho_F_m )/delta_y;
                rhou_inv_flux[I1D(i,j,k)] += ( rhou_F_p - rhou_F_m )/delta_y;
                rhov_inv_flux[I1D(i,j,k)] += ( rhov_F_p - rhov_F_m )/delta_y;
                rhow_inv_flux[I1D(i,j,k)] += ( rhow_F_p - rhow_F_m )/delta_y;
                rhoE_inv_flux[I1D(i,j,k)] += ( rhoE_F_p - rhoE_F_m )/delta_y;
                /// z-direction k+1/2
                index_L = k;                           index_R = k + 1;
                rho_L   = rho_field[I1D(i,j,index_L)]; rho_R   = rho_field[I1D(i,j,index_R)];
                u_L     = w_field[I1D(i,j,index_L)];   u_R     = w_field[I1D(i,j,index_R)];
                v_L     = v_field[I1D(i,j,index_L)];   v_R     = v_field[I1D(i,j,index_R)];
                w_L     = u_field[I1D(i,j,index_L)];   w_R     = u_field[I1D(i,j,index_R)];
                E_L     = E_field[I1D(i,j,index_L)];   E_R     = E_field[I1D(i,j,index_R)];
                P_L     = P_field[I1D(i,j,index_L)];   P_R     = P_field[I1D(i,j,index_R)];
                /// rho, rhou, rhov, rhow, rhoE
                rho_F_p  = 0.5*( rho_L*u_L + rho_R*u_R );
                rhou_F_p = 0.5*( rho_L*u_L*w_L + rho_R*u_R*w_R );
                rhov_F_p = 0.5*( rho_L*u_L*v_L + rho_R*u_R*v_R );
                rhow_F_p = 0.5*( rho_L*u_L*u_L + P_L + rho_R*u_R*u_R + P_R );
                rhoE_F_p = 0.5*( rho_L*u_L*E_L + u_L*P_L + rho_R*u_R*E_R + u_R*P_R );
		/// z-direction k-1/2
                index_L = k - 1;                       index_R = k;
                rho_L   = rho_field[I1D(i,j,index_L)]; rho_R   = rho_field[I1D(i,j,index_R)];
                u_L     = w_field[I1D(i,j,index_L)];   u_R     = w_field[I1D(i,j,index_R)];
                v_L     = v_field[I1D(i,j,index_L)];   v_R     = v_field[I1D(i,j,index_R)];
                w_L     = u_field[I1D(i,j,index_L)];   w_R     = u_field[I1D(i,j,index_R)];
                E_L     = E_field[I1D(i,j,index_L)];   E_R     = E_field[I1D(i,j,index_R)];
                P_L     = P_field[I1D(i,j,index_L)];   P_R     = P_field[I1D(i,j,index_R)];
                /// rho, rhou, rhov, rhow, rhoE
                rho_F_m  = 0.5*( rho_L*u_L + rho_R*u_R );
                rhou_F_m = 0.5*( rho_L*u_L*w_L + rho_R*u_R*w_R );
                rhov_F_m = 0.5*( rho_L*u_L*v_L + rho_R*u_R*v_R );
                rhow_F_m = 0.5*( rho_L*u_L*u_L + P_L + rho_R*u_R*u_R + P_R );
                rhoE_F_m = 0.5*( rho_L*u_L*E_L + u_L*P_L + rho_R*u_R*E_R + u_R*P_R );
		/// Fluxes z-direction
                rho_inv_flux[I1D(i,j,k)]  += ( rho_F_p - rho_F_m )/delta_z;
                rhou_inv_flux[I1D(i,j,k)] += ( rhou_F_p - rhou_F_m )/delta_z;
                rhov_inv_flux[I1D(i,j,k)] += ( rhov_F_p - rhov_F_m )/delta_z;
                rhow_inv_flux[I1D(i,j,k)] += ( rhow_F_p - rhow_F_m )/delta_z;
                rhoE_inv_flux[I1D(i,j,k)] += ( rhoE_F_p - rhoE_F_m )/delta_z;
            }
        }
    }

};


////////// MAIN //////////
int main(int argc, char** argv) {

    /// Initialize MPI
    MPI_Init(&argc, &argv);

    /// Process command line arguments
    string configuration_file;
    if( argc >= 2 ) {
        configuration_file = argv[1];
    } else {
        cout << "Proper usage: RHEA.exe configuration_file.yaml" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /// Construct my RHEA
    myRHEA my_RHEA( configuration_file );

    /// Execute my RHEA
    my_RHEA.execute();

    /// Destruct my RHEA ... destructor is called automatically

    /// Finalize MPI
    MPI_Finalize();

    /// Return exit code of program
    return 0;

}

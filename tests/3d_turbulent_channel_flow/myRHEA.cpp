#include "myRHEA.hpp"

using namespace std;


////////// COMPILATION DIRECTIVES //////////
#define _DYNAMICALLY_ADJUST_PRESSURE_GRADIENT_ 0		/// Based on averaged velocity, adjust pressure gradient on the fly

/// Pi number
//const double pi = 2.0*asin( 1.0 );

/// PROBLEM PARAMETERS ///
//const double R_specific = 287.058;				/// Specific gas constant
const double gamma_0    = 1.4;					/// Heat capacity ratio
//const double c_p        = gamma_0*R_specific/( gamma_0 - 1.0 );	/// Isobaric heat capacity
const double delta      = 1.0;					/// Channel half-height
const double Re_tau     = 180.0;				/// Friction Reynolds number
const double Ma         = 3.0e-1;				/// Mach number
//const double Pr         = 0.71;					/// Prandtl number
const double rho_0      = 1.0;					/// Reference density	
const double u_tau      = 1.0;					/// Friction velocity
const double tau_w      = rho_0*u_tau*u_tau;			/// Wall shear stress
//const double mu         = rho_0*u_tau*delta/Re_tau;		/// Dynamic viscosity	
const double nu         = u_tau*delta/Re_tau;			/// Kinematic viscosity	
//const double kappa      = c_p*mu/Pr;				/// Thermal conductivity	
const double Re_b       = pow( Re_tau/0.09, 1.0/0.88 );		/// Bulk (approximated) Reynolds number
const double u_b        = nu*Re_b/( 2.0*delta );		/// Bulk (approximated) velocity
const double P_0        = rho_0*u_b*u_b/( gamma_0*Ma*Ma );	/// Reference pressure
//const double L_x       = 4.0*pi*delta;				/// Streamwise length
//const double L_y       = 2.0*delta;				/// Wall-normal height
//const double L_z       = 4.0*pi*delta/3.0;			/// Spanwise width
const double kappa_vK   = 0.41;					/// von Kármán constant
const double y_0        = nu/( 9.0*u_tau );			/// Smooth-wall roughness
const double u_0        = ( u_tau/kappa_vK )*( log( delta/y_0 ) + ( y_0/delta ) - 1.0 );	/// Volume average of a log-law velocity profile
const double alpha      = 0.5;					/// Magnitude of velocity perturbation
const double u_avg_ref  = 15.67867;				/// Reference averaged streamwise velocity

////////// myRHEA CLASS //////////

void myRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    srand( my_rank );

    /// All (inner, halo, boundary): u, v, w, P and T
    double random_number, y_dist;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                random_number       = 2.0*( (double) rand()/( RAND_MAX ) ) - 1.0;
                y_dist              = min( mesh->y[j], 2.0*delta - mesh->y[j] );
                u_field[I1D(i,j,k)] = ( 2.0*u_0*y_dist/delta ) + alpha*u_0*random_number;
                v_field[I1D(i,j,k)] = 0.0;
                w_field[I1D(i,j,k)] = 0.0;
                P_field[I1D(i,j,k)] = P_0;
                T_field[I1D(i,j,k)] = thermodynamics->calculateTemperatureFromPressureDensity( P_field[I1D(i,j,k)], rho_0 );
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

#if _DYNAMICALLY_ADJUST_PRESSURE_GRADIENT_
    /// Calculate local delta_avg & u_avg values
    double delta_avg_local = 0.0;
    double u_avg_local     = 0.0;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                double delta_y   = 0.5*( mesh->y[j+1] - mesh->y[j-1] );
		delta_avg_local += delta_y;
		u_avg_local     += delta_y*u_field[I1D(i,j,k)];
            }
        }
    }

    /// Calculate global delta_avg & u_avg values and u_avg
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    double delta_avg_global, u_avg_global;
    MPI_Allreduce(&delta_avg_local, &delta_avg_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&u_avg_local, &u_avg_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double u_avg = u_avg_global/delta_avg_global;
#endif

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
#if _DYNAMICALLY_ADJUST_PRESSURE_GRADIENT_
                f_rhou_field[I1D(i,j,k)] = ( u_avg_ref/u_avg )*( tau_w/delta );		    
#else
		f_rhou_field[I1D(i,j,k)] = tau_w/delta;
#endif
                f_rhov_field[I1D(i,j,k)] = 0.0;
                f_rhow_field[I1D(i,j,k)] = 0.0;
                f_rhoE_field[I1D(i,j,k)] = ( -1.0 )*( f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    //f_rhou_field.update();
    //f_rhov_field.update();
    //f_rhow_field.update();
    //f_rhoE_field.update();

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

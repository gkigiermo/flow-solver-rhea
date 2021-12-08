#include "myRHEA.hpp"

using namespace std;

/// Pi number
const double pi = 2.0*asin( 1.0 );

/// PROBLEM PARAMETERS ///
const double rho_b      = 785.49;                               /// Density bottom wall [kg/m3]
const double mu_b       = 0.00011004;				/// Dynamic viscosity bottom wall [Pa s]
const double T_w_b      = 100.0;				/// Temperature bottom wall [K]
const double T_w_t      = 300.0;				/// Temperature top wall [K]
const double P_0        = 4.0e6;                                /// Reference pressure
const double u_tau_b    = 1.0;					/// Friction velocity
const double Re_tau_b   = 100.0;                                /// Friction Reynolds number
const double Ma         = 3.0e-1;                               /// Mach number
const double delta      = ( mu_b*Re_tau_b )/( rho_b*u_tau_b );  /// Channel half-height [m]
const double Re_b       = pow( Re_tau_b/0.09, 1.0/0.88 );	/// Bulk (approximated) Reynolds number
const double u_b        = ( mu_b/rho_b )*Re_b/( 2.0*delta );	/// Bulk (approximated) velocity
const double tau_w_b    = rho_b*u_tau_b*u_tau_b;                /// Wall shear stress
const double nu_b       = mu_b/rho_b; 	                        /// Kinematic viscosity 
const double kappa_vK   = 0.41;                                 /// von Kármán constant
const double y_0        = nu_b/( 9.0*u_tau_b );                 /// Smooth-wall roughness
const double u_0        = ( u_tau_b/kappa_vK )*( log( delta/y_0 ) + ( y_0/delta ) - 1.0 );        /// Volume average of a log-law velocity profile
const double alpha      = 0.25;                                 /// Magnitude of perturbations

/// Pressure-gradient force
double modulation_ratio_old = 1.0;                              /// Previous modulation ratio
double blending_weight      = 0.25;                             /// Blending weight ( new value )


////////// myRHEA CLASS //////////

void myRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    srand( my_rank );

    /// All (inner, halo, boundary): u, v, w, P and T
    double random_number, y_dist;
    //double random_number;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                random_number       = 2.0*( (double) rand()/( RAND_MAX ) ) - 1.0;
                y_dist              = min( mesh->y[j], 2.0*delta - mesh->y[j] );
		u_field[I1D(i,j,k)] = ( 2.0*u_0*y_dist/delta ) + alpha*u_0*random_number;
                v_field[I1D(i,j,k)] = 0.0;
                w_field[I1D(i,j,k)] = 0.0;
                P_field[I1D(i,j,k)] = P_0;
                T_field[I1D(i,j,k)] = T_w_b + ( mesh->y[j]/( 2.0*delta ) )*( T_w_t - T_w_b );
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

    /// Evaluate numerical shear stress at bottom wall

    /// Calculate local values
    double local_sum_u_inner_b = 0.0;
    double local_number_grid_points_b  = 0.0;
    for(int i = topo->iter_bound[_SOUTH_][_INIX_]; i <= topo->iter_bound[_SOUTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_SOUTH_][_INIY_]; j <= topo->iter_bound[_SOUTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_SOUTH_][_INIZ_]; k <= topo->iter_bound[_SOUTH_][_ENDZ_]; k++) {
                /// Sum inner values
                local_sum_u_inner_b += u_field[I1D(i,j+1,k)];
                /// Sum number grid points
                local_number_grid_points_b += 1.0;
            }
        }
    }

    /// Communicate local values to obtain global & average values
    double global_sum_u_inner_b;
    MPI_Allreduce(&local_sum_u_inner_b, &global_sum_u_inner_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_number_grid_points_b;
    MPI_Allreduce(&local_number_grid_points_b, &global_number_grid_points_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_avg_u_inner_b = global_sum_u_inner_b/global_number_grid_points_b;

    /// Calculate delta_y
    double delta_y = 0.5*( mesh->getGloby(1) - mesh->getGloby(0) );

    /// Calculate tau_w_b_numerical
    double tau_w_b_numerical = mu_b*global_avg_u_inner_b/delta_y;
    
    /// Calculate modulation_ratio & update modulation_ratio_old
    double modulation_ratio = blending_weight*( tau_w_b/( tau_w_b_numerical + 1.0e-15 ) ) + ( 1.0 - blending_weight )*modulation_ratio_old;
    modulation_ratio_old    = modulation_ratio;

    //int my_rank, world_size;
    //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //if( my_rank == 0 ) cout << tau_w_b << "  " << tau_w_b_numerical << "  " << modulation_ratio << endl;

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                //f_rhou_field[I1D(i,j,k)] = tau_w_b/delta;
                f_rhou_field[I1D(i,j,k)] = modulation_ratio*tau_w_b/delta;
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

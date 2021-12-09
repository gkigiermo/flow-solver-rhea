#include "myRHEA.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
const double epsilon = 1.0e-15;					/// Small epsilon number (fixed)
const double pi      = 2.0*asin(1.0);				/// pi number (fixed)
const int cout_presicion = 5;		                	/// Output precision (fixed)

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

/// Estimated uniform body force to drive the flow
double controller_output = tau_w_b/delta;		        /// Initialize controller output
double controller_error  = 0.0;			        	/// Initialize controller error
double controller_K_p    = 1.0e-1;		        	/// Controller proportional gain


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
    double local_sum_u_inner_b    = 0.0;
    double local_sum_u_boundary_b = 0.0;
    double local_number_grid_points_b  = 0.0;
    for(int i = topo->iter_bound[_SOUTH_][_INIX_]; i <= topo->iter_bound[_SOUTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_SOUTH_][_INIY_]; j <= topo->iter_bound[_SOUTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_SOUTH_][_INIZ_]; k <= topo->iter_bound[_SOUTH_][_ENDZ_]; k++) {
                /// Sum inner values
                local_sum_u_inner_b += u_field[I1D(i,j+1,k)];
                /// Sum boundary values
                local_sum_u_boundary_b += u_field[I1D(i,j,k)];
                /// Sum number grid points
                local_number_grid_points_b += 1.0;
            }
        }
    }

    /// Communicate local values to obtain global & average values
    double global_sum_u_inner_b;
    MPI_Allreduce(&local_sum_u_inner_b, &global_sum_u_inner_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_sum_u_boundary_b;
    MPI_Allreduce(&local_sum_u_boundary_b, &global_sum_u_boundary_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_number_grid_points_b;
    MPI_Allreduce(&local_number_grid_points_b, &global_number_grid_points_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_avg_u_inner_b    = global_sum_u_inner_b/global_number_grid_points_b;
    double global_avg_u_boundary_b = global_sum_u_boundary_b/global_number_grid_points_b;

    /// Calculate delta_y
    double delta_y = mesh->getGloby(1) - mesh->getGloby(0);

    /// Calculate tau_w_b_numerical
    double tau_w_b_numerical = mu_b*( global_avg_u_inner_b - global_avg_u_boundary_b )/delta_y;
    
    /// Update controller variables
    controller_error   = ( tau_w_b - tau_w_b_numerical )/delta;
    controller_output += controller_K_p*controller_error;

    //int my_rank, world_size;
    //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //if( my_rank == 0 ) cout << tau_w_b << "  " << tau_w_b_numerical << "  " << controller_output << "  " << controller_error << endl;

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                //f_rhou_field[I1D(i,j,k)] = tau_w_b/delta;
                f_rhou_field[I1D(i,j,k)] = controller_output;
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

void myRHEA::execute() {
    
    /// Start timer: execute
    timers->start( "execute" );

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /// Set output (cout) precision
    cout.precision( cout_presicion );

    /// Start RHEA simulation
    if( my_rank == 0 ) cout << "RHEA (v" << version_number << "): START SIMULATION" << endl;

    /// Initialize variables from restart file or by setting initial conditions
    if( use_restart ) {

        /// Initialize from restart file
        this->initializeFromRestart();


    } else {

        /// Set initial conditions
        this->setInitialConditions();
    
        /// Initialize thermodynamics
        this->initializeThermodynamics();

        /// Calculate transport coefficients
        this->calculateTransportCoefficients();

    }

    /// Calculate conserved variables from primitive variables
    this->primitiveToConservedVariables();

    /// Update previous state of conserved variables
    this->updatePreviousStateConservedVariables();    
    
    /// Start timer: time_iteration_loop
    timers->start( "time_iteration_loop" );

    /// Iterate flow solver RHEA in time
    for(int time_iter = current_time_iter; time_iter < final_time_iter; time_iter++) {

        /// Start timer: calculate_time_step
        timers->start( "calculate_time_step" );

        /// Calculate time step
        this->calculateTimeStep();
        if( ( current_time + delta_t ) > final_time ) delta_t = final_time - current_time;

        /// Stop timer: calculate_time_step
        timers->stop( "calculate_time_step" );

        /// Stop timer: execute
        timers->stop( "execute" );

        /// Start timer: output_solver_state
        timers->start( "output_solver_state" );

        /// Print time iteration information (if criterion satisfied)
        if( ( current_time_iter%print_frequency_iter == 0 ) and ( my_rank == 0 ) ) {
            cout << "Time iteration " << current_time_iter << ": " 
                 << "time = " << scientific << current_time << " [s], "
                 << "time-step = " << scientific << delta_t << " [s], "
                 << "wall-clock time = " << scientific << timers->getAccumulatedMaxTime( "execute" )/3600.0 << " [h]" << endl;
        }

        /// Output current state data to file (if criterion satisfied)
        if( current_time_iter%output_frequency_iter == 0 ) this->outputCurrentStateData();

        /// Stop timer: output_solver_state
        timers->stop( "output_solver_state" );

        /// Start timer: execute
        timers->start( "execute" );

        /// Start timer: rk_iteration_loop
        timers->start( "rk_iteration_loop" );

        /// Runge-Kutta time-integration steps
        for(int rk_time_stage = 1; rk_time_stage <= rk_number_stages; rk_time_stage++) {

            /// Start timer: calculate_thermophysical_properties
            timers->start( "calculate_thermophysical_properties" );

            /// Calculate transport coefficients
            this->calculateTransportCoefficients();

            /// Stop timer: calculate_thermophysical_properties
            timers->stop( "calculate_thermophysical_properties" );

            /// Start timer: calculate_inviscid_fluxes
            timers->start( "calculate_inviscid_fluxes" );

            /// Calculate inviscid fluxes ( x-direction is required to be first! )
            this->calculateInviscidFluxesX();
            this->calculateInviscidFluxesY();
            this->calculateInviscidFluxesZ();

            /// Stop timer: calculate_inviscid_fluxes
            timers->stop( "calculate_inviscid_fluxes" );

            /// Start timer: calculate_viscous_fluxes
            timers->start( "calculate_viscous_fluxes" );

            /// Calculate viscous fluxes
            this->calculateViscousFluxes();

            /// Stop timer: calculate_viscous_fluxes
            timers->stop( "calculate_viscous_fluxes" );

            /// Start timer: calculate_source_terms
            timers->start( "calculate_source_terms" );

            /// Calculate source terms
            this->calculateSourceTerms();

            /// Stop timer: calculate_source_terms
            timers->stop( "calculate_source_terms" );

            /// Start timer: time_advance_conserved_variables
            timers->start( "time_advance_conserved_variables" );

            /// Advance conserved variables in time
            this->timeAdvanceConservedVariables(rk_time_stage);

            /// Stop timer: time_advance_conserved_variables
            timers->stop( "time_advance_conserved_variables" );

            /// Start timer: conserved_to_primitive_variables
            timers->start( "conserved_to_primitive_variables" );

            /// Calculate primitive variables from conserved variables
            this->conservedToPrimitiveVariables();

            /// Stop timer: conserved_to_primitive_variables
            timers->stop( "conserved_to_primitive_variables" );

            /// Start timer: calculate_thermodynamics_from_primitive_variables
            timers->start( "calculate_thermodynamics_from_primitive_variables" );

            /// Calculate thermodynamics from primitive variables
            this->calculateThermodynamicsFromPrimitiveVariables();

            /// Stop timer: calculate_thermodynamics_from_primitive_variables
            timers->stop( "calculate_thermodynamics_from_primitive_variables" );

            /// Energy source/sink to maintain thermodynamic pressure to P_0
	    double rho, ke, e, c_v, c_p;
            for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
                for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                    for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                        if( x_field[I1D(i,j,k)] < ( x_0 + L_x/num_grid_x ) ) {
                            P_field[I1D(i,j,k)] = P_0;
                            thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho, e, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)] );
                            rho_field[I1D(i,j,k)] = rho; 
                            ke = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) );
                            E_field[I1D(i,j,k)] = e + ke;
                            sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                            thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                            c_v_field[I1D(i,j,k)] = c_v;
                            c_p_field[I1D(i,j,k)] = c_p;		       	
                            rhou_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*u_field[I1D(i,j,k)]; 
                            rhov_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*v_field[I1D(i,j,k)]; 
                            rhow_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*w_field[I1D(i,j,k)]; 
                            rhoE_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*E_field[I1D(i,j,k)]; 
		        }
	            }
                }
            }

            /// Start timer: update_boundaries
            timers->start( "update_boundaries" );

            /// Update boundary values
            this->updateBoundaries();

            /// Stop timer: update_boundaries
            timers->stop( "update_boundaries" );

        }

        /// Stop timer: rk_iteration_loop
        timers->stop( "rk_iteration_loop" );

        /// Start timer: update_time_averaged_quantities
        timers->start( "update_time_averaged_quantities" );

        /// Update time-averaged quantities
        if( time_averaging_active ) this->updateTimeAveragedQuantities();

        /// Stop timer: update_time_averaged_quantities
        timers->stop( "update_time_averaged_quantities" );

        /// Start timer: update_previous_state_conserved_variables
        timers->start( "update_previous_state_conserved_variables" );

        /// Update previous state of conserved variables
        this->updatePreviousStateConservedVariables();

        /// Update time and time iteration
        current_time += delta_t;
        current_time_iter += 1;

        /// Check if simulation is completed: current_time > final_time
        if( current_time >= final_time ) break;

        /// Stop timer: update_previous_state_conserved_variables
        timers->stop( "update_previous_state_conserved_variables" );

    }

    /// Stop timer: time_iteration_loop
    timers->stop( "time_iteration_loop" );

    /// Print timers information
    if( print_timers ) timers->printTimers( timers_information_file );

    /// Print time advancement information
    if( my_rank == 0 ) {
        cout << "Time advancement completed -> " 
             << "iteration = " << current_time_iter << ", "
             << "time = " << scientific << current_time << " [s]" << endl;
        }

    /// Output current state data to file
    this->outputCurrentStateData();

    /// End RHEA simulation
    if( my_rank == 0 ) cout << "RHEA (v" << version_number << "): END SIMULATION" << endl;
    
    /// Stop timer: execute
    timers->stop( "execute" );

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

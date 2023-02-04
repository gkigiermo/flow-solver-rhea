#include "myRHEA.hpp"

#ifdef _OPENACC
#include <openacc.h>
#endif

using namespace std;

////////// FIXED PARAMETERS //////////
const double epsilon = 1.0e-15;					/// Small epsilon number (fixed)
const double pi      = 2.0*asin(1.0);				/// pi number (fixed)
const int cout_precision = 5;		                	/// Output precision (fixed)

/// PROBLEM PARAMETERS ///
const double T_c        = 126.192;				/// Nitrogen critical temperature [K]
const double P_c        = 3395800.0;				/// Nitrogen critical pressure [Pa]
const double T_bw       = 0.75*T_c;                             /// Temperature bottom wall [K]
const double T_tw       = 1.5*T_c;                              /// Temperature top wall [K]
const double P_b        = 2.0*P_c;                              /// Bulk pressure [Pa]
const double delta      = 100.0e-6;				/// Channel half-height [m]
const double Re_tau_bw  = 100.0;                                /// Friction Reynolds number bottom wall
const double rho_bw     = 839.39;                               /// Density bottom wall [kg/m3]
const double mu_bw      = 0.00016;				/// Dynamic viscosity bottom wall [Pa s]
const double nu_bw      = mu_bw/rho_bw; 	                /// Kinematic viscosity bottom wall [m2/s]
const double u_tau_bw   = ( mu_bw*Re_tau_bw )/( rho_bw*delta );	/// Friction velocity bottom wall [m/s]
const double tau_bw     = rho_bw*u_tau_bw*u_tau_bw;             /// Wall shear stress bottom wall [kg/(m s2)]
const double kappa_vK   = 0.41;                                 /// von Kármán constant
const double y_0        = nu_bw/( 9.0*u_tau_bw );               /// Smooth-wall roughness bottom wall [m]
const double u_0        = ( u_tau_bw/kappa_vK )*( log( delta/y_0 ) + ( y_0/delta ) - 1.0 );   /// Volume-average of a log-law velocity profile [m/s]
const double alpha      = 0.25;                                 /// Magnitude of perturbations

/// Estimated uniform body force to drive the flow
double controller_output = tau_bw/delta;		        /// Initialize controller output
double controller_error  = 0.0;			        	/// Initialize controller error
double controller_K_p    = 1.0e-1;		        	/// Controller proportional gain

/// Control bulk pressure to maintain the fixed target value P_b
bool control_bulk_pressure = true;

/// Control temperature to T_b_w <= T <= T_t_w
bool control_temperature = true;


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
                P_field[I1D(i,j,k)] = P_b;
                T_field[I1D(i,j,k)] = T_bw + ( mesh->y[j]/( 2.0*delta ) )*( T_tw - T_bw );
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
    double local_sum_u_boundary_bw = 0.0;
    double local_sum_u_inner_bw    = 0.0;
    double local_number_grid_points_bw  = 0.0;
    for(int i = topo->iter_bound[_SOUTH_][_INIX_]; i <= topo->iter_bound[_SOUTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_SOUTH_][_INIY_]; j <= topo->iter_bound[_SOUTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_SOUTH_][_INIZ_]; k <= topo->iter_bound[_SOUTH_][_ENDZ_]; k++) {
//		if( abs( avg_u_field[I1D(i,j,k)] ) > 0.0 ) {
//                    /// Sum boundary values
//                    local_sum_u_boundary_bw += avg_u_field[I1D(i,j,k)];
//                    /// Sum inner values
//                    local_sum_u_inner_bw    += avg_u_field[I1D(i,j+1,k)];
//		} else {
                    /// Sum boundary values
                    local_sum_u_boundary_bw += u_field[I1D(i,j,k)];
                    /// Sum inner values
                    local_sum_u_inner_bw    += u_field[I1D(i,j+1,k)];
//		}
                /// Sum number grid points
                local_number_grid_points_bw += 1.0;
            }
        }
    }

    /// Communicate local values to obtain global & average values
    double global_sum_u_boundary_bw;
    MPI_Allreduce(&local_sum_u_boundary_bw, &global_sum_u_boundary_bw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_sum_u_inner_bw;
    MPI_Allreduce(&local_sum_u_inner_bw, &global_sum_u_inner_bw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_number_grid_points_bw;
    MPI_Allreduce(&local_number_grid_points_bw, &global_number_grid_points_bw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_avg_u_boundary_bw = global_sum_u_boundary_bw/global_number_grid_points_bw;
    double global_avg_u_inner_bw    = global_sum_u_inner_bw/global_number_grid_points_bw;

    /// Calculate delta_y
    double delta_y = mesh->getGloby(1) - mesh->getGloby(0);

    /// Calculate tau_bw_numerical
    double tau_bw_numerical = mu_bw*( global_avg_u_inner_bw - global_avg_u_boundary_bw )/delta_y;
    
    /// Update controller variables
    controller_error   = ( tau_bw - tau_bw_numerical )/delta;
    controller_output += controller_K_p*controller_error;

    //int my_rank, world_size;
    //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //if( my_rank == 0 ) cout << tau_bw << "  " << tau_w_b_numerical << "  " << controller_output << "  " << controller_error << endl;

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                //f_rhou_field[I1D(i,j,k)] = tau_bw/delta;
                f_rhou_field[I1D(i,j,k)] = controller_output;
                f_rhov_field[I1D(i,j,k)] = 0.0;
                f_rhow_field[I1D(i,j,k)] = 0.0;
                f_rhoE_field[I1D(i,j,k)] = 0.0;
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
    cout.precision( cout_precision );

    /// Start RHEA simulation
    if( my_rank == 0 ) cout << "RHEA (v" << version_number << "): START SIMULATION" << endl;

    /// Initialize variables from restart file or by setting initial conditions
    if( use_restart ) {

        /// Initialize from restart file
        this->initializeFromRestart();

        if( artificial_compressibility_method ) {

            /// Calculate thermodynamic (bulk) pressure
            P_thermo = this->calculateVolumeAveragedPressure();

	    /// Calculate artificially modified thermodynamics
            this->calculateArtificiallyModifiedThermodynamics();

	    /// Calculate artificially modified transport coefficients
            this->calculateArtificiallyModifiedTransportCoefficients();

	}

    } else {

        /// Set initial conditions
        this->setInitialConditions();

        /// Initialize thermodynamics
        this->initializeThermodynamics();

        if( artificial_compressibility_method ) {

            /// Calculate thermodynamic (bulk) pressure
            P_thermo = this->calculateVolumeAveragedPressure();

	    /// Calculate artificially modified thermodynamics
            this->calculateArtificiallyModifiedThermodynamics();	    

	    /// Calculate artificially modified transport coefficients
            this->calculateArtificiallyModifiedTransportCoefficients();

	} else {

            /// Calculate transport coefficients
            this->calculateTransportCoefficients();

	}

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

            if( artificial_compressibility_method ) {

	        /// Calculate artificially modified transport coefficients
                this->calculateArtificiallyModifiedTransportCoefficients();

	    } else {

                /// Calculate transport coefficients
                this->calculateTransportCoefficients();

	    }

            /// Stop timer: calculate_thermophysical_properties
            timers->stop( "calculate_thermophysical_properties" );

            /// Start timer: calculate_inviscid_fluxes
            timers->start( "calculate_inviscid_fluxes" );

            /// Calculate inviscid fluxes
            this->calculateInviscidFluxes();

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

            if( artificial_compressibility_method ) {

                /// Calculate thermodynamic (bulk) pressure
                P_thermo = this->calculateVolumeAveragedPressure();

                /// Calculate artificially modified thermodynamics
                this->calculateArtificiallyModifiedThermodynamics();	    

	    }

            /// Stop timer: calculate_thermodynamics_from_primitive_variables
            timers->stop( "calculate_thermodynamics_from_primitive_variables" );

            /// Start: bulk pressure explicitly modified to maintain the fixed target value P_b
            if( control_bulk_pressure ) {

                /// Calculate local values
                double local_sum_PV = 0.0;
                double local_sum_V  = 0.0;
                double volume       = 0.0;
                double delta_x, delta_y, delta_z;
                for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
                    for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                        for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                            /// Geometric stuff
                            delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                            delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                            delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                            /// Calculate volume
                            volume = delta_x*delta_y*delta_z; 
                            /// Sum P*V values
                            local_sum_PV += P_field[I1D(i,j,k)]*volume;
                            /// Sum V values
                            local_sum_V += volume;
                        }
                    }
                }

                /// Communicate local values to obtain global & average values
                double global_sum_PV;
                MPI_Allreduce(&local_sum_PV, &global_sum_PV, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                double global_sum_V;
                MPI_Allreduce(&local_sum_V, &global_sum_V, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                double global_avg_P = global_sum_PV/global_sum_V;

                /// Modify P values
                double ratio_P_b_target_P_b_numerical = P_b/( global_avg_P + epsilon );
                for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
                    for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                        for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                            /// Sum P*V values
                            P_field[I1D(i,j,k)] *= ratio_P_b_target_P_b_numerical;
                        }
                    }
                }

            }
            /// Stop: bulk pressure explicitly modified to maintain the fixed target value P_b
            
	    /// Start: temperature explicitly modified to T_b_w <= T <= T_t_w
            if( control_temperature ) {
	  
	        /// Mantain T to T_b_w <= T <= T_t_w
                for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
                    for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                        for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
			    T_field[I1D(i,j,k)] = max( T_bw, min( T_tw, T_field[I1D(i,j,k)] ) );
			}
		    }
		}

            } 
	    /// Stop: temperature explicitly modified to T_b_w <= T <= T_t_w

	    /// Start: recalculate rho, E, sos, c_v and c_p from P and T
            if( control_bulk_pressure || control_temperature ) {

                double u, v, w, P, T;
                double rho, e, ke, E;
                double c_v, c_p;
                for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
                    for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                        for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
		            /// Obtain primitive variables
                            u = u_field[I1D(i,j,k)];
                            v = v_field[I1D(i,j,k)];
                            w = w_field[I1D(i,j,k)];
                            P = P_field[I1D(i,j,k)];
                            T = T_field[I1D(i,j,k)];
		            /// Update rho, e, ke and E
                            thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho, e, P, T );
                            ke = 0.5*( u*u + v*v + w*w );
                            E  = e + ke;
		            /// Update conserved variables
			    rho_field[I1D(i,j,k)]  = rho;
                            rhou_field[I1D(i,j,k)] = rho*u;
                            rhov_field[I1D(i,j,k)] = rho*v;
                            rhow_field[I1D(i,j,k)] = rho*w;
                            rhoE_field[I1D(i,j,k)] = rho*E;
		            /// Update thermodynamics
		            if( artificial_compressibility_method ) {
                                sos_field[I1D(i,j,k)] = ( 1.0/( alpha + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T, rho );
                                thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T, rho );
			    } else {
                                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P, T, rho );
                                thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P, T, rho );
			    }
                            c_v_field[I1D(i,j,k)] = c_v;
                            c_p_field[I1D(i,j,k)] = c_p;
			}
		    }
		}

	    }
	    /// Stop: recalculate rho, E, sos, c_v and c_p from P and T

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

#ifdef _OPENACC
    /// OpenACC distribution on multiple accelerators (GPU)
    acc_device_t my_device_type;
    int num_devices, gpuId, local_rank;
    MPI_Comm shmcomm;    

    MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm );
    MPI_Comm_rank( shmcomm, &local_rank );           
    my_device_type = acc_get_device_type();                      
    num_devices = acc_get_num_devices( my_device_type );
    gpuId = local_rank % num_devices;
    acc_set_device_num( gpuId, my_device_type );
//    /// OpenACC distribution on multiple accelerators (GPU)
//    acc_device_t device_type = acc_get_device_type();
//    if ( acc_device_nvidia == device_type ) {
//       int ngpus = acc_get_num_devices( acc_device_nvidia );
//       int devicenum = atoi( getenv( "OMPI_COMM_WORLD_LOCAL_RANK" ) );
//       acc_set_device_num( devicenum, acc_device_nvidia );
//    }
//    acc_init(device_type);
#endif

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

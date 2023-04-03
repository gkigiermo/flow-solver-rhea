#include "myRHEA.hpp"

#ifdef _OPENACC
#include <openacc.h>
#endif

using namespace std;

/// Pi number
const double epsilon = 1.0e-15;                                 /// Small epsilon number (fixed)
const double pi      = 2.0*asin(1.0);                           /// pi number (fixed)
const int cout_precision = 5;                                   /// Output precision (fixed)

/// PROBLEM PARAMETERS ///
const double R_specific = 287.0;      	                        /// Specific gas constant [J/(kg K)]
const double gamma_0    = 1.4;                                  /// Heat capacity ratio [-]
const double c_p        = gamma_0*R_specific/( gamma_0 - 1.0 ); /// Isobaric heat capacity [J/(kg K)]
const double Pr         = 0.71;                                 /// Prandtl number [-]
const double eta        = 0.6;                                  /// Temperature difference ratio [-]
const double T_0        = 600.0;                                /// Reference temperature [K]
const double T_hw       = T_0*( 1.0 + eta );                    /// Hot-wall temperature [K]
const double T_cw       = T_0*( 1.0 - eta );                    /// Cold-wall temperature [K]
const double P_0        = 101325.0;                             /// Reference pressure [Pa]
const double L          = 1.0;		                        /// Cavity size [m]
const double g          = 9.81;                                 /// Gravitational constant [m/s2]
const double T_star     = 273.0;                 	        /// Sutherland's reference temperature [K]
const double S          = 110.5;                 	        /// Sutherland's dynamic viscosity constant [K]
const double mu_star    = 1.68e-5;                 	        /// Sutherland's reference dynamic viscosity [Pa s]
const double alpha      = 0.01;                                 /// Magnitude of temperature perturbation

/// Control bulk pressure to maintain the fixed target value P_b
bool control_bulk_pressure = true;

/// Control temperature to T_cw <= T <= T_hw
bool control_temperature = true;

/// Control w-velocity w = 0.0
bool control_w_velocity = true;

////////// myRHEA CLASS //////////

void myRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    srand( my_rank );

    /// All (inner, halo, boundary): u, v, w, P and T
    double random_number;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                    u_field[I1D(i,j,k)] = 0.0;
                    v_field[I1D(i,j,k)] = 0.0;
                    w_field[I1D(i,j,k)] = 0.0;
                    P_field[I1D(i,j,k)] = P_0;
                    random_number       = 2.0*( (double) rand()/( RAND_MAX ) ) - 1.0;
                    T_field[I1D(i,j,k)] = ( ( T_cw - T_hw )/L )*mesh->x[i] + T_hw + alpha*random_number*T_0;
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
                f_rhov_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*( -1.0 )*g;
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

void myRHEA::calculateTransportCoefficients() { 
    
    /// All (inner, halo, boundary) points: mu and kappa
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
		mu_field[I1D(i,j,k)]    = mu_star*pow( T_field[I1D(i,j,k)]/T_star, 3.0/2.0 )*( ( T_star + S )/( T_field[I1D(i,j,k)] + S ) );
                kappa_field[I1D(i,j,k)] = ( mu_field[I1D(i,j,k)]*gamma_0*R_specific )/( ( gamma_0 - 1.0 )*Pr );
            }
        }
    }

    /// Update halo values
    //mu_field.update();
    //kappa_field.update();

};

void myRHEA::calculateArtificiallyModifiedTransportCoefficients() {

    /// All (inner, halo, boundary) points: mu and kappa
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
		mu_field[I1D(i,j,k)]    = mu_star*pow( T_field[I1D(i,j,k)]/T_star, 3.0/2.0 )*( ( T_star + S )/( T_field[I1D(i,j,k)] + S ) );
                kappa_field[I1D(i,j,k)] = ( mu_field[I1D(i,j,k)]*gamma_0*R_specific )/( ( gamma_0 - 1.0 )*Pr );
            }
        }
    }

    /// Update halo values
    //mu_field.update();
    //kappa_field.update();

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

            /// Calculate alpha value of artificial compressibility method
            alpha_acm = this->calculateAlphaArtificialCompressibilityMethod();

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

            /// Calculate alpha value of artificial compressibility method
            alpha_acm = this->calculateAlphaArtificialCompressibilityMethod();

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

                /// Calculate alpha value of artificial compressibility method
                alpha_acm = this->calculateAlphaArtificialCompressibilityMethod();

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
                double ratio_P_0_target_P_0_numerical = P_0/( global_avg_P + epsilon );
                for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
                    for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                        for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                            /// Sum P*V values
                            P_field[I1D(i,j,k)] *= ratio_P_0_target_P_0_numerical;
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
                            T_field[I1D(i,j,k)] = max( T_cw, min( T_hw, T_field[I1D(i,j,k)] ) );
                        }
                    }
                }

            }
            /// Stop: temperature explicitly modified to T_b_w <= T <= T_t_w

            /// Start: w-velocity explicitly modified to w = 0.0
            if( control_w_velocity ) {

                /// Mantain T to T_b_w <= T <= T_t_w
                for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
                    for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                        for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                            w_field[I1D(i,j,k)] = 0.0;
                        }
                    }
                }

            }
            /// Stop: w-velocity explicitly modified to w = 0.0

            /// Start: recalculate rho, E, sos, c_v and c_p from P and T
            if( control_bulk_pressure || control_temperature || control_w_velocity ) {

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
                                sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T, rho );
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

#include "FlowSolverRHEA.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
const double rk_order = 3;			// Time-integration Runge-Kutta order
const double epsilon  = 1.0e-15;		// Small epsilon number
const double Pi       = 2.0*asin(1.0);		// Pi number

////////// FlowSolverRHEA CLASS //////////
FlowSolverRHEA::FlowSolverRHEA() {};

FlowSolverRHEA::FlowSolverRHEA(const string name_configuration_file) : configuration_file(name_configuration_file) {

    /// Read configuration file
    // ... work in progress (YAML language)	

    // The lines below need to be introduced via configuration (input) file !!
    const int _DEBUG_    = 0;
    const double Re_tau  = 100.0;				// Friction Reynolds number [-]
    const double delta   = 1.0;					// Channel half-height [m]
    const double u_tau   = 1.0;					// Friction velocity [m/s]
    const double rho_ref = 1.0;					// Reference density [kg/m3]
    const double P_ref   = 101325.0;				// Reference pressure [Pa]

    R_specific        = 287.058;
    gamma             = 1.4;
    mu                = rho_ref*u_tau*delta/Re_tau;
    kappa             = 0.0;
    x_0               = 0.0;
    y_0               = 0.0;
    z_0               = 0.0;
    L_x               = 4.0*Pi*delta;
    L_y               = 2.0*delta;
    L_z               = 4.0*Pi*delta/3.0;
    initial_time      = 0.0;
    current_time      = 0.0;
    final_time        = 1.0e3;
    num_grid_x        = 64;
    num_grid_y        = 64;
    num_grid_z        = 64;
    CFL               = 0.9;
    initial_time_iter = 0;
    current_time_iter = 0;
    final_time_iter   = 1e6;
    output_iter       = 1e2;
    bocos[_WEST_]     = _PERIODIC_;
    bocos[_EAST_]     = _PERIODIC_;
    bocos[_SOUTH_]    = _PERIODIC_;
    bocos[_NORTH_]    = _PERIODIC_;
    bocos[_BACK_]     = _PERIODIC_;
    bocos[_FRONT_]    = _PERIODIC_;
    np_x              = 2;
    np_y              = 1;
    np_z              = 1;
    // The lines above need to be introduced via configuration (input) file !!

    /// Initialize (construct) computational domain
    dom = new domain(L_x, L_y, L_z, x_0, y_0, z_0, RHEA_NX, RHEA_NY, RHEA_NZ);

    /// Add boundary conditions to computational domain
    dom->updateBocos(bocos);

    /// Initialize (construct) communication scheme
    topo = new comm_scheme(dom, np_x, np_y, np_z);
    //if(topo->getRank() == 0) dom->printDomain();
    //for(int p = 0; p < np_x*np_y*np_z; p++) topo->printCommSchemeToFile(p);

    // The lines below are temporary ... will need to be removed!
    _lNx_ = topo->getlNx();
    _lNy_ = topo->getlNy();
    _lNz_ = topo->getlNz();

    /// Set parallel topology primitive, conserved and thermodynamic variables	
    rho_field.setTopology(topo);
    u_field.setTopology(topo);
    v_field.setTopology(topo);
    w_field.setTopology(topo);
    E_field.setTopology(topo);
    rhou_field.setTopology(topo);
    rhov_field.setTopology(topo);
    rhow_field.setTopology(topo);
    rhoE_field.setTopology(topo);
    P_field.setTopology(topo);
    T_field.setTopology(topo);
    sos_field.setTopology(topo);

    /// Set parallel topology time-integration variables	
    rho_0_field.setTopology(topo);
    rhou_0_field.setTopology(topo);
    rhov_0_field.setTopology(topo);
    rhow_0_field.setTopology(topo);
    rhoE_0_field.setTopology(topo);    

    /// Set parallel topology time-integration fluxes	
    rho_rk1_flux.setTopology(topo);
    rho_rk2_flux.setTopology(topo);
    rho_rk3_flux.setTopology(topo);
    rhou_rk1_flux.setTopology(topo);
    rhou_rk2_flux.setTopology(topo);
    rhou_rk3_flux.setTopology(topo);    
    rhov_rk1_flux.setTopology(topo);
    rhov_rk2_flux.setTopology(topo);
    rhov_rk3_flux.setTopology(topo);    
    rhow_rk1_flux.setTopology(topo);
    rhow_rk2_flux.setTopology(topo);
    rhow_rk3_flux.setTopology(topo);
    rhoE_rk1_flux.setTopology(topo);
    rhoE_rk2_flux.setTopology(topo);
    rhoE_rk3_flux.setTopology(topo);

    /// Set parallel topology inviscid fluxes	
    rho_inv_flux.setTopology(topo);
    rhou_inv_flux.setTopology(topo);
    rhov_inv_flux.setTopology(topo);
    rhow_inv_flux.setTopology(topo);
    rhoE_inv_flux.setTopology(topo);

    /// Set parallel topology viscous fluxes	
    rhou_vis_flux.setTopology(topo);
    rhov_vis_flux.setTopology(topo);
    rhow_vis_flux.setTopology(topo);
    rhoE_vis_flux.setTopology(topo);

    /// Set parallel topology volumetric & external forces
    f_rhou_field.setTopology(topo);
    f_rhov_field.setTopology(topo);
    f_rhow_field.setTopology(topo);
    f_rhoE_field.setTopology(topo);

};

FlowSolverRHEA::~FlowSolverRHEA() {

    /// Free boundary conditions, domain and topology
    if(bocos != NULL) free(bocos);	
    if(dom != NULL) free(dom);	
    if(topo != NULL) free(topo);

};

void FlowSolverRHEA::setInitialConditions() {

    /// All (inner & boundary) points: u, v, w, P and T
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_INIX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_INIY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_INIZ_]; k++) {
                u_field[I1D(i,j,k)] = 0.0;
                v_field[I1D(i,j,k)] = 0.0;
                w_field[I1D(i,j,k)] = 0.0;
                P_field[I1D(i,j,k)] = 0.0;
                T_field[I1D(i,j,k)] = 0.0;
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

void FlowSolverRHEA::updateBoundaries() {

    /// West points: rho, rhou, rhov, rhow and rhoE
    for(int i = topo->iter_bound[_WEST_][_INIX_]; i <= topo->iter_bound[_WEST_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_WEST_][_INIY_]; j <= topo->iter_bound[_WEST_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_WEST_][_INIZ_]; k <= topo->iter_bound[_WEST_][_ENDZ_]; k++) {
                // ... work in progress: implement different types of boundaries!
                rho_field[I1D(i,j,k)]  = 0.0;
                rhou_field[I1D(i,j,k)] = 0.0;
                rhov_field[I1D(i,j,k)] = 0.0;
                rhow_field[I1D(i,j,k)] = 0.0;
                rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

    /// East points: rho, rhou, rhov, rhow and rhoE
    for(int i = topo->iter_bound[_EAST_][_INIX_]; i <= topo->iter_bound[_EAST_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_EAST_][_INIY_]; j <= topo->iter_bound[_EAST_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_EAST_][_INIZ_]; k <= topo->iter_bound[_EAST_][_ENDZ_]; k++) {
                // ... work in progress: implement different types of boundaries!
                rho_field[I1D(i,j,k)]  = 0.0;
                rhou_field[I1D(i,j,k)] = 0.0;
                rhov_field[I1D(i,j,k)] = 0.0;
                rhow_field[I1D(i,j,k)] = 0.0;
                rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

    /// South points: rho, rhou, rhov, rhow and rhoE
    for(int i = topo->iter_bound[_SOUTH_][_INIX_]; i <= topo->iter_bound[_SOUTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_SOUTH_][_INIY_]; j <= topo->iter_bound[_SOUTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_SOUTH_][_INIZ_]; k <= topo->iter_bound[_SOUTH_][_ENDZ_]; k++) {
                // ... work in progress: implement different types of boundaries!
                rho_field[I1D(i,j,k)]  = 0.0;
                rhou_field[I1D(i,j,k)] = 0.0;
                rhov_field[I1D(i,j,k)] = 0.0;
                rhow_field[I1D(i,j,k)] = 0.0;
                rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

    /// North points: rho, rhou, rhov, rhow and rhoE
    for(int i = topo->iter_bound[_NORTH_][_INIX_]; i <= topo->iter_bound[_NORTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_NORTH_][_INIY_]; j <= topo->iter_bound[_NORTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_NORTH_][_INIZ_]; k <= topo->iter_bound[_NORTH_][_ENDZ_]; k++) {
                // ... work in progress: implement different types of boundaries!
                rho_field[I1D(i,j,k)]  = 0.0;
                rhou_field[I1D(i,j,k)] = 0.0;
                rhov_field[I1D(i,j,k)] = 0.0;
                rhow_field[I1D(i,j,k)] = 0.0;
                rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

    /// Back points: rho, rhou, rhov, rhow and rhoE
    for(int i = topo->iter_bound[_BACK_][_INIX_]; i <= topo->iter_bound[_BACK_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_BACK_][_INIY_]; j <= topo->iter_bound[_BACK_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_BACK_][_INIZ_]; k <= topo->iter_bound[_BACK_][_ENDZ_]; k++) {
                // ... work in progress: implement different types of boundaries!
                rho_field[I1D(i,j,k)]  = 0.0;
                rhou_field[I1D(i,j,k)] = 0.0;
                rhov_field[I1D(i,j,k)] = 0.0;
                rhow_field[I1D(i,j,k)] = 0.0;
                rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

    /// Front points: rho, rhou, rhov, rhow and rhoE
    for(int i = topo->iter_bound[_FRONT_][_INIX_]; i <= topo->iter_bound[_FRONT_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_FRONT_][_INIY_]; j <= topo->iter_bound[_FRONT_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_FRONT_][_INIZ_]; k <= topo->iter_bound[_FRONT_][_ENDZ_]; k++) {
                // ... work in progress: implement different types of boundaries!
                rho_field[I1D(i,j,k)]  = 0.0;
                rhou_field[I1D(i,j,k)] = 0.0;
                rhov_field[I1D(i,j,k)] = 0.0;
                rhow_field[I1D(i,j,k)] = 0.0;
                rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

};

void FlowSolverRHEA::initializeThermodynamics() {

    /// Ideal-gas model:
    /// rho = P/(R_specific*T) is density
    /// e = P/(rho*(gamma - 1)) is specific internal energy
    /// ke = (u*u + v*v + w*w)/2 is specific kinetic energy
    /// E = e + ke is total energy
    /// sos = sqrt(gamma*P/rho) is speed of sound

    /// All (inner & boundary) points: rho, e, ke, E and sos
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_INIX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_INIY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_INIZ_]; k++) {
                rho_field[I1D(i,j,k)] = ( 1.0/( R_specific*T_field[I1D(i,j,k)] ) )*P_field[I1D(i,j,k)];
                double e              = ( 1.0/( rho_field[I1D(i,j,k)]*( gamma - 1.0 ) ) )*P_field[I1D(i,j,k)];
                double ke             = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) );
                E_field[I1D(i,j,k)]   = e + ke;
                sos_field[I1D(i,j,k)] = sqrt( gamma*( ( 1.0/rho_field[I1D(i,j,k)] )*P_field[I1D(i,j,k)] ) );
            }
        }
    }

    /// Update halo values
    rho_field.update();
    E_field.update();
    sos_field.update();

};

void FlowSolverRHEA::primitiveToConservedVariables() {

    /// All (inner & boundary) points: rhou, rhov, rhow and rhoE
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_INIX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_INIY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_INIZ_]; k++) {
                rhou_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*u_field[I1D(i,j,k)]; 
                rhov_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*v_field[I1D(i,j,k)]; 
                rhow_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*w_field[I1D(i,j,k)]; 
                rhoE_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*E_field[I1D(i,j,k)]; 
            }
        }
    }

    /// Update halo values
    rhou_field.update();
    rhov_field.update();
    rhow_field.update();
    rhoE_field.update();

};

void FlowSolverRHEA::updatePreviousStateConservedVariables() {

    /// All (inner & boundary) points: rho_0, rhou_0 rhov_0, rhow_0 and rhoE_0
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_INIX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_INIY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_INIZ_]; k++) {
                rho_0_field[I1D(i,j,k)]  = rho_field[I1D(i,j,k)]; 
                rhou_0_field[I1D(i,j,k)] = rhou_field[I1D(i,j,k)]; 
                rhov_0_field[I1D(i,j,k)] = rhov_field[I1D(i,j,k)]; 
                rhow_0_field[I1D(i,j,k)] = rhow_field[I1D(i,j,k)]; 
                rhoE_0_field[I1D(i,j,k)] = rhoE_field[I1D(i,j,k)]; 
            }
        }
    }

    /// Update halo values
    rho_0_field.update();
    rhou_0_field.update();
    rhov_0_field.update();
    rhow_0_field.update();
    rhoE_0_field.update();

};

void FlowSolverRHEA::calculateTimeStep() {

    /// Inviscid time step size for explicit schemes:
    /// E. F. Toro.
    /// Riemann solvers and numerical methods for fluid dynamics.
    /// Springer, 2009.

    /// Viscous time step size for explicit schemes:
    /// E. Turkel, R.C. Swanson, V. N. Vatsa, J.A. White.
    /// Multigrid for hypersonic viscous two- and three-dimensional flows.
    /// NASA Contractor Report 187603, 1991.

    /// Calculate specific heat capacities
    const double c_v = R_specific/( gamma - 1.0 );
    const double c_p = c_v*gamma;

    /// Calculate Prandtl (Pr) number
    double Pr = gamma;
    if(kappa > epsilon) Pr = c_p*mu/kappa;

    /// Initialize to largest double value
    double local_delta_t = numeric_limits<double>::max();

    /// Inner points: find minimum (local) delta_t
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_INIX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_INIY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_INIZ_]; k++) {
                /// Geometric stuff
                double delta_x = 0.5*( dom->x[i+1] - dom->x[i-1] ); 
                double delta_y = 0.5*( dom->y[j+1] - dom->y[j-1] ); 
                double delta_z = 0.5*( dom->z[k+1] - dom->z[k-1] );                
                /// x-direction inviscid & viscous terms
                double S_x = abs( u_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_x/S_x );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_x, 2.0 )/( mu*gamma + epsilon ) );
                /// y-direction inviscid & viscous terms
                double S_y = abs( v_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_y/S_y );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_y, 2.0 )/( mu*gamma + epsilon ) );
                /// z-direction inviscid & viscous terms
                double S_z = abs( w_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_z/S_z );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_z, 2.0 )/( mu*gamma + epsilon ) );
            }
        }
    }

    /// Find minimum (global) delta_t
    double global_delta_t;
    MPI_Allreduce(&local_delta_t, &global_delta_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    /// Set new time step
    delta_t = global_delta_t;

};




////////// MAIN //////////
int main(int argc, char** argv) {

    /// Initialize MPI
    int world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    /// Start RHEA simulation
    if(world_rank == 0) cout << "RHEA: START SIMULATION" << endl;

    /// Construct flow solver RHEA
    FlowSolverRHEA flow_solver_RHEA("configuration_file.yaml");

    /// Set initial conditions
    flow_solver_RHEA.setInitialConditions();

    /// Initialize thermodynamics
    flow_solver_RHEA.initializeThermodynamics();

    /// Calculate conserved variables from primitive variables
    flow_solver_RHEA.primitiveToConservedVariables();

    /// Update previous state of conserved variables
    flow_solver_RHEA.updatePreviousStateConservedVariables();

    /// Iterate flow solver RHEA in time
    for(int t = flow_solver_RHEA.getCurrentTimeIteration(); t <= flow_solver_RHEA.getFinalTimeIteration(); t++) {

        /// Calculate time step
        flow_solver_RHEA.calculateTimeStep();
        if( ( flow_solver_RHEA.getCurrentTime() + flow_solver_RHEA.getDeltaTime() ) > flow_solver_RHEA.getFinalTime() ) {
            flow_solver_RHEA.setDeltaTime( flow_solver_RHEA.getFinalTime() - flow_solver_RHEA.getCurrentTime() );
        }

        /// Print time iteration information
        if(world_rank == 0) {
            cout << "Time iteration " << t << ":" 
                 << "t = " << flow_solver_RHEA.getCurrentTime() << " [s], "
                 << "delta_t = " << flow_solver_RHEA.getDeltaTime() << " [s]" << endl;
        }

        /// Update time
        flow_solver_RHEA.setCurrentTime( flow_solver_RHEA.getCurrentTime() + flow_solver_RHEA.getDeltaTime() );

        /// Check if simulation is completed: current_time > final_time
        if(flow_solver_RHEA.getCurrentTime() >= flow_solver_RHEA.getFinalTime() ) break;

    }



    /// Destruct flow solver RHEA ... destructor is called automatically

    /// End RHEA simulation
    if(world_rank == 0) cout << "RHEA: END SIMULATION" << endl;

    /// Finalize MPI
    MPI_Finalize();

}

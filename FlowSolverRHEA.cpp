#include "FlowSolverRHEA.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
const double epsilon = 1.0e-15;			// Small epsilon number (fixed)
const double Pi      = 2.0*asin(1.0);		// Pi number (fixed)

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

    R_specific          = 287.058;
    gamma               = 1.4;
    mu                  = rho_ref*u_tau*delta/Re_tau;
    kappa               = 0.0;
    x_0                 = 0.0;
    y_0                 = 0.0;
    z_0                 = 0.0;
    L_x                 = 4.0*Pi*delta;
    L_y                 = 2.0*delta;
    L_z                 = 4.0*Pi*delta/3.0;
    current_time        = 0.0;
    final_time          = 1.0e3;
    output_data_file    = "output_state_data.csv";
    num_grid_x          = 64;
    num_grid_y          = 64;
    num_grid_z          = 64;
    CFL                 = 0.9;
    current_time_iter   = 0;
    final_time_iter     = 1e6;
    output_iter         = 1e2;
    rk_order            = 3;
    bocos_type[_WEST_]  = _PERIODIC_;
    bocos_type[_EAST_]  = _PERIODIC_;
    bocos_type[_SOUTH_] = _PERIODIC_;
    bocos_type[_NORTH_] = _PERIODIC_;
    bocos_type[_BACK_]  = _PERIODIC_;
    bocos_type[_FRONT_] = _PERIODIC_;
    bocos_u[_WEST_]     = 0.0;
    bocos_u[_EAST_]     = 0.0;
    bocos_u[_SOUTH_]    = 0.0;
    bocos_u[_NORTH_]    = 0.0;
    bocos_u[_BACK_]     = 0.0;
    bocos_u[_FRONT_]    = 0.0;
    bocos_v[_WEST_]     = 0.0;
    bocos_v[_EAST_]     = 0.0;
    bocos_v[_SOUTH_]    = 0.0;
    bocos_v[_NORTH_]    = 0.0;
    bocos_v[_BACK_]     = 0.0;
    bocos_v[_FRONT_]    = 0.0;
    bocos_w[_WEST_]     = 0.0;
    bocos_w[_EAST_]     = 0.0;
    bocos_w[_SOUTH_]    = 0.0;
    bocos_w[_NORTH_]    = 0.0;
    bocos_w[_BACK_]     = 0.0;
    bocos_w[_FRONT_]    = 0.0;
    bocos_P[_WEST_]     = 0.0;
    bocos_P[_EAST_]     = 0.0;
    bocos_P[_SOUTH_]    = 0.0;
    bocos_P[_NORTH_]    = 0.0;
    bocos_P[_BACK_]     = 0.0;
    bocos_P[_FRONT_]    = 0.0;
    bocos_T[_WEST_]     = 0.0;
    bocos_T[_EAST_]     = 0.0;
    bocos_T[_SOUTH_]    = 0.0;
    bocos_T[_NORTH_]    = 0.0;
    bocos_T[_BACK_]     = 0.0;
    bocos_T[_FRONT_]    = 0.0;
    np_x                = 2;
    np_y                = 1;
    np_z                = 1;
    // The lines above need to be introduced via configuration (input) file !!

    /// Construct (initialize) computational domain
    mesh = new domain(L_x, L_y, L_z, x_0, y_0, z_0, RHEA_NX, RHEA_NY, RHEA_NZ);

    /// Add boundary conditions to computational domain
    mesh->updateBocos(bocos_type);

    /// Construct (initialize) communication scheme
    topo = new comm_scheme(mesh, np_x, np_y, np_z);
    //if(topo->getRank() == 0) mesh->printDomain();
    //for(int p = 0; p < np_x*np_y*np_z; p++) topo->printCommSchemeToFile(p);

    // The lines below are temporary ... will need to be removed!
    _lNx_ = topo->getlNx();
    _lNy_ = topo->getlNy();
    _lNz_ = topo->getlNz();

    /// Set parallel topology of primitive, conserved and thermodynamic variables	
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

    /// Set parallel topology of time-integration variables	
    rho_0_field.setTopology(topo);
    rhou_0_field.setTopology(topo);
    rhov_0_field.setTopology(topo);
    rhow_0_field.setTopology(topo);
    rhoE_0_field.setTopology(topo);    

    /// Set parallel topology of time-integration fluxes	
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

    /// Set parallel topology of inviscid fluxes	
    rho_inv_flux.setTopology(topo);
    rhou_inv_flux.setTopology(topo);
    rhov_inv_flux.setTopology(topo);
    rhow_inv_flux.setTopology(topo);
    rhoE_inv_flux.setTopology(topo);

    /// Set parallel topology of viscous fluxes	
    rhou_vis_flux.setTopology(topo);
    rhov_vis_flux.setTopology(topo);
    rhow_vis_flux.setTopology(topo);
    rhoE_vis_flux.setTopology(topo);

    /// Set parallel topology of source terms
    f_rhou_field.setTopology(topo);
    f_rhov_field.setTopology(topo);
    f_rhow_field.setTopology(topo);
    f_rhoE_field.setTopology(topo);

};

FlowSolverRHEA::~FlowSolverRHEA() {

    /// Free mesh, topo and bocos
    if(mesh != NULL) free(mesh);	
    if(topo != NULL) free(topo);
    free(bocos_type);	
    free(bocos_u);	
    free(bocos_v);	
    free(bocos_w);	
    free(bocos_P);	
    free(bocos_T);	

};

void FlowSolverRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// All (inner & boundary) points: u, v, w, P and T
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
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

void FlowSolverRHEA::initializeThermodynamics() {

    /// Ideal-gas model:
    /// rho = P/(R_specific*T) is density
    /// e = P/(rho*(gamma - 1)) is specific internal energy
    /// ke = (u*u + v*v + w*w)/2 is specific kinetic energy
    /// E = e + ke is total energy
    /// sos = sqrt(gamma*P/rho) is speed of sound

    /// All (inner & boundary) points: rho, e, ke, E and sos
    double e, ke;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( R_specific*T_field[I1D(i,j,k)] );
                e                     = P_field[I1D(i,j,k)]/( rho_field[I1D(i,j,k)]*( gamma - 1.0 ) );
                ke                    = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) );
                E_field[I1D(i,j,k)]   = e + ke;
                sos_field[I1D(i,j,k)] = sqrt( gamma*P_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)] );
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
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
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

void FlowSolverRHEA::conservedToPrimitiveVariables() {

    /// All (inner & boundary) points: u, v, w and E
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                u_field[I1D(i,j,k)] = rhou_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)]; 
                v_field[I1D(i,j,k)] = rhov_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)]; 
                w_field[I1D(i,j,k)] = rhow_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)]; 
                E_field[I1D(i,j,k)] = rhoE_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)]; 
            }
        }
    }

    /// Update halo values
    u_field.update();
    v_field.update();
    w_field.update();
    E_field.update();

};

void FlowSolverRHEA::calculateThermodynamicsFromPrimitiveVariables() {

    /// Ideal-gas model:
    /// c_v = R_specific/(gamma - 1) is specific heat capcity at constant volume
    /// ke = (u*u + v*v + w*w)/2 is specific kinetic energy
    /// e = E - ke is specific internal energy
    /// P = e*rho*(gamma - 1) is pressure
    /// T = e/c_v is temperature
    /// sos = sqrt(gamma*P/rho) is speed of sound

    /// Specific heat capacity at constant volume
    const double c_v = R_specific/( gamma - 1.0 );

    /// All (inner & boundary) points: ke, e, P, T and sos
    double ke, e;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                ke                    = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) ); 
                e                     = E_field[I1D(i,j,k)] - ke;
                P_field[I1D(i,j,k)]   = e*rho_field[I1D(i,j,k)]*( gamma - 1.0 ); 
                T_field[I1D(i,j,k)]   = e/c_v; 
                sos_field[I1D(i,j,k)] = sqrt( gamma*P_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    P_field.update();
    T_field.update();
    sos_field.update();

};

void FlowSolverRHEA::updateBoundaries() {

    /// General form: w_g*phi_g + w_in*phi_in = phi_b
    /// phi_g is ghost cell value
    /// phi_in is inner cell value
    /// phi_b is boundary value/flux
    /// w_g is ghost cell weight
    /// w_in is inner cell weight

    /// Declare variables & weights
    double rho_b, rhou_b, rhov_b, rhov_b, e_b, ke_b, rhoE_b, w_g, w_in;

    /// West boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_WEST_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_WEST_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGlobx[1] - mesh->getGlobx[0] );
        w_in = ( -1.0 )/( mesh->getGlobx[1] - mesh->getGlobx[0] );
    }
    rho_b  = bocos_P[_WEST_]/( R_specific*bocos_T[_WEST_] + epsilon );
    rhou_b = rho_b*bocos_u[_WEST_];
    rhov_b = rho_b*bocos_v[_WEST_];
    rhow_b = rho_b*bocos_w[_WEST_];
    e_b    = bocos_P[_WEST_]/( rho_b*( gamma - 1.0 ) + epsilon );
    ke_b   = 0.5*( pow( bocos_u[_WEST_], 2.0 ) + pow( bocos_v[_WEST_], 2.0 ) + pow( bocos_w[_WEST_], 2.0 ) );
    rhoE_b = rho_b*( e_b + ke_b );
    for(int i = topo->iter_bound[_WEST_][_INIX_]; i <= topo->iter_bound[_WEST_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_WEST_][_INIY_]; j <= topo->iter_bound[_WEST_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_WEST_][_INIZ_]; k <= topo->iter_bound[_WEST_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)]  = ( rho_b - w_in*rho_field[I1D(i+1,j,k)] )/w_g;
                rhou_field[I1D(i,j,k)] = ( rhou_b - w_in*rhou_field[I1D(i+1,j,k)] )/w_g;
                rhov_field[I1D(i,j,k)] = ( rhov_b - w_in*rhov_field[I1D(i+1,j,k)] )/w_g;
                rhow_field[I1D(i,j,k)] = ( rhow_b - w_in*rhow_field[I1D(i+1,j,k)] )/w_g;
                rhoE_field[I1D(i,j,k)] = ( rhoE_b - w_in*rhoE_field[I1D(i+1,j,k)] )/w_g;
            }
        }
    }

    /// East boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_EAST_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_EAST_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGlobx[mesh->getGNx()+1] - mesh->getGlobx[mesh->getGNx()] );
        w_in = ( -1.0 )/( mesh->getGlobx[mesh->getGNx()+1] - mesh->getGlobx[mesh->getGNx()] );
    }
    rho_b  = bocos_P[_EAST_]/( R_specific*bocos_T[_EAST_] + epsilon );
    rhou_b = rho_b*bocos_u[_EAST_];
    rhov_b = rho_b*bocos_v[_EAST_];
    rhow_b = rho_b*bocos_w[_EAST_];
    e_b    = bocos_P[_EAST_]/( rho_b*( gamma - 1.0 ) + epsilon );
    ke_b   = 0.5*( pow( bocos_u[_EAST_], 2.0 ) + pow( bocos_v[_EAST_], 2.0 ) + pow( bocos_w[_EAST_], 2.0 ) );
    rhoE_b = rho_b*( e_b + ke_b );
    for(int i = topo->iter_bound[_EAST_][_INIX_]; i <= topo->iter_bound[_EAST_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_EAST_][_INIY_]; j <= topo->iter_bound[_EAST_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_EAST_][_INIZ_]; k <= topo->iter_bound[_EAST_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)]  = ( rho_b - w_in*rho_field[I1D(i-1,j,k)] )/w_g;
                rhou_field[I1D(i,j,k)] = ( rhou_b - w_in*rhou_field[I1D(i-1,j,k)] )/w_g;
                rhov_field[I1D(i,j,k)] = ( rhov_b - w_in*rhov_field[I1D(i-1,j,k)] )/w_g;
                rhow_field[I1D(i,j,k)] = ( rhow_b - w_in*rhow_field[I1D(i-1,j,k)] )/w_g;
                rhoE_field[I1D(i,j,k)] = ( rhoE_b - w_in*rhoE_field[I1D(i-1,j,k)] )/w_g;
            }
        }
    }

    /// South boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_SOUTH_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_SOUTH_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGloby[1] - mesh->getGloby[0] );
        w_in = ( -1.0 )/( mesh->getGloby[1] - mesh->getGloby[0] );
    }
    rho_b  = bocos_P[_SOUTH_]/( R_specific*bocos_T[_SOUTH_] + epsilon );
    rhou_b = rho_b*bocos_u[_SOUTH_];
    rhov_b = rho_b*bocos_v[_SOUTH_];
    rhow_b = rho_b*bocos_w[_SOUTH_];
    e_b    = bocos_P[_SOUTH_]/( rho_b*( gamma - 1.0 ) + epsilon );
    ke_b   = 0.5*( pow( bocos_u[_SOUTH_], 2.0 ) + pow( bocos_v[_SOUTH_], 2.0 ) + pow( bocos_w[_SOUTH_], 2.0 ) );
    rhoE_b = rho_b*( e_b + ke_b );
    for(int i = topo->iter_bound[_SOUTH_][_INIX_]; i <= topo->iter_bound[_SOUTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_SOUTH_][_INIY_]; j <= topo->iter_bound[_SOUTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_SOUTH_][_INIZ_]; k <= topo->iter_bound[_SOUTH_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)]  = ( rho_b - w_in*rho_field[I1D(i,j+1,k)] )/w_g;
                rhou_field[I1D(i,j,k)] = ( rhou_b - w_in*rhou_field[I1D(i,j+1,k)] )/w_g;
                rhov_field[I1D(i,j,k)] = ( rhov_b - w_in*rhov_field[I1D(i,j+1,k)] )/w_g;
                rhow_field[I1D(i,j,k)] = ( rhow_b - w_in*rhow_field[I1D(i,j+1,k)] )/w_g;
                rhoE_field[I1D(i,j,k)] = ( rhoE_b - w_in*rhoE_field[I1D(i,j+1,k)] )/w_g;
            }
        }
    }

    /// North boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_NORTH_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_NORTH_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGloby[mesh->getGNy()+1] - mesh->getGloby[mesh->getGNy()] );
        w_in = ( -1.0 )/( mesh->getGloby[mesh->getGNy()+1] - mesh->getGloby[mesh->getGNy()] );
    }
    rho_b  = bocos_P[_NORTH_]/( R_specific*bocos_T[_NORTH_] + epsilon );
    rhou_b = rho_b*bocos_u[_NORTH_];
    rhov_b = rho_b*bocos_v[_NORTH_];
    rhow_b = rho_b*bocos_w[_NORTH_];
    e_b    = bocos_P[_NORTH_]/( rho_b*( gamma - 1.0 ) + epsilon );
    ke_b   = 0.5*( pow( bocos_u[_NORTH_], 2.0 ) + pow( bocos_v[_NORTH_], 2.0 ) + pow( bocos_w[_NORTH_], 2.0 ) );
    rhoE_b = rho_b*( e_b + ke_b );
    for(int i = topo->iter_bound[_NORTH_][_INIX_]; i <= topo->iter_bound[_NORTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_NORTH_][_INIY_]; j <= topo->iter_bound[_NORTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_NORTH_][_INIZ_]; k <= topo->iter_bound[_NORTH_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)]  = ( rho_b - w_in*rho_field[I1D(i,j-1,k)] )/w_g;
                rhou_field[I1D(i,j,k)] = ( rhou_b - w_in*rhou_field[I1D(i,j-1,k)] )/w_g;
                rhov_field[I1D(i,j,k)] = ( rhov_b - w_in*rhov_field[I1D(i,j-1,k)] )/w_g;
                rhow_field[I1D(i,j,k)] = ( rhow_b - w_in*rhow_field[I1D(i,j-1,k)] )/w_g;
                rhoE_field[I1D(i,j,k)] = ( rhoE_b - w_in*rhoE_field[I1D(i,j-1,k)] )/w_g;
            }
        }
    }

    /// Back boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_BACK_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_BACK_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGlobz[1] - mesh->getGlobz[0] );
        w_in = ( -1.0 )/( mesh->getGlobz[1] - mesh->getGlobz[0] );
    }
    rho_b  = bocos_P[_BACK_]/( R_specific*bocos_T[_BACK_] + epsilon );
    rhou_b = rho_b*bocos_u[_BACK_];
    rhov_b = rho_b*bocos_v[_BACK_];
    rhow_b = rho_b*bocos_w[_BACK_];
    e_b    = bocos_P[_BACK_]/( rho_b*( gamma - 1.0 ) + epsilon );
    ke_b   = 0.5*( pow( bocos_u[_BACK_], 2.0 ) + pow( bocos_v[_BACK_], 2.0 ) + pow( bocos_w[_BACK_], 2.0 ) );
    rhoE_b = rho_b*( e_b + ke_b );
    for(int i = topo->iter_bound[_BACK_][_INIX_]; i <= topo->iter_bound[_BACK_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_BACK_][_INIY_]; j <= topo->iter_bound[_BACK_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_BACK_][_INIZ_]; k <= topo->iter_bound[_BACK_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)]  = ( rho_b - w_in*rho_field[I1D(i,j,k+1)] )/w_g;
                rhou_field[I1D(i,j,k)] = ( rhou_b - w_in*rhou_field[I1D(i,j,k+1)] )/w_g;
                rhov_field[I1D(i,j,k)] = ( rhov_b - w_in*rhov_field[I1D(i,j,k+1)] )/w_g;
                rhow_field[I1D(i,j,k)] = ( rhow_b - w_in*rhow_field[I1D(i,j,k+1)] )/w_g;
                rhoE_field[I1D(i,j,k)] = ( rhoE_b - w_in*rhoE_field[I1D(i,j,k+1)] )/w_g;
            }
        }
    }

    /// Front boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_FRONT_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_FRONT_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGlobz[mesh->getGNz()+1] - mesh->getGlobz[mesh->getGNz()] );
        w_in = ( -1.0 )/( mesh->getGlobz[mesh->getGNz()+1] - mesh->getGlobz[mesh->getGNz()] );
    }
    rho_b  = bocos_P[_FRONT_]/( R_specific*bocos_T[_FRONT_] + epsilon );
    rhou_b = rho_b*bocos_u[_FRONT_];
    rhov_b = rho_b*bocos_v[_FRONT_];
    rhow_b = rho_b*bocos_w[_FRONT_];
    e_b    = bocos_P[_FRONT_]/( rho_b*( gamma - 1.0 ) + epsilon );
    ke_b   = 0.5*( pow( bocos_u[_FRONT_], 2.0 ) + pow( bocos_v[_FRONT_], 2.0 ) + pow( bocos_w[_FRONT_], 2.0 ) );
    rhoE_b = rho_b*( e_b + ke_b );
    for(int i = topo->iter_bound[_FRONT_][_INIX_]; i <= topo->iter_bound[_FRONT_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_FRONT_][_INIY_]; j <= topo->iter_bound[_FRONT_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_FRONT_][_INIZ_]; k <= topo->iter_bound[_FRONT_][_ENDZ_]; k++) {
                rho_field[I1D(i,j,k)]  = ( rho_b - w_in*rho_field[I1D(i,j,k-1)] )/w_g;
                rhou_field[I1D(i,j,k)] = ( rhou_b - w_in*rhou_field[I1D(i,j,k-1)] )/w_g;
                rhov_field[I1D(i,j,k)] = ( rhov_b - w_in*rhov_field[I1D(i,j,k-1)] )/w_g;
                rhow_field[I1D(i,j,k)] = ( rhow_b - w_in*rhow_field[I1D(i,j,k-1)] )/w_g;
                rhoE_field[I1D(i,j,k)] = ( rhoE_b - w_in*rhoE_field[I1D(i,j,k-1)] )/w_g;
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

void FlowSolverRHEA::updatePreviousStateConservedVariables() {

    /// All (inner & boundary) points: rho_0, rhou_0 rhov_0, rhow_0 and rhoE_0
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
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
    double delta_x, delta_y, delta_z;
    double S_x, S_y, S_z;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = 0.5*( mesh->x[i+1] - mesh->x[i-1] ); 
                delta_y = 0.5*( mesh->y[j+1] - mesh->y[j-1] ); 
                delta_z = 0.5*( mesh->z[k+1] - mesh->z[k-1] );                
                /// x-direction inviscid & viscous terms
                S_x           = abs( u_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_x/S_x );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_x, 2.0 )/( mu*gamma + epsilon ) );
                /// y-direction inviscid & viscous terms
                S_y           = abs( v_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_y/S_y );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_y, 2.0 )/( mu*gamma + epsilon ) );
                /// z-direction inviscid & viscous terms
                S_z           = abs( w_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
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

void FlowSolverRHEA::calculateSourceTerms() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                f_rhou_field[I1D(i,j,k)] = 0.0;
                f_rhov_field[I1D(i,j,k)] = 0.0;
                f_rhow_field[I1D(i,j,k)] = 0.0;
                f_rhoE_field[I1D(i,j,k)] = 0.0;
                //f_rhoE_field[I1D(i,j,k)] = ( -1.0 )*( f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)] )
            }
        }
    }

    /// Update halo values
    f_rhou_field.update();
    f_rhov_field.update();
    f_rhow_field.update();
    f_rhoE_field.update();

};

void FlowSolverRHEA::calculateWavesSpeed(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, double &S_L, double &S_R) {

    /// HLLC approximate Riemann solver:
    /// E. F. Toro.
    /// Riemann solvers and numerical methods for fluid dynamics.
    /// Springer, 2009.

    double rho_bar = 0.5*( rho_L + rho_R );
    double a_bar   = 0.5*( a_L + a_R );
    double P_pvrs  = 0.5*( P_L + P_R ) - 0.5*( u_R - u_L )*rho_bar*a_bar;
    double P_star  = max( 0.0, P_pvrs );
    double q_L     = 1.0;
    if(P_star > P_L) q_L = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_L ) - 1.0 ) );
    double q_R     = 1.0;
    if(P_star > P_R) q_R = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_R ) - 1.0 ) );
    S_L = u_L - a_L*q_L;
    S_R = u_R + a_R*q_R;

};

double FlowSolverRHEA::calculateHllcFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// HLLC approximate Riemann solver:
    /// E. F. Toro.
    /// Riemann solvers and numerical methods for fluid dynamics.
    /// Springer, 2009.

    double S_L, S_R;
    this->calculateWavesSpeed( rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R, S_L, S_R );
    double S_star   = ( P_R - P_L + rho_L*u_L*( S_L - u_L ) - rho_R*u_R*( S_R - u_R ) )/( rho_L*( S_L - u_L ) - rho_R*( S_R - u_R ) );
    double U_star_L = rho_L*( ( S_L - u_L )/( S_L - S_star ) );
    double U_star_R = rho_R*( ( S_R - u_R )/( S_R - S_star ) );
    if(var_type == 0) {
        U_star_L *= 1.0;
        U_star_R *= 1.0;       
    } else if(var_type == 1) {
        U_star_L *= S_star;
        U_star_R *= S_star;
    } else if(var_type == 2) {
        U_star_L *= v_L;
        U_star_R *= v_R;
    } else if(var_type == 3) {
        U_star_L *= w_L;
        U_star_R *= w_R;
    } else if(var_type == 4) {
        U_star_L *= ( E_L + ( S_star - u_L )*( S_star + P_L/( rho_L*( S_L - u_L ) ) ) );
        U_star_R *= ( E_R + ( S_star - u_R )*( S_star + P_R/( rho_R*( S_R - u_R ) ) ) );
    }
    double F_star_L = F_L + S_L*( U_star_L - U_L );
    double F_star_R = F_R + S_R*( U_star_R - U_R );
    double F = 0.0;
    if(0.0 <= S_L) {
        F = F_L;
    } else if(( S_L <= 0.0 ) && ( 0.0 <= S_star )) {
        F = F_star_L;
    } else if(( S_star <= 0.0 ) && ( 0.0 <= S_R )) {
        F = F_star_R;
    } else if(0.0 >= S_R) {
        F = F_R;
    }

    return( F );

};

void FlowSolverRHEA::calculateInviscidFluxes() {

    /// First-order Godunov-type unsplit method for Euler equations:
    /// E. F. Toro.
    /// Riemann solvers and numerical methods for fluid dynamics.
    /// Springer, 2009.

    /// Inner points: rho, rhou, rhov, rhow and rhoE
    int index_L, index_R, var_type;
    double delta_x, delta_y, delta_z;
    double rho_L, u_L, v_L, w_L, E_L, P_L, a_L;
    double rho_R, u_R, v_R, w_R, E_R, P_R, a_R;
    double rho_F_L, rho_U_L, rho_F_R, rho_U_R, rho_F_p, rho_F_m;
    double rhou_F_L, rhou_U_L, rhou_F_R, rhou_U_R, rhou_F_p, rhou_F_m;
    double rhov_F_L, rhov_U_L, rhov_F_R, rhov_U_R, rhov_F_p, rhov_F_m;
    double rhow_F_L, rhow_U_L, rhow_F_R, rhow_U_R, rhow_F_p, rhow_F_m;
    double rhoE_F_L, rhoE_U_L, rhoE_F_R, rhoE_U_R, rhoE_F_p, rhoE_F_m;
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
                a_L     = sos_field[I1D(index_L,j,k)]; a_R     = sos_field[I1D(index_R,j,k)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)];
                rho_F_R  = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)];
                rho_U_L  = rho_field[I1D(index_L,j,k)];
                rho_U_R  = rho_field[I1D(index_R,j,k)];
                rho_F_p = this->calculateHllcFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 1;
                rhou_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)] + P_field[I1D(index_L,j,k)];
                rhou_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)] + P_field[I1D(index_R,j,k)];
                rhou_U_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)];
                rhou_U_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)];
                rhou_F_p = this->calculateHllcFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*v_field[I1D(index_L,j,k)];
                rhov_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*v_field[I1D(index_R,j,k)];
                rhov_U_L = rho_field[I1D(index_L,j,k)]*v_field[I1D(index_L,j,k)];
                rhov_U_R = rho_field[I1D(index_R,j,k)]*v_field[I1D(index_R,j,k)];
                rhov_F_p = this->calculateHllcFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*w_field[I1D(index_L,j,k)];
                rhow_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*w_field[I1D(index_R,j,k)];
                rhow_U_L = rho_field[I1D(index_L,j,k)]*w_field[I1D(index_L,j,k)];
                rhow_U_R = rho_field[I1D(index_R,j,k)]*w_field[I1D(index_R,j,k)];
                rhow_F_p = this->calculateHllcFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*E_field[I1D(index_L,j,k)] + u_field[I1D(index_L,j,k)]*P_field[I1D(index_L,j,k)];
                rhoE_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*E_field[I1D(index_R,j,k)] + u_field[I1D(index_R,j,k)]*P_field[I1D(index_R,j,k)];
                rhoE_U_L = rho_field[I1D(index_L,j,k)]*E_field[I1D(index_L,j,k)];
                rhoE_U_R = rho_field[I1D(index_R,j,k)]*E_field[I1D(index_R,j,k)];
                rhoE_F_p = this->calculateHllcFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// x-direction i-1/2
                index_L = i - 1;                       index_R = i;
                rho_L   = rho_field[I1D(index_L,j,k)]; rho_R   = rho_field[I1D(index_R,j,k)];
                u_L     = u_field[I1D(index_L,j,k)];   u_R     = u_field[I1D(index_R,j,k)];
                v_L     = v_field[I1D(index_L,j,k)];   v_R     = v_field[I1D(index_R,j,k)];
                w_L     = w_field[I1D(index_L,j,k)];   w_R     = w_field[I1D(index_R,j,k)];
                E_L     = E_field[I1D(index_L,j,k)];   E_R     = E_field[I1D(index_R,j,k)];
                P_L     = P_field[I1D(index_L,j,k)];   P_R     = P_field[I1D(index_R,j,k)];
                a_L     = sos_field[I1D(index_L,j,k)]; a_R     = sos_field[I1D(index_R,j,k)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)];
                rho_F_R  = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)];
                rho_U_L  = rho_field[I1D(index_L,j,k)];
                rho_U_R  = rho_field[I1D(index_R,j,k)];
                rho_F_m = this->calculateHllcFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 1;
                rhou_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)] + P_field[I1D(index_L,j,k)];
                rhou_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)] + P_field[I1D(index_R,j,k)];
                rhou_U_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)];
                rhou_U_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)];
                rhou_F_m = this->calculateHllcFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*v_field[I1D(index_L,j,k)];
                rhov_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*v_field[I1D(index_R,j,k)];
                rhov_U_L = rho_field[I1D(index_L,j,k)]*v_field[I1D(index_L,j,k)];
                rhov_U_R = rho_field[I1D(index_R,j,k)]*v_field[I1D(index_R,j,k)];
                rhov_F_m = this->calculateHllcFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*w_field[I1D(index_L,j,k)];
                rhow_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*w_field[I1D(index_R,j,k)];
                rhow_U_L = rho_field[I1D(index_L,j,k)]*w_field[I1D(index_L,j,k)];
                rhow_U_R = rho_field[I1D(index_R,j,k)]*w_field[I1D(index_R,j,k)];
                rhow_F_m = this->calculateHllcFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_field[I1D(index_L,j,k)]*u_field[I1D(index_L,j,k)]*E_field[I1D(index_L,j,k)] + u_field[I1D(index_L,j,k)]*P_field[I1D(index_L,j,k)];
                rhoE_F_R = rho_field[I1D(index_R,j,k)]*u_field[I1D(index_R,j,k)]*E_field[I1D(index_R,j,k)] + u_field[I1D(index_R,j,k)]*P_field[I1D(index_R,j,k)];
                rhoE_U_L = rho_field[I1D(index_L,j,k)]*E_field[I1D(index_L,j,k)];
                rhoE_U_R = rho_field[I1D(index_R,j,k)]*E_field[I1D(index_R,j,k)];
                rhoE_F_m = this->calculateHllcFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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
                a_L     = sos_field[I1D(i,index_L,k)]; a_R     = sos_field[I1D(i,index_R,k)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)];
                rho_F_R  = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)];
                rho_U_L  = rho_field[I1D(i,index_L,k)];
                rho_U_R  = rho_field[I1D(i,index_R,k)];
                rho_F_p = this->calculateHllcFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 2;
                rhou_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*u_field[I1D(i,index_L,k)];
                rhou_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*u_field[I1D(i,index_R,k)];
                rhou_U_L = rho_field[I1D(i,index_L,k)]*u_field[I1D(i,index_L,k)];
                rhou_U_R = rho_field[I1D(i,index_R,k)]*u_field[I1D(i,index_R,k)];
                rhou_F_p = this->calculateHllcFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 1;
                rhov_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)] + P_field[I1D(i,index_L,k)];
                rhov_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)] + P_field[I1D(i,index_R,k)];
                rhov_U_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)];
                rhov_U_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)];
                rhov_F_p = this->calculateHllcFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*w_field[I1D(i,index_L,k)];
                rhow_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*w_field[I1D(i,index_R,k)];
                rhow_U_L = rho_field[I1D(i,index_L,k)]*w_field[I1D(i,index_L,k)];
                rhow_U_R = rho_field[I1D(i,index_R,k)]*w_field[I1D(i,index_R,k)];
                rhow_F_p = this->calculateHllcFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*E_field[I1D(i,index_L,k)] + v_field[I1D(i,index_L,k)]*P_field[I1D(i,index_L,k)];
                rhoE_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*E_field[I1D(i,index_R,k)] + v_field[I1D(i,index_R,k)]*P_field[I1D(i,index_R,k)];
                rhoE_U_L = rho_field[I1D(i,index_L,k)]*E_field[I1D(i,index_L,k)];
                rhoE_U_R = rho_field[I1D(i,index_R,k)]*E_field[I1D(i,index_R,k)];
                rhoE_F_p = this->calculateHllcFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// y-direction j-1/2
                index_L = j - 1;                       index_R = j;
                rho_L   = rho_field[I1D(i,index_L,k)]; rho_R   = rho_field[I1D(i,index_R,k)];
                u_L     = v_field[I1D(i,index_L,k)];   u_R     = v_field[I1D(i,index_R,k)];
                v_L     = u_field[I1D(i,index_L,k)];   v_R     = u_field[I1D(i,index_R,k)];
                w_L     = w_field[I1D(i,index_L,k)];   w_R     = w_field[I1D(i,index_R,k)];
                E_L     = E_field[I1D(i,index_L,k)];   E_R     = E_field[I1D(i,index_R,k)];
                P_L     = P_field[I1D(i,index_L,k)];   P_R     = P_field[I1D(i,index_R,k)];
                a_L     = sos_field[I1D(i,index_L,k)]; a_R     = sos_field[I1D(i,index_R,k)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)];
                rho_F_R  = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)];
                rho_U_L  = rho_field[I1D(i,index_L,k)];
                rho_U_R  = rho_field[I1D(i,index_R,k)];
                rho_F_m = this->calculateHllcFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 2;
                rhou_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*u_field[I1D(i,index_L,k)];
                rhou_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*u_field[I1D(i,index_R,k)];
                rhou_U_L = rho_field[I1D(i,index_L,k)]*u_field[I1D(i,index_L,k)];
                rhou_U_R = rho_field[I1D(i,index_R,k)]*u_field[I1D(i,index_R,k)];
                rhou_F_m = this->calculateHllcFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 1;
                rhov_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)] + P_field[I1D(i,index_L,k)];
                rhov_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)] + P_field[I1D(i,index_R,k)];
                rhov_U_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)];
                rhov_U_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)];
                rhov_F_m = this->calculateHllcFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*w_field[I1D(i,index_L,k)];
                rhow_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*w_field[I1D(i,index_R,k)];
                rhow_U_L = rho_field[I1D(i,index_L,k)]*w_field[I1D(i,index_L,k)];
                rhow_U_R = rho_field[I1D(i,index_R,k)]*w_field[I1D(i,index_R,k)];
                rhow_F_m = this->calculateHllcFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_field[I1D(i,index_L,k)]*v_field[I1D(i,index_L,k)]*E_field[I1D(i,index_L,k)] + v_field[I1D(i,index_L,k)]*P_field[I1D(i,index_L,k)];
                rhoE_F_R = rho_field[I1D(i,index_R,k)]*v_field[I1D(i,index_R,k)]*E_field[I1D(i,index_R,k)] + v_field[I1D(i,index_R,k)]*P_field[I1D(i,index_R,k)];
                rhoE_U_L = rho_field[I1D(i,index_L,k)]*E_field[I1D(i,index_L,k)];
                rhoE_U_R = rho_field[I1D(i,index_R,k)]*E_field[I1D(i,index_R,k)];
                rhoE_F_m = this->calculateHllcFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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
                a_L     = sos_field[I1D(i,j,index_L)]; a_R     = sos_field[I1D(i,j,index_R)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)];
                rho_F_R  = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)];
                rho_U_L  = rho_field[I1D(i,j,index_L)];
                rho_U_R  = rho_field[I1D(i,j,index_R)];
                rho_F_p = this->calculateHllcFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 3;
                rhou_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*u_field[I1D(i,j,index_L)];
                rhou_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*u_field[I1D(i,j,index_R)];
                rhou_U_L = rho_field[I1D(i,j,index_L)]*u_field[I1D(i,j,index_L)];
                rhou_U_R = rho_field[I1D(i,j,index_R)]*u_field[I1D(i,j,index_R)];
                rhou_F_p = this->calculateHllcFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*v_field[I1D(i,j,index_L)];
                rhov_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*v_field[I1D(i,j,index_R)];
                rhov_U_L = rho_field[I1D(i,j,index_L)]*v_field[I1D(i,j,index_L)];
                rhov_U_R = rho_field[I1D(i,j,index_R)]*v_field[I1D(i,j,index_R)];
                rhov_F_p = this->calculateHllcFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 1;
                rhow_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)] + P_field[I1D(i,j,index_L)];
                rhow_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)] + P_field[I1D(i,j,index_R)];
                rhow_U_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)];
                rhow_U_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)];
                rhow_F_p = this->calculateHllcFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*E_field[I1D(i,j,index_L)] + w_field[I1D(i,j,index_L)]*P_field[I1D(i,j,index_L)];
                rhoE_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*E_field[I1D(i,j,index_R)] + w_field[I1D(i,j,index_R)]*P_field[I1D(i,j,index_R)];
                rhoE_U_L = rho_field[I1D(i,j,index_L)]*E_field[I1D(i,j,index_L)];
                rhoE_U_R = rho_field[I1D(i,j,index_R)]*E_field[I1D(i,j,index_R)];
                rhoE_F_p = this->calculateHllcFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// z-direction k-1/2
                index_L = k - 1;                       index_R = k;
                rho_L   = rho_field[I1D(i,j,index_L)]; rho_R   = rho_field[I1D(i,j,index_R)];
                u_L     = w_field[I1D(i,j,index_L)];   u_R     = w_field[I1D(i,j,index_R)];
                v_L     = v_field[I1D(i,j,index_L)];   v_R     = v_field[I1D(i,j,index_R)];
                w_L     = u_field[I1D(i,j,index_L)];   w_R     = u_field[I1D(i,j,index_R)];
                E_L     = E_field[I1D(i,j,index_L)];   E_R     = E_field[I1D(i,j,index_R)];
                P_L     = P_field[I1D(i,j,index_L)];   P_R     = P_field[I1D(i,j,index_R)];
                a_L     = sos_field[I1D(i,j,index_L)]; a_R     = sos_field[I1D(i,j,index_R)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)];
                rho_F_R  = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)];
                rho_U_L  = rho_field[I1D(i,j,index_L)];
                rho_U_R  = rho_field[I1D(i,j,index_R)];
                rho_F_m = this->calculateHllcFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 3;
                rhou_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*u_field[I1D(i,j,index_L)];
                rhou_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*u_field[I1D(i,j,index_R)];
                rhou_U_L = rho_field[I1D(i,j,index_L)]*u_field[I1D(i,j,index_L)];
                rhou_U_R = rho_field[I1D(i,j,index_R)]*u_field[I1D(i,j,index_R)];
                rhou_F_m = this->calculateHllcFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*v_field[I1D(i,j,index_L)];
                rhov_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*v_field[I1D(i,j,index_R)];
                rhov_U_L = rho_field[I1D(i,j,index_L)]*v_field[I1D(i,j,index_L)];
                rhov_U_R = rho_field[I1D(i,j,index_R)]*v_field[I1D(i,j,index_R)];
                rhov_F_m = this->calculateHllcFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 1;
                rhow_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)] + P_field[I1D(i,j,index_L)];
                rhow_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)] + P_field[I1D(i,j,index_R)];
                rhow_U_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)];
                rhow_U_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)];
                rhow_F_m = this->calculateHllcFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_field[I1D(i,j,index_L)]*w_field[I1D(i,j,index_L)]*E_field[I1D(i,j,index_L)] + w_field[I1D(i,j,index_L)]*P_field[I1D(i,j,index_L)];
                rhoE_F_R = rho_field[I1D(i,j,index_R)]*w_field[I1D(i,j,index_R)]*E_field[I1D(i,j,index_R)] + w_field[I1D(i,j,index_R)]*P_field[I1D(i,j,index_R)];
                rhoE_U_L = rho_field[I1D(i,j,index_L)]*E_field[I1D(i,j,index_L)];
                rhoE_U_R = rho_field[I1D(i,j,index_R)]*E_field[I1D(i,j,index_R)];
                rhoE_F_m = this->calculateHllcFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// Fluxes z-direction
                rho_inv_flux[I1D(i,j,k)]  += ( rho_F_p - rho_F_m )/delta_z;
                rhou_inv_flux[I1D(i,j,k)] += ( rhou_F_p - rhou_F_m )/delta_z;
                rhov_inv_flux[I1D(i,j,k)] += ( rhov_F_p - rhov_F_m )/delta_z;
                rhow_inv_flux[I1D(i,j,k)] += ( rhow_F_p - rhow_F_m )/delta_z;
                rhoE_inv_flux[I1D(i,j,k)] += ( rhoE_F_p - rhoE_F_m )/delta_z;
            }
        }
    }

    /// Update halo values
    rho_inv_flux.update();
    rhou_inv_flux.update();
    rhov_inv_flux.update();
    rhow_inv_flux.update();
    rhoE_inv_flux.update();

};

void FlowSolverRHEA::calculateViscousFluxes() {

    /// Second-order central finite differences for derivatives:
    /// P. Moin.
    /// Fundamentals of engineering numerical analysis.
    /// Cambridge University Press, 2010.

    /// Inner points: rhou, rhov, rhow and rhoE
    double delta_x, delta_y, delta_z;
    double div_tau_xx, div_tau_xy, div_tau_xz, div_tau_yx, div_tau_yy, div_tau_yz, div_tau_zx, div_tau_zy, div_tau_zz;
    double div_q, vel_p, vel_m, div_tauuvw_x, div_tauuvw_y, div_tauuvw_z, div_tauuvw, f_rhouvw; 
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = 0.5*( mesh->x[i+1] - mesh->x[i-1] ); 
                delta_y = 0.5*( mesh->y[j+1] - mesh->y[j-1] ); 
                delta_z = 0.5*( mesh->z[k+1] - mesh->z[k-1] );
                /// Divergence of tau tensor terms
                div_tau_xx  = ( 2.0/delta_x )*( 2.0*mu*( ( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                       - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xx -= ( 2.0/delta_x )*( ( 2.0/3.0 )*mu*( ( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                               - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xx -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( ( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i+1,j-1,k)] )/delta_y )
                                                               - ( ( v_field[I1D(i-1,j+1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                div_tau_xx -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( ( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i+1,j,k-1)] )/delta_z )
                                                               - ( ( w_field[I1D(i-1,j,k+1)] - w_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                div_tau_xy  = ( 2.0/delta_x )*( mu*( ( ( v_field[I1D(i+1,j,k)] - v_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) ) 
                                                   - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xy += ( 1.0/delta_x )*( mu*( ( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i+1,j-1,k)] )/delta_y )
                                                   - ( ( u_field[I1D(i-1,j+1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                div_tau_xz  = ( 2.0/delta_x )*( mu*( ( ( w_field[I1D(i+1,j,k)] - w_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                   - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xz += ( 1.0/delta_x )*( mu*( ( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i+1,j,k-1)] )/delta_z )
                                                   - ( ( u_field[I1D(i-1,j,k+1)] - u_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                div_tau_yx  = ( 2.0/delta_y )*( mu*( ( ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                   - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yx += ( 1.0/delta_y )*( mu*( ( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i-1,j+1,k)] )/delta_x )
                                                   - ( ( v_field[I1D(i+1,j-1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                div_tau_yy  = ( 2.0/delta_y )*( 2.0*mu*( ( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                       - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yy -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( ( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i-1,j+1,k)] )/delta_x )
                                                               - ( ( u_field[I1D(i+1,j-1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                div_tau_yy -= ( 2.0/delta_y )*( ( 2.0/3.0 )*mu*( ( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                               - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yy -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( ( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j+1,k-1)] )/delta_z )
                                                               - ( ( w_field[I1D(i,j-1,k+1)] - w_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                div_tau_yz  = ( 2.0/delta_y )*( mu*( ( ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                   - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yz += ( 1.0/delta_y )*( mu*( ( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j+1,k-1)] )/delta_z )
                                                   - ( ( v_field[I1D(i,j-1,k+1)] - v_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                div_tau_zx  = ( 2.0/delta_z )*( mu*( ( ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                   - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                div_tau_zx += ( 1.0/delta_z )*( mu*( ( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i-1,j,k+1)] )/delta_x )
                                                   - ( ( w_field[I1D(i+1,j,k-1)] - w_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                div_tau_zy  = ( 2.0/delta_z )*( mu*( ( ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                   - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                div_tau_zy += ( 1.0/delta_z )*( mu*( ( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j-1,k+1)] )/delta_y )
                                                   - ( ( w_field[I1D(i,j+1,k-1)] - w_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                div_tau_zz  = ( 2.0/delta_z )*( 2.0*mu*( ( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                       - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                div_tau_zz -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( ( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i-1,j,k+1)] )/delta_x )
                                                               - ( ( u_field[I1D(i+1,j,k-1)] - u_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                div_tau_zz -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( ( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j-1,k+1)] )/delta_y )
                                                               - ( ( v_field[I1D(i,j+1,k-1)] - v_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                div_tau_zz -= ( 2.0/delta_z )*( ( 2.0/3.0 )*mu*( ( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                               - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                /// Fourier term
                div_q = ( -1.0 )*( ( 2.0/delta_x )*( kappa*( ( ( T_field[I1D(i+1,j,k)] - T_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                           - ( ( T_field[I1D(i,j,k)] - T_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) )
                                 + ( 2.0/delta_y )*( kappa*( ( ( T_field[I1D(i,j+1,k)] - T_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                           - ( ( T_field[I1D(i,j,k)] - T_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) )
                                 + ( 2.0/delta_z )*( kappa*( ( ( T_field[I1D(i,j,k+1)] - T_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                           - ( ( T_field[I1D(i,j,k)] - T_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) ) );
                /// Divergence of tau*velocity terms
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x  = ( 2.0/delta_x )*( 2.0*mu*( vel_p*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                         - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x -= ( 2.0/delta_x )*( ( 2.0/3.0 )*mu*( vel_p*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                 - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i+1,j,k)] ); vel_m = 0.5*( u_field[I1D(i-1,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( vel_p*( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i+1,j-1,k)] )/delta_y )
                                                                 - vel_m*( ( v_field[I1D(i-1,j+1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i+1,j,k)] ); vel_m = 0.5*( u_field[I1D(i-1,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu*( vel_p*( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i+1,j,k-1)] )/delta_z )
                                                                 - vel_m*( ( w_field[I1D(i-1,j,k+1)] - w_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( v_field[I1D(i+1,j,k)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 2.0/delta_x )*( mu*( vel_p*( ( v_field[I1D(i+1,j,k)] - v_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                     - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i+1,j,k)] + v_field[I1D(i+1,j,k)] ); vel_m = 0.5*( v_field[I1D(i-1,j,k)] + v_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 1.0/delta_x )*( mu*( vel_p*( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i+1,j-1,k)] )/delta_y )
                                                     - vel_m*( ( u_field[I1D(i-1,j+1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                vel_p = 0.5*( w_field[I1D(i+1,j,k)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 2.0/delta_x )*( mu*( vel_p*( ( w_field[I1D(i+1,j,k)] - w_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                     - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( w_field[I1D(i+1,j,k)] + w_field[I1D(i+1,j,k)] ); vel_m = 0.5*( w_field[I1D(i-1,j,k)] + w_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 1.0/delta_x )*( mu*( vel_p*( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i+1,j,k-1)] )/delta_z )
                                                     - vel_m*( ( u_field[I1D(i-1,j,k+1)] - u_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( u_field[I1D(i,j+1,k)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i,j-1,k)] );
                div_tauuvw_y  = ( 2.0/delta_y )*( mu*( vel_p*( ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                     - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i,j+1,k)] + u_field[I1D(i,j+1,k)] ); vel_m = 0.5*( u_field[I1D(i,j-1,k)] + u_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 1.0/delta_y )*( mu*( vel_p*( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i-1,j+1,k)] )/delta_x )
                                                     - vel_m*( ( v_field[I1D(i+1,j-1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 2.0/delta_y )*( 2.0*mu*( vel_p*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                         - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j+1,k)] ); vel_m = 0.5*( v_field[I1D(i,j-1,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( vel_p*( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i-1,j+1,k)] )/delta_x )
                                                                 - vel_m*( ( u_field[I1D(i+1,j-1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y -= ( 2.0/delta_y )*( ( 2.0/3.0 )*mu*( vel_p*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                 - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j+1,k)] ); vel_m = 0.5*( v_field[I1D(i,j-1,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu*( vel_p*( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j+1,k-1)] )/delta_z )
                                                                 - vel_m*( ( w_field[I1D(i,j-1,k+1)] - w_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( w_field[I1D(i,j+1,k)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 2.0/delta_y )*( mu*( vel_p*( ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                     - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( w_field[I1D(i,j+1,k)] + w_field[I1D(i,j+1,k)] ); vel_m = 0.5*( w_field[I1D(i,j-1,k)] + w_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 1.0/delta_y )*( mu*( vel_p*( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j+1,k-1)] )/delta_z )
                                                     - vel_m*( ( v_field[I1D(i,j-1,k+1)] - v_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( u_field[I1D(i,j,k+1)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i,j,k-1)] );
                div_tauuvw_z  = ( 2.0/delta_z )*( mu*( vel_p*( ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                     - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i,j,k+1)] + u_field[I1D(i,j,k+1)] ); vel_m = 0.5*( u_field[I1D(i,j,k-1)] + u_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 1.0/delta_z )*( mu*( vel_p*( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i-1,j,k+1)] )/delta_x )
                                                     - vel_m*( ( w_field[I1D(i+1,j,k-1)] - w_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                vel_p = 0.5*( v_field[I1D(i,j,k+1)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 2.0/delta_z )*( mu*( vel_p*( ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                     - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i,j,k+1)] + v_field[I1D(i,j,k+1)] ); vel_m = 0.5*( v_field[I1D(i,j,k-1)] + v_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 1.0/delta_z )*( mu*( vel_p*( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j-1,k+1)] )/delta_y )
                                                     - vel_m*( ( w_field[I1D(i,j+1,k-1)] - w_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 2.0/delta_z )*( 2.0*mu*( vel_p*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                         - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k+1)] ); vel_m = 0.5*( w_field[I1D(i,j,k-1)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( vel_p*( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i-1,j,k+1)] )/delta_x )
                                                                 - vel_m*( ( u_field[I1D(i+1,j,k-1)] - u_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k+1)] ); vel_m = 0.5*( w_field[I1D(i,j,k-1)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu*( vel_p*( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j-1,k+1)] )/delta_y )
                                                                 - vel_m*( ( v_field[I1D(i,j+1,k-1)] - v_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z -= ( 2.0/delta_z )*( ( 2.0/3.0 )*mu*( vel_p*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                 - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                /// Work of viscous stresses
                div_tauuvw = div_tauuvw_x + div_tauuvw_y + div_tauuvw_z; 
                /// Work of sources
                f_rhouvw = f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)];
                /// Viscous fluxes
                rhou_vis_flux[I1D(i,j,k)] = div_tau_xx + div_tau_yx + div_tau_zx + f_rhou_field[I1D(i,j,k)];
                rhov_vis_flux[I1D(i,j,k)] = div_tau_xy + div_tau_yy + div_tau_zy + f_rhov_field[I1D(i,j,k)];
                rhow_vis_flux[I1D(i,j,k)] = div_tau_xz + div_tau_yz + div_tau_zz + f_rhow_field[I1D(i,j,k)];
                rhoE_vis_flux[I1D(i,j,k)] = ( -1.0 )*div_q + div_tauuvw + f_rhouvw + f_rhoE_field[I1D(i,j,k)];
            }
        }
    }

    /// Update halo values
    rhou_vis_flux.update();
    rhov_vis_flux.update();
    rhow_vis_flux.update();
    rhoE_vis_flux.update();

};

void FlowSolverRHEA::sumInviscidViscousFluxes() {

    /// First Runge-Kutta step
    if(rk_step == 1) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    rho_rk1_flux[I1D(i,j,k)]  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                    rhou_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)]; 
                    rhov_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)]; 
                    rhow_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)]; 
                    rhoE_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)]; 
                }
            }
        }

        /// Update halo values
        rho_rk1_flux.update();
        rhou_rk1_flux.update();
        rhov_rk1_flux.update();
        rhow_rk1_flux.update();
        rhoE_rk1_flux.update();

    }

    /// Second Runge-Kutta step
    if(rk_step == 2) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    rho_rk2_flux[I1D(i,j,k)]  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                    rhou_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)]; 
                    rhov_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)]; 
                    rhow_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)]; 
                    rhoE_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)];
                }
            }
        }

        /// Update halo values
        rho_rk2_flux.update();
        rhou_rk2_flux.update();
        rhov_rk2_flux.update();
        rhow_rk2_flux.update();
        rhoE_rk2_flux.update();

    }

    /// Third Runge-Kutta step
    if(rk_step == 3) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    rho_rk3_flux[I1D(i,j,k)]  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                    rhou_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)]; 
                    rhov_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)]; 
                    rhow_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)]; 
                    rhoE_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)];
                }
            }
        }

        /// Update halo values
        rho_rk3_flux.update();
        rhou_rk3_flux.update();
        rhov_rk3_flux.update();
        rhow_rk3_flux.update();
        rhoE_rk3_flux.update();

    }

};

void FlowSolverRHEA::timeAdvanceConservedVariables() {

    /// Explicit third-order strong-stability-preserving Runge-Kutta (SSP-RK3) method:
    /// S. Gottlieb, C.-W. Shu & E. Tadmor.
    /// Strong stability-preserving high-order time discretization methods.
    /// SIAM Review 43, 89-112, 2001.

    /// First Runge-Kutta step
    if(rk_step == 1) {
        /// Inner points: rho, rhou, rhov, rhow and rhoE
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    rho_field[I1D(i,j,k)]  = rho_0_field[I1D(i,j,k)] + delta_t*rho_rk1_flux[I1D(i,j,k)];
                    rhou_field[I1D(i,j,k)] = rhou_0_field[I1D(i,j,k)] + delta_t*rhou_rk1_flux[I1D(i,j,k)];
                    rhov_field[I1D(i,j,k)] = rhov_0_field[I1D(i,j,k)] + delta_t*rhov_rk1_flux[I1D(i,j,k)];
                    rhow_field[I1D(i,j,k)] = rhow_0_field[I1D(i,j,k)] + delta_t*rhow_rk1_flux[I1D(i,j,k)];
                    rhoE_field[I1D(i,j,k)] = rhoE_0_field[I1D(i,j,k)] + delta_t*rhoE_rk1_flux[I1D(i,j,k)];
                }
            }
        }
    }

    /// Second Runge-Kutta step
    if(rk_step == 2) {
        /// Inner points: rho, rhou, rhov, rhow and rhoE
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    rho_field[I1D(i,j,k)]  = rho_0_field[I1D(i,j,k)] + ( delta_t/4.0 )*( rho_rk1_flux[I1D(i,j,k)] + rho_rk2_flux[I1D(i,j,k)] );
                    rhou_field[I1D(i,j,k)] = rhou_0_field[I1D(i,j,k)] + ( delta_t/4.0 )*( rhou_rk1_flux[I1D(i,j,k)] + rhou_rk2_flux[I1D(i,j,k)] );
                    rhov_field[I1D(i,j,k)] = rhov_0_field[I1D(i,j,k)] + ( delta_t/4.0 )*( rhov_rk1_flux[I1D(i,j,k)] + rhov_rk2_flux[I1D(i,j,k)] );
                    rhow_field[I1D(i,j,k)] = rhow_0_field[I1D(i,j,k)] + ( delta_t/4.0 )*( rhow_rk1_flux[I1D(i,j,k)] + rhow_rk2_flux[I1D(i,j,k)] );
                    rhoE_field[I1D(i,j,k)] = rhoE_0_field[I1D(i,j,k)] + ( delta_t/4.0 )*( rhoE_rk1_flux[I1D(i,j,k)] + rhoE_rk2_flux[I1D(i,j,k)] );
                }
            }
        }
    }

    /// Third Runge-Kutta step
    if(rk_step == 3) {
        /// Inner points: rho, rhou, rhov, rhow and rhoE
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    rho_field[I1D(i,j,k)]  = rho_0_field[I1D(i,j,k)] + ( delta_t/6.0 )*( rho_rk1_flux[I1D(i,j,k)] + rho_rk2_flux[I1D(i,j,k)] + 4.0*rho_rk3_flux[I1D(i,j,k)] );
                    rhou_field[I1D(i,j,k)] = rhou_0_field[I1D(i,j,k)] + ( delta_t/6.0 )*( rhou_rk1_flux[I1D(i,j,k)] + rhou_rk2_flux[I1D(i,j,k)] + 4.0*rhou_rk3_flux[I1D(i,j,k)] );
                    rhov_field[I1D(i,j,k)] = rhov_0_field[I1D(i,j,k)] + ( delta_t/6.0 )*( rhov_rk1_flux[I1D(i,j,k)] + rhov_rk2_flux[I1D(i,j,k)] + 4.0*rhov_rk3_flux[I1D(i,j,k)] );
                    rhow_field[I1D(i,j,k)] = rhow_0_field[I1D(i,j,k)] + ( delta_t/6.0 )*( rhow_rk1_flux[I1D(i,j,k)] + rhow_rk2_flux[I1D(i,j,k)] + 4.0*rhow_rk3_flux[I1D(i,j,k)] );
                    rhoE_field[I1D(i,j,k)] = rhoE_0_field[I1D(i,j,k)] + ( delta_t/6.0 )*( rhoE_rk1_flux[I1D(i,j,k)] + rhoE_rk2_flux[I1D(i,j,k)] + 4.0*rhoE_rk3_flux[I1D(i,j,k)] );
                }
            }
        }
    }

    /// Update halo values
    rho_field.update();
    rhou_field.update();
    rhov_field.update();
    rhow_field.update();
    rhoE_field.update();

};

void FlowSolverRHEA::outputCurrentStateData() {

    // This method needs to be implemented using HDF5!
   
    /// Obtain MPI information
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
    /// Open output data file
    ofstream output_file(output_data_file);

    /// Header string
    if(my_rank == 0) output_file << "#x [m],y [m],z [m],rho [kg/m3],u [m/s],v [m/s],w [m/s],E [J/kg],P [Pa],T [K],sos [m/s]" << endl;

    /// All (inner & boundary) points: x, y, z, rho, u, v, w, E, P, T, sos
    for(int p = 0; p < world_size; p++) {
        /// Serialize writing to file
        if(my_rank == p) {
            for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
                for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                    for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                        output_file << mesh->x[i] << ","; 
                        output_file << mesh->y[j] << ","; 
                        output_file << mesh->z[k] << ","; 
                        output_file << rho_field[I1D(i,j,k)] << ","; 
                        output_file << u_field[I1D(i,j,k)] << ","; 
                        output_file << v_field[I1D(i,j,k)] << ","; 
                        output_file << w_field[I1D(i,j,k)] << ","; 
                        output_file << E_field[I1D(i,j,k)] << ","; 
                        output_file << P_field[I1D(i,j,k)] << ","; 
                        output_file << T_field[I1D(i,j,k)] << ","; 
                        output_file << sos_field[I1D(i,j,k)] << endl; 
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /// Close output data file
    output_file.close();

};



////////// MAIN //////////
int main(int argc, char** argv) {

    /// Initialize MPI
    MPI_Init(&argc, &argv);
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /// Start RHEA simulation
    if(my_rank == 0) cout << "RHEA: START SIMULATION" << endl;

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
    for(int time_iter = flow_solver_RHEA.getCurrentTimeIteration(); time_iter <= flow_solver_RHEA.getFinalTimeIteration(); time_iter++) {

        /// Calculate time step
        flow_solver_RHEA.calculateTimeStep();
        if( ( flow_solver_RHEA.getCurrentTime() + flow_solver_RHEA.getTimeStep() ) > flow_solver_RHEA.getFinalTime() ) {
            flow_solver_RHEA.setTimeStep( flow_solver_RHEA.getFinalTime() - flow_solver_RHEA.getCurrentTime() );
        }

        /// Print time iteration information
        if(my_rank == 0) {
            cout << "Time iteration " << flow_solver_RHEA.getCurrentTimeIteration() << ":" 
                 << "time = " << flow_solver_RHEA.getCurrentTime() << " [s], "
                 << "time-step = " << flow_solver_RHEA.getTimeStep() << " [s]" << endl;
        }

        /// Output current state data to file (if criterion satisfied)
        if(flow_solver_RHEA.getCurrentTimeIteration()%flow_solver_RHEA.getOutputIteration() == 0) flow_solver_RHEA.outputCurrentStateData();

        /// Runge-Kutta time-integration steps
        flow_solver_RHEA.setCurrentRungeKuttaStep( 1 );
        for(int rk_step = flow_solver_RHEA.getCurrentRungeKuttaStep(); rk_step <= flow_solver_RHEA.getRungeKuttaOrder(); rk_step++) {

            /// Calculate source terms
            flow_solver_RHEA.calculateSourceTerms();

            /// Calculate inviscid & viscous (right-hand side) fluxes
            flow_solver_RHEA.calculateInviscidFluxes();
            flow_solver_RHEA.calculateViscousFluxes();
            flow_solver_RHEA.sumInviscidViscousFluxes();

            /// Advance conserved variables in time
            flow_solver_RHEA.timeAdvanceConservedVariables();

            /// Update boundary values
            flow_solver_RHEA.updateBoundaries();

            /// Calculate primitive variables from conserved variables
            flow_solver_RHEA.conservedToPrimitiveVariables();

            /// Calculate thermodynamics from primitive variables
            flow_solver_RHEA.calculateThermodynamicsFromPrimitiveVariables();

            /// Update Runge-Kutta step
            flow_solver_RHEA.setCurrentRungeKuttaStep( flow_solver_RHEA.getCurrentRungeKuttaStep() + 1 );

        }

        /// Update previous state of conserved variables
        flow_solver_RHEA.updatePreviousStateConservedVariables();

        /// Update time and time iteration
        flow_solver_RHEA.setCurrentTime( flow_solver_RHEA.getCurrentTime() + flow_solver_RHEA.getTimeStep() );
        flow_solver_RHEA.setCurrentTimeIteration( flow_solver_RHEA.getCurrentTimeIteration() + 1 );

        /// Check if simulation is completed: current_time > final_time
        if(flow_solver_RHEA.getCurrentTime() >= flow_solver_RHEA.getFinalTime() ) break;

    }

    /// Print time iteration information
    if(my_rank == 0) {
        cout << "Time iteration " << flow_solver_RHEA.getCurrentTimeIteration() << ":" 
             << "time = " << flow_solver_RHEA.getCurrentTime() << " [s]" << endl;
    }

    /// Output current state data to file
    flow_solver_RHEA.outputCurrentStateData();

    /// Destruct flow solver RHEA ... destructor is called automatically

    /// End RHEA simulation
    if(my_rank == 0) cout << "RHEA: END SIMULATION" << endl;

    /// Finalize MPI
    MPI_Finalize();

}

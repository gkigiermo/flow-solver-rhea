#include "FlowSolverRHEA.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
const double epsilon     = 1.0e-15;			// Small epsilon number (fixed)
const double Pi          = 2.0*asin(1.0);		// Pi number (fixed)
const int cout_presicion = 5;		                // Output precision (fixed)

////////// FlowSolverRHEA CLASS //////////
FlowSolverRHEA::FlowSolverRHEA() {};

FlowSolverRHEA::FlowSolverRHEA(const string name_configuration_file) : configuration_file(name_configuration_file) {

    /// Read configuration (input) file written in YAML language
    this->readConfigurationFile();
	
    /// Set value of fixed variables
    rk_order = 3;

    /// Construct (initialize) computational domain
    mesh = new domain(L_x, L_y, L_z, x_0, y_0, z_0, A_x, A_y, A_z, RHEA_NX, RHEA_NY, RHEA_NZ);

    /// Add boundary conditions to computational domain
    mesh->updateBocos(bocos_type);

    /// Construct (initialize) communication scheme
    topo = new comm_scheme(mesh, np_x, np_y, np_z);
    //if(topo->getRank() == 0) mesh->printDomain();
    //for(int p = 0; p < np_x*np_y*np_z; p++) topo->printCommSchemeToFile(p);

    // The lines below are temporary ... will need to be removed!
    _DEBUG_ = 0;
    _lNx_ = topo->getlNx();
    _lNy_ = topo->getlNy();
    _lNz_ = topo->getlNz();

    /// Set parallel topology of mesh coordinates
    x_field.setTopology(topo,"x [m]");
    y_field.setTopology(topo,"y [m]");
    z_field.setTopology(topo,"z [m]");
	
    /// Set parallel topology of primitive, conserved, thermodynamic and thermophysical variables	
    rho_field.setTopology(topo,"rho [kg/m^3]");
    u_field.setTopology(topo,"u [m/s]");
    v_field.setTopology(topo,"v [m/s]");
    w_field.setTopology(topo,"w [m/s]");
    E_field.setTopology(topo,"E [J/kg]");
    rhou_field.setTopology(topo,"rhou [kg/(m^2⋅s)]");
    rhov_field.setTopology(topo,"rhov [kg/(m^2⋅s)]");
    rhow_field.setTopology(topo,"rhow [kg/(m^2⋅s)]");
    rhoE_field.setTopology(topo,"rhoE [J/m^3]");
    P_field.setTopology(topo,"P [Pa]");
    T_field.setTopology(topo,"T [K]");
    sos_field.setTopology(topo,"sos [m/s]");
    mu_field.setTopology(topo,"mu [Pa·s]");
    kappa_field.setTopology(topo,"kappa [W/(m⋅K)]");

    /// Set parallel topology of time-integration variables	
    rho_0_field.setTopology(topo,"rho_0 [kg/m^3]");
    rhou_0_field.setTopology(topo,"rhou_0 [kg/(m^2⋅s)]");
    rhov_0_field.setTopology(topo,"rhov_0 [kg/(m^2⋅s)]");
    rhow_0_field.setTopology(topo,"rhow_0 [kg/(m^2⋅s)]");
    rhoE_0_field.setTopology(topo,"rhoE_0 [J/m^3]");    

    /// Set parallel topology of time-integration fluxes	
    rho_rk1_flux.setTopology(topo,"rho_rk1 [kg/(m^3⋅s)]");
    rho_rk2_flux.setTopology(topo,"rho_rk2 [kg/(m^3⋅s)]");
    rho_rk3_flux.setTopology(topo,"rho_rk3 [kg/(m^3⋅s)]");
    rhou_rk1_flux.setTopology(topo,"rhou_rk1 [kg/(m^2⋅s^2)]");
    rhou_rk2_flux.setTopology(topo,"rhou_rk2 [kg/(m^2⋅s^2)]");
    rhou_rk3_flux.setTopology(topo,"rhou_rk3 [kg/(m^2⋅s^2)]");    
    rhov_rk1_flux.setTopology(topo,"rhov_rk1 [kg/(m^2⋅s^2)]");
    rhov_rk2_flux.setTopology(topo,"rhov_rk2 [kg/(m^2⋅s^2)]");
    rhov_rk3_flux.setTopology(topo,"rhov_rk3 [kg/(m^2⋅s^2)]");    
    rhow_rk1_flux.setTopology(topo,"rhow_rk1 [kg/(m^2⋅s^2)]");
    rhow_rk2_flux.setTopology(topo,"rhow_rk2 [kg/(m^2⋅s^2)]");
    rhow_rk3_flux.setTopology(topo,"rhow_rk3 [kg/(m^2⋅s^2)]");
    rhoE_rk1_flux.setTopology(topo,"rhoE_rk1 [J/(m^3⋅s)]");
    rhoE_rk2_flux.setTopology(topo,"rhoE_rk2 [J/(m^3⋅s)]");
    rhoE_rk3_flux.setTopology(topo,"rhoE_rk3 [J/(m^3⋅s)]");

    /// Set parallel topology of inviscid fluxes	
    rho_inv_flux.setTopology(topo,"rho_inv [kg/(m^3⋅s)]");
    rhou_inv_flux.setTopology(topo,"rhou_inv [kg/(m^2⋅s^2)]");
    rhov_inv_flux.setTopology(topo,"rhov_inv [kg/(m^2⋅s^2)]");
    rhow_inv_flux.setTopology(topo,"rhow_inv [kg/(m^2⋅s^2)]");
    rhoE_inv_flux.setTopology(topo,"rhoE_inv [J/(m^3⋅s)]");

    /// Set parallel topology of viscous fluxes	
    rhou_vis_flux.setTopology(topo,"rhou_vis [kg/(m^2⋅s^2)]");
    rhov_vis_flux.setTopology(topo,"rhov_vis [kg/(m^2⋅s^2)]");
    rhow_vis_flux.setTopology(topo,"rhow_vis [kg/(m^2⋅s^2)]");
    rhoE_vis_flux.setTopology(topo,"rhoE_vis [J/(m^3⋅s)]");

    /// Set parallel topology of source terms
    f_rhou_field.setTopology(topo,"f_rhou [kg/(m^2⋅s^2)]");
    f_rhov_field.setTopology(topo,"f_rhov [kg/(m^2⋅s^2)]");
    f_rhow_field.setTopology(topo,"f_rhow [kg/(m^2⋅s^2)]");
    f_rhoE_field.setTopology(topo,"f_rhoE [J/(m^3⋅s)]");

    /// Fill x, y and z fields
    this->fillMeshCoordinateFields();

    /// Construct (initialize) HDF5 writer/reader
    char char_array[output_data_file.length()+1]; 
    strcpy(char_array,output_data_file.c_str());
    hdf5_data = new printer(topo,char_array);
    hdf5_data->addField(&x_field);
    hdf5_data->addField(&y_field);
    hdf5_data->addField(&z_field);
    hdf5_data->addField(&rho_field);
    hdf5_data->addField(&u_field);
    hdf5_data->addField(&v_field);
    hdf5_data->addField(&w_field);
    hdf5_data->addField(&E_field);
    hdf5_data->addField(&P_field);
    hdf5_data->addField(&T_field);
    hdf5_data->addField(&sos_field);
    hdf5_data->addField(&mu_field);
    hdf5_data->addField(&kappa_field);

};

FlowSolverRHEA::~FlowSolverRHEA() {

    /// Free mesh, topo, hdf5_data and bocos
    if(mesh != NULL) free(mesh);	
    if(topo != NULL) free(topo);
    if(hdf5_data != NULL) free(hdf5_data);
    free(bocos_type);	
    free(bocos_u);	
    free(bocos_v);	
    free(bocos_w);	
    free(bocos_P);	
    free(bocos_T);	

};

void FlowSolverRHEA::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid properties
    const YAML::Node & fluid_properties = configuration["fluid_properties"];
    R_specific = fluid_properties["R_specific"].as<double>();
    gamma      = fluid_properties["gamma"].as<double>();
    mu         = fluid_properties["mu"].as<double>();
    kappa      = fluid_properties["kappa"].as<double>();

    /// Problem parameters
    const YAML::Node & problem_parameters = configuration["problem_parameters"];
    x_0          = problem_parameters["x_0"].as<double>();
    y_0          = problem_parameters["y_0"].as<double>();
    z_0          = problem_parameters["z_0"].as<double>();
    L_x          = problem_parameters["L_x"].as<double>();
    L_y          = problem_parameters["L_y"].as<double>();
    L_z          = problem_parameters["L_z"].as<double>();
    current_time = problem_parameters["current_time"].as<double>();
    final_time   = problem_parameters["final_time"].as<double>();

    /// Computational parameters
    const YAML::Node & computational_parameters = configuration["computational_parameters"];
    num_grid_x        = computational_parameters["num_grid_x"].as<int>();
    num_grid_y        = computational_parameters["num_grid_y"].as<int>();
    num_grid_z        = computational_parameters["num_grid_z"].as<int>();
    A_x               = computational_parameters["A_x"].as<double>();
    A_y               = computational_parameters["A_y"].as<double>();
    A_z               = computational_parameters["A_z"].as<double>();
    CFL               = computational_parameters["CFL"].as<double>();
    current_time_iter = computational_parameters["current_time_iter"].as<int>();
    final_time_iter   = computational_parameters["final_time_iter"].as<int>();

    /// Boundary conditions
    string dummy_boco;
    const YAML::Node & boundary_conditons = configuration["boundary_conditions"];
    /// West
    dummy_boco = boundary_conditons["west_bc_type"].as<string>();
    if( dummy_boco == "DIRICHLET" ) {
        bocos_type[_WEST_] = _DIRICHLET_;
    } else if( dummy_boco == "NEUMANN" ) {
        bocos_type[_WEST_] = _NEUMANN_;
    } else if( dummy_boco == "PERIODIC" ) {
        bocos_type[_WEST_] = _PERIODIC_;
    } else {
        cout << "West boundary condition not available!" << endl;
    }
    /// East
    dummy_boco = boundary_conditons["east_bc_type"].as<string>();
    if( dummy_boco == "DIRICHLET" ) {
        bocos_type[_EAST_] = _DIRICHLET_;
    } else if( dummy_boco == "NEUMANN" ) {
        bocos_type[_EAST_] = _NEUMANN_;
    } else if( dummy_boco == "PERIODIC" ) {
        bocos_type[_EAST_] = _PERIODIC_;
    } else {
        cout << "East boundary condition not available!" << endl;
    }
    /// South
    dummy_boco = boundary_conditons["south_bc_type"].as<string>();
    if( dummy_boco == "DIRICHLET" ) {
        bocos_type[_SOUTH_] = _DIRICHLET_;
    } else if( dummy_boco == "NEUMANN" ) {
        bocos_type[_SOUTH_] = _NEUMANN_;
    } else if( dummy_boco == "PERIODIC" ) {
        bocos_type[_SOUTH_] = _PERIODIC_;
    } else {
        cout << "South boundary condition not available!" << endl;
    }
    /// North
    dummy_boco = boundary_conditons["north_bc_type"].as<string>();
    if( dummy_boco == "DIRICHLET" ) {
        bocos_type[_NORTH_] = _DIRICHLET_;
    } else if( dummy_boco == "NEUMANN" ) {
        bocos_type[_NORTH_] = _NEUMANN_;
    } else if( dummy_boco == "PERIODIC" ) {
        bocos_type[_NORTH_] = _PERIODIC_;
    } else {
        cout << "North boundary condition not available!" << endl;
    }
    /// Back
    dummy_boco = boundary_conditons["back_bc_type"].as<string>();
    if( dummy_boco == "DIRICHLET" ) {
        bocos_type[_BACK_] = _DIRICHLET_;
    } else if( dummy_boco == "NEUMANN" ) {
        bocos_type[_BACK_] = _NEUMANN_;
    } else if( dummy_boco == "PERIODIC" ) {
        bocos_type[_BACK_] = _PERIODIC_;
    } else {
        cout << "Back boundary condition not available!" << endl;
    }
    /// Front
    dummy_boco = boundary_conditons["front_bc_type"].as<string>();
    if( dummy_boco == "DIRICHLET" ) {
        bocos_type[_FRONT_] = _DIRICHLET_;
    } else if( dummy_boco == "NEUMANN" ) {
        bocos_type[_FRONT_] = _NEUMANN_;
    } else if( dummy_boco == "PERIODIC" ) {
        bocos_type[_FRONT_] = _PERIODIC_;
    } else {
        cout << "Front boundary condition not available!" << endl;
    }
    bocos_u[_WEST_]  = boundary_conditons["west_bc_u"].as<double>();
    bocos_u[_EAST_]  = boundary_conditons["east_bc_u"].as<double>();
    bocos_u[_SOUTH_] = boundary_conditons["south_bc_u"].as<double>();
    bocos_u[_NORTH_] = boundary_conditons["north_bc_u"].as<double>();
    bocos_u[_BACK_]  = boundary_conditons["back_bc_u"].as<double>();
    bocos_u[_FRONT_] = boundary_conditons["front_bc_u"].as<double>();
    bocos_v[_WEST_]  = boundary_conditons["west_bc_v"].as<double>();
    bocos_v[_EAST_]  = boundary_conditons["east_bc_v"].as<double>();
    bocos_v[_SOUTH_] = boundary_conditons["south_bc_v"].as<double>();
    bocos_v[_NORTH_] = boundary_conditons["north_bc_v"].as<double>();
    bocos_v[_BACK_]  = boundary_conditons["back_bc_v"].as<double>();
    bocos_v[_FRONT_] = boundary_conditons["front_bc_v"].as<double>();
    bocos_w[_WEST_]  = boundary_conditons["west_bc_w"].as<double>();
    bocos_w[_EAST_]  = boundary_conditons["east_bc_w"].as<double>();
    bocos_w[_SOUTH_] = boundary_conditons["south_bc_w"].as<double>();
    bocos_w[_NORTH_] = boundary_conditons["north_bc_w"].as<double>();
    bocos_w[_BACK_]  = boundary_conditons["back_bc_w"].as<double>();
    bocos_w[_FRONT_] = boundary_conditons["front_bc_w"].as<double>();
    bocos_P[_WEST_]  = boundary_conditons["west_bc_P"].as<double>();
    bocos_P[_EAST_]  = boundary_conditons["east_bc_P"].as<double>();
    bocos_P[_SOUTH_] = boundary_conditons["south_bc_P"].as<double>();
    bocos_P[_NORTH_] = boundary_conditons["north_bc_P"].as<double>();
    bocos_P[_BACK_]  = boundary_conditons["back_bc_P"].as<double>();
    bocos_P[_FRONT_] = boundary_conditons["front_bc_P"].as<double>();
    bocos_T[_WEST_]  = boundary_conditons["west_bc_T"].as<double>();
    bocos_T[_EAST_]  = boundary_conditons["east_bc_T"].as<double>();
    bocos_T[_SOUTH_] = boundary_conditons["south_bc_T"].as<double>();
    bocos_T[_NORTH_] = boundary_conditons["north_bc_T"].as<double>();
    bocos_T[_BACK_]  = boundary_conditons["back_bc_T"].as<double>();
    bocos_T[_FRONT_] = boundary_conditons["front_bc_T"].as<double>();

    /// Write/read file parameters
    const YAML::Node & write_read_parameters = configuration["write_read_parameters"];
    output_data_file       = write_read_parameters["output_data_file"].as<string>();
    output_frequency_iter  = write_read_parameters["output_frequency_iter"].as<int>();
    use_restart            = write_read_parameters["use_restart"].as<bool>();
    restart_data_file_iter = write_read_parameters["restart_data_file_iter"].as<int>();

    /// Parallelization scheme
    const YAML::Node & parallelization_scheme = configuration["parallelization_scheme"];
    np_x = parallelization_scheme["np_x"].as<int>();
    np_y = parallelization_scheme["np_y"].as<int>();
    np_z = parallelization_scheme["np_z"].as<int>();

};

void FlowSolverRHEA::fillMeshCoordinateFields() {

    /// All (inner & boundary) points: x, y and z
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                x_field[I1D(i,j,k)] = mesh->x[i];
                y_field[I1D(i,j,k)] = mesh->y[j];
                z_field[I1D(i,j,k)] = mesh->z[k];
            }
        }
    }

    /// Update halo values
    //x_field.update();
    //y_field.update();
    //z_field.update();

};

void FlowSolverRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// All (inner & boundary) points: u, v, w, P and T
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                double random_number = (double)rand()/(RAND_MAX + 1.0);
                u_field[I1D(i,j,k)] = random_number*1.0*sin( mesh->x[i] );
                v_field[I1D(i,j,k)] = random_number*1.0*sin( mesh->y[j] );
                w_field[I1D(i,j,k)] = random_number*1.0*sin( mesh->z[k] );
                P_field[I1D(i,j,k)] = 101325.0;
                T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( 1.0*R_specific );
            }
        }
    }

    /// Update halo values
    //u_field.update();
    //v_field.update();
    //w_field.update();
    //P_field.update();
    //T_field.update();

};

void FlowSolverRHEA::initializeFromRestart() {

    /// Read from file to restart solver (HDF5 data)
    hdf5_data->read(restart_data_file_iter);

    /// Update halo values
    x_field.update();
    y_field.update();
    z_field.update();
    rho_field.update();
    u_field.update();
    v_field.update();
    w_field.update();
    E_field.update();
    P_field.update();
    T_field.update();
    sos_field.update();

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
    //rho_field.update();
    //E_field.update();
    //sos_field.update();

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
    //rhou_field.update();
    //rhov_field.update();
    //rhow_field.update();
    //rhoE_field.update();

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
    //u_field.update();
    //v_field.update();
    //w_field.update();
    //E_field.update();

};

void FlowSolverRHEA::calculateThermodynamicsFromPrimitiveVariables() {

    /// Ideal-gas model:
    /// c_v = R_specific/(gamma - 1) is specific heat capacity at constant volume
    /// ke = (u*u + v*v + w*w)/2 is specific kinetic energy
    /// e = E - ke is specific internal energy
    /// P = e*rho*(gamma - 1) is pressure
    /// T = e/c_v is temperature
    /// sos = sqrt(gamma*P/rho) is speed of sound

    /// All (inner & boundary) points: ke, e, P, T and sos
    double c_v, ke, e;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                c_v                   = R_specific/( gamma - 1.0 );
                ke                    = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) ); 
                e                     = E_field[I1D(i,j,k)] - ke;
                P_field[I1D(i,j,k)]   = e*rho_field[I1D(i,j,k)]*( gamma - 1.0 ); 
                T_field[I1D(i,j,k)]   = e/c_v; 
                sos_field[I1D(i,j,k)] = sqrt( gamma*P_field[I1D(i,j,k)]/rho_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    //P_field.update();
    //T_field.update();
    //sos_field.update();

};

void FlowSolverRHEA::calculatePointDensityInternalEnergyFromPressureTemperature(const double &P, const double &T, double &rho, double &e) {

    /// Ideal-gas model:
    /// rho = P/(R_specific*T) is density
    /// e = P/(rho*((c_p/c_v) - 1)) is specific internal energy

    double c_v, c_p;
    this->calculateSpecificHeatCapacities( c_v, c_p );

    rho = P/( R_specific*T );
    e   = P/( rho*( ( c_p/c_v ) - 1.0 ) );

};

void FlowSolverRHEA::calculateSpecificHeatCapacities(double &c_v, double &c_p) {

    /// Ideal-gas model:
    /// c_v = R_specific/(gamma - 1)
    /// c_p = c_v*gamma

    c_v = R_specific/( gamma - 1.0 );
    c_p = c_v*gamma;

};

void FlowSolverRHEA::updateBoundaries() {

    /// General form: w_g*phi_g + w_in*phi_in = phi_b
    /// phi_g is ghost cell value
    /// phi_in is inner cell value
    /// phi_b is boundary value/flux
    /// w_g is ghost cell weight
    /// w_in is inner cell weight

    /// Declare variables & weights
    double rho_b, rhou_b, rhov_b, rhow_b, e_b, ke_b, rhoE_b, w_g, w_in;

    /// West boundary points: rho, rhou, rhov, rhow and rhoE
    if( bocos_type[_WEST_] == _DIRICHLET_ ) {
        w_g  = ( 1.0 )*( 1.0/2.0 );
        w_in = ( 1.0 )*( 1.0/2.0 );
    }
    if( bocos_type[_WEST_] == _NEUMANN_ ) {
        w_g  = (  1.0 )/( mesh->getGlobx(1) - mesh->getGlobx(0) );
        w_in = ( -1.0 )/( mesh->getGlobx(1) - mesh->getGlobx(0) );
    }
    this->calculatePointDensityInternalEnergyFromPressureTemperature( bocos_P[_WEST_], bocos_T[_WEST_], rho_b, e_b );
    rhou_b = rho_b*bocos_u[_WEST_];
    rhov_b = rho_b*bocos_v[_WEST_];
    rhow_b = rho_b*bocos_w[_WEST_];
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
        w_g  = (  1.0 )/( mesh->getGlobx(mesh->getGNx()+1) - mesh->getGlobx(mesh->getGNx()) );
        w_in = ( -1.0 )/( mesh->getGlobx(mesh->getGNx()+1) - mesh->getGlobx(mesh->getGNx()) );
    }
    this->calculatePointDensityInternalEnergyFromPressureTemperature( bocos_P[_EAST_], bocos_T[_EAST_], rho_b, e_b );
    rhou_b = rho_b*bocos_u[_EAST_];
    rhov_b = rho_b*bocos_v[_EAST_];
    rhow_b = rho_b*bocos_w[_EAST_];
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
        w_g  = (  1.0 )/( mesh->getGloby(1) - mesh->getGloby(0) );
        w_in = ( -1.0 )/( mesh->getGloby(1) - mesh->getGloby(0) );
    }
    this->calculatePointDensityInternalEnergyFromPressureTemperature( bocos_P[_SOUTH_], bocos_T[_SOUTH_], rho_b, e_b );
    rhou_b = rho_b*bocos_u[_SOUTH_];
    rhov_b = rho_b*bocos_v[_SOUTH_];
    rhow_b = rho_b*bocos_w[_SOUTH_];
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
        w_g  = (  1.0 )/( mesh->getGloby(mesh->getGNy()+1) - mesh->getGloby(mesh->getGNy()) );
        w_in = ( -1.0 )/( mesh->getGloby(mesh->getGNy()+1) - mesh->getGloby(mesh->getGNy()) );
    }
    this->calculatePointDensityInternalEnergyFromPressureTemperature( bocos_P[_NORTH_], bocos_T[_NORTH_], rho_b, e_b );
    rhou_b = rho_b*bocos_u[_NORTH_];
    rhov_b = rho_b*bocos_v[_NORTH_];
    rhow_b = rho_b*bocos_w[_NORTH_];
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
        w_g  = (  1.0 )/( mesh->getGlobz(1) - mesh->getGlobz(0) );
        w_in = ( -1.0 )/( mesh->getGlobz(1) - mesh->getGlobz(0) );
    }
    this->calculatePointDensityInternalEnergyFromPressureTemperature( bocos_P[_BACK_], bocos_T[_BACK_], rho_b, e_b );
    rhou_b = rho_b*bocos_u[_BACK_];
    rhov_b = rho_b*bocos_v[_BACK_];
    rhow_b = rho_b*bocos_w[_BACK_];
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
        w_g  = (  1.0 )/( mesh->getGlobz(mesh->getGNz()+1) - mesh->getGlobz(mesh->getGNz()) );
        w_in = ( -1.0 )/( mesh->getGlobz(mesh->getGNz()+1) - mesh->getGlobz(mesh->getGNz()) );
    }
    this->calculatePointDensityInternalEnergyFromPressureTemperature( bocos_P[_FRONT_], bocos_T[_FRONT_], rho_b, e_b );
    rhou_b = rho_b*bocos_u[_FRONT_];
    rhov_b = rho_b*bocos_v[_FRONT_];
    rhow_b = rho_b*bocos_w[_FRONT_];
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
    //rho_0_field.update();
    //rhou_0_field.update();
    //rhov_0_field.update();
    //rhow_0_field.update();
    //rhoE_0_field.update();

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

    /// Initialize to largest double value
    double local_delta_t = numeric_limits<double>::max();

    /// Inner points: find minimum (local) delta_t
    double c_v, c_p, Pr;
    double delta_x, delta_y, delta_z;
    double S_x, S_y, S_z;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                this->calculateSpecificHeatCapacities( c_v, c_p );
                /// Prandtl (Pr) number
                Pr = c_p/c_v;
                if(kappa_field[I1D(i,j,k)] > epsilon) Pr = c_p*mu_field[I1D(i,j,k)]/kappa_field[I1D(i,j,k)];
                /// Geometric stuff
                delta_x = 0.5*( mesh->x[i+1] - mesh->x[i-1] ); 
                delta_y = 0.5*( mesh->y[j+1] - mesh->y[j-1] ); 
                delta_z = 0.5*( mesh->z[k+1] - mesh->z[k-1] );                
                /// x-direction inviscid & viscous terms
                S_x           = abs( u_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_x/S_x );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_x, 2.0 )/( mu_field[I1D(i,j,k)]*( c_p/c_v ) + epsilon ) );
                /// y-direction inviscid & viscous terms
                S_y           = abs( v_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_y/S_y );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_y, 2.0 )/( mu_field[I1D(i,j,k)]*( c_p/c_v ) + epsilon ) );
                /// z-direction inviscid & viscous terms
                S_z           = abs( w_field[I1D(i,j,k)] ) + sos_field[I1D(i,j,k)];
                local_delta_t = min( local_delta_t, CFL*delta_z/S_z );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_z, 2.0 )/( mu_field[I1D(i,j,k)]*( c_p/c_v ) + epsilon ) );
            }
        }
    }

    /// Find minimum (global) delta_t
    double global_delta_t;
    MPI_Allreduce(&local_delta_t, &global_delta_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    /// Set new time step
    delta_t = global_delta_t;

};

void FlowSolverRHEA::calculateThermophysicalProperties() {
    
    /// Constant thermophysical properties model introduced via configuration file

    /// All (inner & boundary) points: mu and kappa
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                mu_field[I1D(i,j,k)]    = mu;
                kappa_field[I1D(i,j,k)] = kappa;
            }
        }
    }

    /// Update halo values
    //mu_field.update();
    //kappa_field.update()

};

void FlowSolverRHEA::calculateSourceTerms() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                //f_rhou_field[I1D(i,j,k)] = 0.0;
                f_rhou_field[I1D(i,j,k)] = 1.0;
                f_rhov_field[I1D(i,j,k)] = 0.0;
                f_rhow_field[I1D(i,j,k)] = 0.0;
                //f_rhoE_field[I1D(i,j,k)] = 0.0;
                f_rhoE_field[I1D(i,j,k)] = ( -1.0 )*( f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)] );
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

    double c_v, c_p;
    this->calculateSpecificHeatCapacities( c_v, c_p );

    double rho_bar = 0.5*( rho_L + rho_R );
    double a_bar   = 0.5*( a_L + a_R );
    double P_pvrs  = 0.5*( P_L + P_R ) - 0.5*( u_R - u_L )*rho_bar*a_bar;
    double P_star  = max( 0.0, P_pvrs );
    double q_L     = 1.0;
    if(P_star > P_L) q_L = sqrt( 1.0 + ( ( ( c_p/c_v ) + 1.0 )/( 2.0*( c_p/c_v ) ) )*( ( P_star/P_L ) - 1.0 ) );
    double q_R     = 1.0;
    if(P_star > P_R) q_R = sqrt( 1.0 + ( ( ( c_p/c_v ) + 1.0 )/( 2.0*( c_p/c_v ) ) )*( ( P_star/P_R ) - 1.0 ) );
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
                div_tau_xx  = ( 2.0/delta_x )*( 2.0*mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                         - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xx -= ( 2.0/delta_x )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                                 - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xx -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i+1,j-1,k)] )/delta_y )
                                                                                 - ( ( v_field[I1D(i-1,j+1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                div_tau_xx -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i+1,j,k-1)] )/delta_z )
                                                                                 - ( ( w_field[I1D(i-1,j,k+1)] - w_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                div_tau_xy  = ( 2.0/delta_x )*( mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i+1,j,k)] - v_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) ) 
                                                                     - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xy += ( 1.0/delta_x )*( mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i+1,j-1,k)] )/delta_y )
                                                                     - ( ( u_field[I1D(i-1,j+1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                div_tau_xz  = ( 2.0/delta_x )*( mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i+1,j,k)] - w_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                     - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                div_tau_xz += ( 1.0/delta_x )*( mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i+1,j,k-1)] )/delta_z )
                                                                     - ( ( u_field[I1D(i-1,j,k+1)] - u_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                div_tau_yx  = ( 2.0/delta_y )*( mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                     - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yx += ( 1.0/delta_y )*( mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i-1,j+1,k)] )/delta_x )
                                                                     - ( ( v_field[I1D(i+1,j-1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                div_tau_yy  = ( 2.0/delta_y )*( 2.0*mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                         - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yy -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i-1,j+1,k)] )/delta_x )
                                                                                 - ( ( u_field[I1D(i+1,j-1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                div_tau_yy -= ( 2.0/delta_y )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                                 - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yy -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j+1,k-1)] )/delta_z )
                                                                                 - ( ( w_field[I1D(i,j-1,k+1)] - w_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                div_tau_yz  = ( 2.0/delta_y )*( mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                     - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                div_tau_yz += ( 1.0/delta_y )*( mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j+1,k-1)] )/delta_z )
                                                                     - ( ( v_field[I1D(i,j-1,k+1)] - v_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                div_tau_zx  = ( 2.0/delta_z )*( mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                     - ( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                div_tau_zx += ( 1.0/delta_z )*( mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i-1,j,k+1)] )/delta_x )
                                                                     - ( ( w_field[I1D(i+1,j,k-1)] - w_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                div_tau_zy  = ( 2.0/delta_z )*( mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                     - ( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                div_tau_zy += ( 1.0/delta_z )*( mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j-1,k+1)] )/delta_y )
                                                                     - ( ( w_field[I1D(i,j+1,k-1)] - w_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                div_tau_zz  = ( 2.0/delta_z )*( 2.0*mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                         - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                div_tau_zz -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i-1,j,k+1)] )/delta_x )
                                                                                 - ( ( u_field[I1D(i+1,j,k-1)] - u_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                div_tau_zz -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j-1,k+1)] )/delta_y )
                                                                                 - ( ( v_field[I1D(i,j+1,k-1)] - v_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                div_tau_zz -= ( 2.0/delta_z )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( ( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                                 - ( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                /// Fourier term
                div_q = ( -1.0 )*( ( 2.0/delta_x )*( kappa_field[I1D(i,j,k)]*( ( ( T_field[I1D(i+1,j,k)] - T_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                             - ( ( T_field[I1D(i,j,k)] - T_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) )
                                 + ( 2.0/delta_y )*( kappa_field[I1D(i,j,k)]*( ( ( T_field[I1D(i,j+1,k)] - T_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                             - ( ( T_field[I1D(i,j,k)] - T_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) )
                                 + ( 2.0/delta_z )*( kappa_field[I1D(i,j,k)]*( ( ( T_field[I1D(i,j,k+1)] - T_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                             - ( ( T_field[I1D(i,j,k)] - T_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) ) );
                /// Divergence of tau*velocity terms
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x  = ( 2.0/delta_x )*( 2.0*mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                           - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x -= ( 2.0/delta_x )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                                   - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i+1,j,k)] ); vel_m = 0.5*( u_field[I1D(i-1,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i+1,j-1,k)] )/delta_y )
                                                                                   - vel_m*( ( v_field[I1D(i-1,j+1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                vel_p = 0.5*( u_field[I1D(i+1,j,k)] + u_field[I1D(i+1,j,k)] ); vel_m = 0.5*( u_field[I1D(i-1,j,k)] + u_field[I1D(i-1,j,k)] );
                div_tauuvw_x -= ( 1.0/delta_x )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i+1,j,k-1)] )/delta_z )
                                                                                   - vel_m*( ( w_field[I1D(i-1,j,k+1)] - w_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( v_field[I1D(i+1,j,k)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 2.0/delta_x )*( mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i+1,j,k)] - v_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                       - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i+1,j,k)] + v_field[I1D(i+1,j,k)] ); vel_m = 0.5*( v_field[I1D(i-1,j,k)] + v_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 1.0/delta_x )*( mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i+1,j-1,k)] )/delta_y )
                                                                       - vel_m*( ( u_field[I1D(i-1,j+1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_y ) ) );
                vel_p = 0.5*( w_field[I1D(i+1,j,k)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 2.0/delta_x )*( mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i+1,j,k)] - w_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] ) )
                                                                       - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) ) ) );
                vel_p = 0.5*( w_field[I1D(i+1,j,k)] + w_field[I1D(i+1,j,k)] ); vel_m = 0.5*( w_field[I1D(i-1,j,k)] + w_field[I1D(i-1,j,k)] );
                div_tauuvw_x += ( 1.0/delta_x )*( mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i+1,j,k-1)] )/delta_z )
                                                                       - vel_m*( ( u_field[I1D(i-1,j,k+1)] - u_field[I1D(i-1,j,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( u_field[I1D(i,j+1,k)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i,j-1,k)] );
                div_tauuvw_y  = ( 2.0/delta_y )*( mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                       - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i,j+1,k)] + u_field[I1D(i,j+1,k)] ); vel_m = 0.5*( u_field[I1D(i,j-1,k)] + u_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 1.0/delta_y )*( mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i-1,j+1,k)] )/delta_x )
                                                                       - vel_m*( ( v_field[I1D(i+1,j-1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 2.0/delta_y )*( 2.0*mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                           - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j+1,k)] ); vel_m = 0.5*( v_field[I1D(i,j-1,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i-1,j+1,k)] )/delta_x )
                                                                                   - vel_m*( ( u_field[I1D(i+1,j-1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_x ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y -= ( 2.0/delta_y )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                                   - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i,j+1,k)] + v_field[I1D(i,j+1,k)] ); vel_m = 0.5*( v_field[I1D(i,j-1,k)] + v_field[I1D(i,j-1,k)] );
                div_tauuvw_y -= ( 1.0/delta_y )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j+1,k-1)] )/delta_z )
                                                                                   - vel_m*( ( w_field[I1D(i,j-1,k+1)] - w_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( w_field[I1D(i,j+1,k)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 2.0/delta_y )*( mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] ) )
                                                                       - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) ) ) );
                vel_p = 0.5*( w_field[I1D(i,j+1,k)] + w_field[I1D(i,j+1,k)] ); vel_m = 0.5*( w_field[I1D(i,j-1,k)] + w_field[I1D(i,j-1,k)] );
                div_tauuvw_y += ( 1.0/delta_y )*( mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j+1,k-1)] )/delta_z )
                                                                       - vel_m*( ( v_field[I1D(i,j-1,k+1)] - v_field[I1D(i,j-1,k-1)] )/delta_z ) ) );
                vel_p = 0.5*( u_field[I1D(i,j,k+1)] + u_field[I1D(i,j,k)] ); vel_m = 0.5*( u_field[I1D(i,j,k)] + u_field[I1D(i,j,k-1)] );
                div_tauuvw_z  = ( 2.0/delta_z )*( mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                       - vel_m*( ( u_field[I1D(i,j,k)] - u_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                vel_p = 0.5*( u_field[I1D(i,j,k+1)] + u_field[I1D(i,j,k+1)] ); vel_m = 0.5*( u_field[I1D(i,j,k-1)] + u_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 1.0/delta_z )*( mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i-1,j,k+1)] )/delta_x )
                                                                       - vel_m*( ( w_field[I1D(i+1,j,k-1)] - w_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                vel_p = 0.5*( v_field[I1D(i,j,k+1)] + v_field[I1D(i,j,k)] ); vel_m = 0.5*( v_field[I1D(i,j,k)] + v_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 2.0/delta_z )*( mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                       - vel_m*( ( v_field[I1D(i,j,k)] - v_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                vel_p = 0.5*( v_field[I1D(i,j,k+1)] + v_field[I1D(i,j,k+1)] ); vel_m = 0.5*( v_field[I1D(i,j,k-1)] + v_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 1.0/delta_z )*( mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j-1,k+1)] )/delta_y )
                                                                       - vel_m*( ( w_field[I1D(i,j+1,k-1)] - w_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z += ( 2.0/delta_z )*( 2.0*mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                           - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k+1)] ); vel_m = 0.5*( w_field[I1D(i,j,k-1)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i-1,j,k+1)] )/delta_x )
                                                                                   - vel_m*( ( u_field[I1D(i+1,j,k-1)] - u_field[I1D(i-1,j,k-1)] )/delta_x ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k+1)] ); vel_m = 0.5*( w_field[I1D(i,j,k-1)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z -= ( 1.0/delta_z )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j-1,k+1)] )/delta_y )
                                                                                   - vel_m*( ( v_field[I1D(i,j+1,k-1)] - v_field[I1D(i,j-1,k-1)] )/delta_y ) ) );
                vel_p = 0.5*( w_field[I1D(i,j,k+1)] + w_field[I1D(i,j,k)] ); vel_m = 0.5*( w_field[I1D(i,j,k)] + w_field[I1D(i,j,k-1)] );
                div_tauuvw_z -= ( 2.0/delta_z )*( ( 2.0/3.0 )*mu_field[I1D(i,j,k)]*( vel_p*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] ) )
                                                                                   - vel_m*( ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) ) );
                /// Work of viscous stresses
                div_tauuvw = div_tauuvw_x + div_tauuvw_y + div_tauuvw_z; 
                /// Viscous fluxes
                rhou_vis_flux[I1D(i,j,k)] = div_tau_xx + div_tau_yx + div_tau_zx;
                rhov_vis_flux[I1D(i,j,k)] = div_tau_xy + div_tau_yy + div_tau_zy;
                rhow_vis_flux[I1D(i,j,k)] = div_tau_xz + div_tau_yz + div_tau_zz;
                rhoE_vis_flux[I1D(i,j,k)] = ( -1.0 )*div_q + div_tauuvw;
            }
        }
    }

    /// Update halo values
    rhou_vis_flux.update();
    rhov_vis_flux.update();
    rhow_vis_flux.update();
    rhoE_vis_flux.update();

};

void FlowSolverRHEA::sumInviscidViscousFluxesSourceTerms() {

    /// First Runge-Kutta step
    if(rk_step == 1) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        double f_rhouvw = 0.0;
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    /// Work of sources
                    f_rhouvw = f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)];
                    /// Sum all fluxes
                    rho_rk1_flux[I1D(i,j,k)]  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                    rhou_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)] + f_rhou_field[I1D(i,j,k)]; 
                    rhov_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]; 
                    rhow_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]; 
                    rhoE_rk1_flux[I1D(i,j,k)] = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)] + f_rhoE_field[I1D(i,j,k)] + f_rhouvw; 
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
        double f_rhouvw = 0.0;
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    /// Work of sources
                    f_rhouvw = f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)];
                    /// Sum all fluxes
                    rho_rk2_flux[I1D(i,j,k)]  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                    rhou_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)] + f_rhou_field[I1D(i,j,k)]; 
                    rhov_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]; 
                    rhow_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]; 
                    rhoE_rk2_flux[I1D(i,j,k)] = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)] + f_rhoE_field[I1D(i,j,k)] + f_rhouvw; 
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
        double f_rhouvw = 0.0;
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    /// Work of sources
                    f_rhouvw = f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)];
                    /// Sum all fluxes
                    rho_rk3_flux[I1D(i,j,k)]  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                    rhou_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)] + f_rhou_field[I1D(i,j,k)]; 
                    rhov_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]; 
                    rhow_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]; 
                    rhoE_rk3_flux[I1D(i,j,k)] = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)] + f_rhoE_field[I1D(i,j,k)] + f_rhouvw; 
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

    /// Write to file current solver state (HDF5 data)
    hdf5_data->write(current_time_iter);

};



////////// MAIN //////////
int main(int argc, char** argv) {

    /// Initialize MPI
    MPI_Init(&argc, &argv);
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /// Set output (cout) precision
    cout.precision( cout_presicion );

    /// Start RHEA simulation
    if(my_rank == 0) cout << "RHEA: START SIMULATION" << endl;

    /// Construct flow solver RHEA
    FlowSolverRHEA flow_solver_RHEA("configuration_file.yaml");

    /// Initialize variables from restart file or by setting initial conditions
    if( flow_solver_RHEA.getUseRestart() ) {

        /// Initialize from restart file
        flow_solver_RHEA.initializeFromRestart();


    } else {

        /// Set initial conditions
        flow_solver_RHEA.setInitialConditions();
    
        /// Initialize thermodynamics
        flow_solver_RHEA.initializeThermodynamics();

        /// Calculate thermophysical properties
        flow_solver_RHEA.calculateThermophysicalProperties();

    }

    /// Calculate conserved variables from primitive variables
    flow_solver_RHEA.primitiveToConservedVariables();

    /// Update previous state of conserved variables
    flow_solver_RHEA.updatePreviousStateConservedVariables();

    /// Iterate flow solver RHEA in time
    for(int time_iter = flow_solver_RHEA.getCurrentTimeIteration(); time_iter < flow_solver_RHEA.getFinalTimeIteration(); time_iter++) {

        /// Calculate time step
        flow_solver_RHEA.calculateTimeStep();
        if( ( flow_solver_RHEA.getCurrentTime() + flow_solver_RHEA.getTimeStep() ) > flow_solver_RHEA.getFinalTime() ) {
            flow_solver_RHEA.setTimeStep( flow_solver_RHEA.getFinalTime() - flow_solver_RHEA.getCurrentTime() );
        }

        /// Print time iteration information
        if(my_rank == 0) {
            cout << "Time iteration " << flow_solver_RHEA.getCurrentTimeIteration() << ": " 
                 << "time = " << scientific << flow_solver_RHEA.getCurrentTime() << " [s], "
                 << "time-step = " << scientific << flow_solver_RHEA.getTimeStep() << " [s]" << endl;
        }

        /// Output current state data to file (if criterion satisfied)
        if(flow_solver_RHEA.getCurrentTimeIteration()%flow_solver_RHEA.getOutputIterationFrequency() == 0) flow_solver_RHEA.outputCurrentStateData();

        /// Runge-Kutta time-integration steps
        flow_solver_RHEA.setCurrentRungeKuttaStep( 1 );
        for(int rk_step = flow_solver_RHEA.getCurrentRungeKuttaStep(); rk_step <= flow_solver_RHEA.getRungeKuttaOrder(); rk_step++) {

            /// Calculate thermophysical properties
            flow_solver_RHEA.calculateThermophysicalProperties();

            /// Calculate source terms
            flow_solver_RHEA.calculateSourceTerms();

            /// Calculate inviscid and viscous fluxes
            flow_solver_RHEA.calculateInviscidFluxes();
            flow_solver_RHEA.calculateViscousFluxes();

            /// Sum inviscid & viscous fluxes, and source terms (right-hand side)
            flow_solver_RHEA.sumInviscidViscousFluxesSourceTerms();

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

    /// Print time advancement information
    if(my_rank == 0) {
        cout << "Time advancement completed -> " 
             << "iteration = " << flow_solver_RHEA.getCurrentTimeIteration() << ", "
             << "time = " << scientific << flow_solver_RHEA.getCurrentTime() << " [s]" << endl;
        }

    /// Output current state data to file
    flow_solver_RHEA.outputCurrentStateData();

    /// Destruct flow solver RHEA ... destructor is called automatically

    /// End RHEA simulation
    if(my_rank == 0) cout << "RHEA: END SIMULATION" << endl;

    /// Finalize MPI
    MPI_Finalize();

}

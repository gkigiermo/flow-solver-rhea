#include "FlowSolverRHEA.hpp"

using namespace std;


////////// COMPILATION DIRECTIVES //////////
#define _PRESSURE_BASED_WAVE_SPEED_ESTIMATES_ 0		/// Select approach for estimating wave speeds

////////// FIXED PARAMETERS //////////
const double epsilon     = 1.0e-15;			/// Small epsilon number (fixed)
const double pi          = 2.0*asin(1.0);		/// pi number (fixed)
const int cout_presicion = 5;		                /// Output precision (fixed)


////////// FlowSolverRHEA CLASS //////////

FlowSolverRHEA::FlowSolverRHEA() {};

FlowSolverRHEA::FlowSolverRHEA(const string name_configuration_file) : configuration_file(name_configuration_file) {

    /// Read configuration (input) file
    this->readConfigurationFile();
	
    /// Set value of selected variables
    current_time      = 0.0;	/// Current time (restart will overwrite it)
    current_time_iter = 0;	/// Current time iteration (restart will overwrite it)
    averaging_time    = 0.0;	/// Current averaging time (restart will overwrite it)

    /// Construct (initialize) thermodynamic model
    if( thermodynamic_model == "IDEAL_GAS" ) {
        thermodynamics = new IdealGasModel( configuration_file );
    } else if( thermodynamic_model == "STIFFENED_GAS" ) {
        thermodynamics = new StiffenedGasModel( configuration_file );
    } else if( thermodynamic_model == "PENG_ROBINSON" ) {
        thermodynamics = new PengRobinsonModel( configuration_file );
    } else {
        cout << "Thermodynamic model not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /// Construct (initialize) transport coefficients model
    if( transport_coefficients_model == "CONSTANT" ) {
        transport_coefficients = new ConstantTransportCoefficients( configuration_file );
    } else if( transport_coefficients_model == "HIGH_PRESSURE" ) {
        transport_coefficients = new HighPressureTransportCoefficients( configuration_file );
    } else {
        cout << "Transport coefficients model not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /// Construct (initialize) Riemann solver
    if( riemann_solver_scheme == "CENTRAL" ) {
        riemann_solver = new CentralFluxApproximateRiemannSolver();
    } else if( riemann_solver_scheme == "MURMAN-ROE" ) {
        riemann_solver = new MurmanRoeFluxApproximateRiemannSolver();
    } else if( riemann_solver_scheme == "KGP" ) {
        riemann_solver = new KgpFluxApproximateRiemannSolver();
    } else if( riemann_solver_scheme == "HLL" ) {
        riemann_solver = new HllApproximateRiemannSolver();
    } else if( riemann_solver_scheme == "HLLC" ) {
        riemann_solver = new HllcApproximateRiemannSolver();
    } else if( riemann_solver_scheme == "HLLC+" ) {
        riemann_solver = new HllcPlusApproximateRiemannSolver();
    } else {
        cout << "Riemann solver not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /// Construct (initialize) computational domain
    mesh = new ComputationalDomain(L_x, L_y, L_z, x_0, y_0, z_0, A_x, A_y, A_z, num_grid_x, num_grid_y, num_grid_z);

    /// Set boundary conditions to computational domain
    mesh->setBocos(bocos_type);

    /// Construct (initialize) parallel topology
    topo = new ParallelTopology(mesh, np_x, np_y, np_z);
    //if(topo->getRank() == 0) mesh->printDomain();
    //for(int p = 0; p < np_x*np_y*np_z; p++) topo->printCommSchemeToFile(p);

    /// Set local mesh values for I1D macro
    _lNx_ = topo->getlNx();
    _lNy_ = topo->getlNy();
    _lNz_ = topo->getlNz();

    /// Set parallel topology of mesh coordinates
    x_field.setTopology(topo,"x");
    y_field.setTopology(topo,"y");
    z_field.setTopology(topo,"z");
	
    /// Set parallel topology of primitive, conserved, thermodynamic and thermophysical variables	
    rho_field.setTopology(topo,"rho");
    u_field.setTopology(topo,"u");
    v_field.setTopology(topo,"v");
    w_field.setTopology(topo,"w");
    E_field.setTopology(topo,"E");
    rhou_field.setTopology(topo,"rhou");
    rhov_field.setTopology(topo,"rhov");
    rhow_field.setTopology(topo,"rhow");
    rhoE_field.setTopology(topo,"rhoE");
    P_field.setTopology(topo,"P");
    T_field.setTopology(topo,"T");
    sos_field.setTopology(topo,"sos");
    mu_field.setTopology(topo,"mu");
    kappa_field.setTopology(topo,"kappa");

    /// Set parallel topology of time-integration variables	
    rho_0_field.setTopology(topo,"rho_0");
    rhou_0_field.setTopology(topo,"rhou_0");
    rhov_0_field.setTopology(topo,"rhov_0");
    rhow_0_field.setTopology(topo,"rhow_0");
    rhoE_0_field.setTopology(topo,"rhoE_0");    

    /// Set parallel topology of time-integration fluxes	
    rho_rk1_flux.setTopology(topo,"rho_rk1");
    rho_rk2_flux.setTopology(topo,"rho_rk2");
    rho_rk3_flux.setTopology(topo,"rho_rk3");
    rhou_rk1_flux.setTopology(topo,"rhou_rk1");
    rhou_rk2_flux.setTopology(topo,"rhou_rk2");
    rhou_rk3_flux.setTopology(topo,"rhou_rk3");    
    rhov_rk1_flux.setTopology(topo,"rhov_rk1");
    rhov_rk2_flux.setTopology(topo,"rhov_rk2");
    rhov_rk3_flux.setTopology(topo,"rhov_rk3");    
    rhow_rk1_flux.setTopology(topo,"rhow_rk1");
    rhow_rk2_flux.setTopology(topo,"rhow_rk2");
    rhow_rk3_flux.setTopology(topo,"rhow_rk3");
    rhoE_rk1_flux.setTopology(topo,"rhoE_rk1");
    rhoE_rk2_flux.setTopology(topo,"rhoE_rk2");
    rhoE_rk3_flux.setTopology(topo,"rhoE_rk3");

    /// Set parallel topology of inviscid fluxes	
    rho_inv_flux.setTopology(topo,"rho_inv");
    rhou_inv_flux.setTopology(topo,"rhou_inv");
    rhov_inv_flux.setTopology(topo,"rhov_inv");
    rhow_inv_flux.setTopology(topo,"rhow_inv");
    rhoE_inv_flux.setTopology(topo,"rhoE_inv");

    /// Set parallel topology of viscous fluxes	
    rhou_vis_flux.setTopology(topo,"rhou_vis");
    rhov_vis_flux.setTopology(topo,"rhov_vis");
    rhow_vis_flux.setTopology(topo,"rhow_vis");
    rhoE_vis_flux.setTopology(topo,"rhoE_vis");

    /// Set parallel topology of source terms
    f_rhou_field.setTopology(topo,"f_rhou");
    f_rhov_field.setTopology(topo,"f_rhov");
    f_rhow_field.setTopology(topo,"f_rhow");
    f_rhoE_field.setTopology(topo,"f_rhoE");

    /// Set parallel topology of time-averaged quantities
    avg_rho_field.setTopology(topo,"avg_rho");
    avg_rhou_field.setTopology(topo,"avg_rhou");
    avg_rhov_field.setTopology(topo,"avg_rhov");
    avg_rhow_field.setTopology(topo,"avg_rhow");
    avg_rhoE_field.setTopology(topo,"avg_rhoE");
    avg_P_field.setTopology(topo,"avg_P");
    avg_T_field.setTopology(topo,"avg_T");
    avg_sos_field.setTopology(topo,"avg_sos");
    rmsf_rho_field.setTopology(topo,"rmsf_rho");
    rmsf_rhou_field.setTopology(topo,"rmsf_rhou");
    rmsf_rhov_field.setTopology(topo,"rmsf_rhov");
    rmsf_rhow_field.setTopology(topo,"rmsf_rhow");
    rmsf_rhoE_field.setTopology(topo,"rmsf_rhoE");
    rmsf_P_field.setTopology(topo,"rmsf_P");
    rmsf_T_field.setTopology(topo,"rmsf_T");
    rmsf_sos_field.setTopology(topo,"rmsf_sos");

    /// Fill x, y and z fields
    this->fillMeshCoordinateFields();

    /// Construct (initialize) writer/reader
    char char_array[ output_data_file_name.length() + 1 ]; 
    strcpy( char_array,output_data_file_name.c_str() );
    writer_reader = new ManagerHDF5( topo, char_array, generate_xdmf );
    writer_reader->addAttributeDouble( "Time");
    writer_reader->addAttributeInt( "Iteration" );
    writer_reader->addAttributeDouble( "AveragingTime" );
    writer_reader->addField(&x_field);
    writer_reader->addField(&y_field);
    writer_reader->addField(&z_field);
    writer_reader->addField(&rho_field);
    writer_reader->addField(&u_field);
    writer_reader->addField(&v_field);
    writer_reader->addField(&w_field);
    writer_reader->addField(&E_field);
    writer_reader->addField(&P_field);
    writer_reader->addField(&T_field);
    writer_reader->addField(&sos_field);
    writer_reader->addField(&mu_field);
    writer_reader->addField(&kappa_field);
    writer_reader->addField(&avg_rho_field);
    writer_reader->addField(&avg_rhou_field);
    writer_reader->addField(&avg_rhov_field);
    writer_reader->addField(&avg_rhow_field);
    writer_reader->addField(&avg_rhoE_field);
    writer_reader->addField(&avg_P_field);
    writer_reader->addField(&avg_T_field);
    writer_reader->addField(&avg_sos_field);
    writer_reader->addField(&rmsf_rho_field);
    writer_reader->addField(&rmsf_rhou_field);
    writer_reader->addField(&rmsf_rhov_field);
    writer_reader->addField(&rmsf_rhow_field);
    writer_reader->addField(&rmsf_rhoE_field);
    writer_reader->addField(&rmsf_P_field);
    writer_reader->addField(&rmsf_T_field);
    writer_reader->addField(&rmsf_sos_field);

    /// Construct (initialize) timers
    timers = new ParallelTimer();
    timers->createTimer( "time_iteration_loop" );
    timers->createTimer( "calculate_time_step" );
    timers->createTimer( "output_solver_state" );
    timers->createTimer( "rk_iteration_loop" );
    timers->createTimer( "calculate_thermophysical_properties" );
    timers->createTimer( "calculate_inviscid_fluxes" );
    timers->createTimer( "calculate_viscous_fluxes" );
    timers->createTimer( "calculate_source_terms" );
    timers->createTimer( "sum_inviscid_viscous_fluxes_source_terms" );
    timers->createTimer( "time_advance_conserved_variables" );
    timers->createTimer( "update_boundaries" );
    timers->createTimer( "conserved_to_primitive_variables" );
    timers->createTimer( "calculate_thermodynamics_from_primitive_variables" );
    timers->createTimer( "update_time_averaged_quantities" );
    timers->createTimer( "update_previous_state_conserved_variables" );

};

FlowSolverRHEA::~FlowSolverRHEA() {

    /// Free thermodynamics, transport_coefficients, mesh, topo, writer_reader and timers
    if( thermodynamics != NULL ) free( thermodynamics );
    if( transport_coefficients != NULL ) free( transport_coefficients );
    if( riemann_solver != NULL ) free( riemann_solver );
    if( mesh != NULL ) free( mesh );	
    if( topo != NULL ) free( topo );
    if( writer_reader != NULL ) free( writer_reader );
    if( timers != NULL ) free( timers );

};

void FlowSolverRHEA::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid properties
    const YAML::Node & fluid_flow_properties = configuration["fluid_flow_properties"];
    thermodynamic_model          = fluid_flow_properties["thermodynamic_model"].as<string>();
    transport_coefficients_model = fluid_flow_properties["transport_coefficients_model"].as<string>();

    /// Problem parameters
    const YAML::Node & problem_parameters = configuration["problem_parameters"];
    x_0          = problem_parameters["x_0"].as<double>();
    y_0          = problem_parameters["y_0"].as<double>();
    z_0          = problem_parameters["z_0"].as<double>();
    L_x          = problem_parameters["L_x"].as<double>();
    L_y          = problem_parameters["L_y"].as<double>();
    L_z          = problem_parameters["L_z"].as<double>();
    final_time   = problem_parameters["final_time"].as<double>();

    /// Computational parameters
    const YAML::Node & computational_parameters = configuration["computational_parameters"];
    num_grid_x            = computational_parameters["num_grid_x"].as<int>();
    num_grid_y            = computational_parameters["num_grid_y"].as<int>();
    num_grid_z            = computational_parameters["num_grid_z"].as<int>();
    A_x                   = computational_parameters["A_x"].as<double>();
    A_y                   = computational_parameters["A_y"].as<double>();
    A_z                   = computational_parameters["A_z"].as<double>();
    CFL                   = computational_parameters["CFL"].as<double>();
    riemann_solver_scheme = computational_parameters["riemann_solver_scheme"].as<string>();
    final_time_iter       = computational_parameters["final_time_iter"].as<int>();

    /// Boundary conditions
    string dummy_type_boco;
    const YAML::Node & boundary_conditions = configuration["boundary_conditions"];
    /// West
    dummy_type_boco = boundary_conditions["west_bc"][0].as<string>();
    if( dummy_type_boco == "DIRICHLET" ) {
        bocos_type[_WEST_] = _DIRICHLET_;
    } else if( dummy_type_boco == "NEUMANN" ) {
        bocos_type[_WEST_] = _NEUMANN_;
    } else if( dummy_type_boco == "PERIODIC" ) {
        bocos_type[_WEST_] = _PERIODIC_;
    } else if( dummy_type_boco == "SUBSONIC_INFLOW" ) {
        bocos_type[_WEST_] = _SUBSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUBSONIC_OUTFLOW" ) {
        bocos_type[_WEST_] = _SUBSONIC_OUTFLOW_;
    } else {
        cout << "West boundary condition not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    bocos_u[_WEST_] = boundary_conditions["west_bc"][1].as<double>();
    bocos_v[_WEST_] = boundary_conditions["west_bc"][2].as<double>();
    bocos_w[_WEST_] = boundary_conditions["west_bc"][3].as<double>();
    bocos_P[_WEST_] = boundary_conditions["west_bc"][4].as<double>();
    bocos_T[_WEST_] = boundary_conditions["west_bc"][5].as<double>();
    /// East
    dummy_type_boco = boundary_conditions["east_bc"][0].as<string>();
    if( dummy_type_boco == "DIRICHLET" ) {
        bocos_type[_EAST_] = _DIRICHLET_;
    } else if( dummy_type_boco == "NEUMANN" ) {
        bocos_type[_EAST_] = _NEUMANN_;
    } else if( dummy_type_boco == "PERIODIC" ) {
        bocos_type[_EAST_] = _PERIODIC_;
    } else if( dummy_type_boco == "SUBSONIC_INFLOW" ) {
        bocos_type[_EAST_] = _SUBSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUBSONIC_OUTFLOW" ) {
        bocos_type[_EAST_] = _SUBSONIC_OUTFLOW_;
    } else {
        cout << "East boundary condition not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    bocos_u[_EAST_] = boundary_conditions["east_bc"][1].as<double>();
    bocos_v[_EAST_] = boundary_conditions["east_bc"][2].as<double>();
    bocos_w[_EAST_] = boundary_conditions["east_bc"][3].as<double>();
    bocos_P[_EAST_] = boundary_conditions["east_bc"][4].as<double>();
    bocos_T[_EAST_] = boundary_conditions["east_bc"][5].as<double>();
    /// South
    dummy_type_boco = boundary_conditions["south_bc"][0].as<string>();
    if( dummy_type_boco == "DIRICHLET" ) {
        bocos_type[_SOUTH_] = _DIRICHLET_;
    } else if( dummy_type_boco == "NEUMANN" ) {
        bocos_type[_SOUTH_] = _NEUMANN_;
    } else if( dummy_type_boco == "PERIODIC" ) {
        bocos_type[_SOUTH_] = _PERIODIC_;
    } else if( dummy_type_boco == "SUBSONIC_INFLOW" ) {
        bocos_type[_SOUTH_] = _SUBSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUBSONIC_OUTFLOW" ) {
        bocos_type[_SOUTH_] = _SUBSONIC_OUTFLOW_;
    } else {
        cout << "South boundary condition not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    bocos_u[_SOUTH_] = boundary_conditions["south_bc"][1].as<double>();
    bocos_v[_SOUTH_] = boundary_conditions["south_bc"][2].as<double>();
    bocos_w[_SOUTH_] = boundary_conditions["south_bc"][3].as<double>();
    bocos_P[_SOUTH_] = boundary_conditions["south_bc"][4].as<double>();
    bocos_T[_SOUTH_] = boundary_conditions["south_bc"][5].as<double>();
    /// North
    dummy_type_boco = boundary_conditions["north_bc"][0].as<string>();
    if( dummy_type_boco == "DIRICHLET" ) {
        bocos_type[_NORTH_] = _DIRICHLET_;
    } else if( dummy_type_boco == "NEUMANN" ) {
        bocos_type[_NORTH_] = _NEUMANN_;
    } else if( dummy_type_boco == "PERIODIC" ) {
        bocos_type[_NORTH_] = _PERIODIC_;
    } else if( dummy_type_boco == "SUBSONIC_INFLOW" ) {
        bocos_type[_NORTH_] = _SUBSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUBSONIC_OUTFLOW" ) {
        bocos_type[_NORTH_] = _SUBSONIC_OUTFLOW_;
    } else {
        cout << "North boundary condition not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    bocos_u[_NORTH_] = boundary_conditions["north_bc"][1].as<double>();
    bocos_v[_NORTH_] = boundary_conditions["north_bc"][2].as<double>();
    bocos_w[_NORTH_] = boundary_conditions["north_bc"][3].as<double>();
    bocos_P[_NORTH_] = boundary_conditions["north_bc"][4].as<double>();
    bocos_T[_NORTH_] = boundary_conditions["north_bc"][5].as<double>();
    /// Back
    dummy_type_boco = boundary_conditions["back_bc"][0].as<string>();
    if( dummy_type_boco == "DIRICHLET" ) {
        bocos_type[_BACK_] = _DIRICHLET_;
    } else if( dummy_type_boco == "NEUMANN" ) {
        bocos_type[_BACK_] = _NEUMANN_;
    } else if( dummy_type_boco == "PERIODIC" ) {
        bocos_type[_BACK_] = _PERIODIC_;
    } else if( dummy_type_boco == "SUBSONIC_INFLOW" ) {
        bocos_type[_BACK_] = _SUBSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUBSONIC_OUTFLOW" ) {
        bocos_type[_BACK_] = _SUBSONIC_OUTFLOW_;
    } else {
        cout << "Back boundary condition not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    bocos_u[_BACK_] = boundary_conditions["back_bc"][1].as<double>();
    bocos_v[_BACK_] = boundary_conditions["back_bc"][2].as<double>();
    bocos_w[_BACK_] = boundary_conditions["back_bc"][3].as<double>();
    bocos_P[_BACK_] = boundary_conditions["back_bc"][4].as<double>();
    bocos_T[_BACK_] = boundary_conditions["back_bc"][5].as<double>();
    /// Front
    dummy_type_boco = boundary_conditions["front_bc"][0].as<string>();
    if( dummy_type_boco == "DIRICHLET" ) {
        bocos_type[_FRONT_] = _DIRICHLET_;
    } else if( dummy_type_boco == "NEUMANN" ) {
        bocos_type[_FRONT_] = _NEUMANN_;
    } else if( dummy_type_boco == "PERIODIC" ) {
        bocos_type[_FRONT_] = _PERIODIC_;
    } else if( dummy_type_boco == "SUBSONIC_INFLOW" ) {
        bocos_type[_FRONT_] = _SUBSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUBSONIC_OUTFLOW" ) {
        bocos_type[_FRONT_] = _SUBSONIC_OUTFLOW_;
    } else {
        cout << "Front boundary condition not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    bocos_u[_FRONT_] = boundary_conditions["front_bc"][1].as<double>();
    bocos_v[_FRONT_] = boundary_conditions["front_bc"][2].as<double>();
    bocos_w[_FRONT_] = boundary_conditions["front_bc"][3].as<double>();
    bocos_P[_FRONT_] = boundary_conditions["front_bc"][4].as<double>();
    bocos_T[_FRONT_] = boundary_conditions["front_bc"][5].as<double>();

    /// Print/Write/Read file parameters
    const YAML::Node & print_write_read_parameters = configuration["print_write_read_parameters"];
    print_frequency_iter  = print_write_read_parameters["print_frequency_iter"].as<int>();
    output_data_file_name = print_write_read_parameters["output_data_file_name"].as<string>();
    output_frequency_iter = print_write_read_parameters["output_frequency_iter"].as<int>();
    generate_xdmf         = print_write_read_parameters["generate_xdmf"].as<bool>();
    use_restart           = print_write_read_parameters["use_restart"].as<bool>();
    restart_data_file     = print_write_read_parameters["restart_data_file"].as<string>();
    time_averaging_active = print_write_read_parameters["time_averaging_active"].as<bool>();
    reset_time_averaging  = print_write_read_parameters["reset_time_averaging"].as<bool>();

    /// Timers information
    const YAML::Node & timers_information = configuration["timers_information"];
    print_timers = timers_information["print_timers"].as<bool>();
    timers_information_file = timers_information["timers_information_file"].as<string>();

    /// Parallelization scheme
    const YAML::Node & parallelization_scheme = configuration["parallelization_scheme"];
    np_x = parallelization_scheme["np_x"].as<int>();
    np_y = parallelization_scheme["np_y"].as<int>();
    np_z = parallelization_scheme["np_z"].as<int>();

};

void FlowSolverRHEA::fillMeshCoordinateFields() {

    /// All (inner, halo, boundary) points: x, y and z
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                x_field[I1D(i,j,k)] = mesh->x[i];
                y_field[I1D(i,j,k)] = mesh->y[j];
                z_field[I1D(i,j,k)] = mesh->z[k];
            }
        }
    }

};

void FlowSolverRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// All (inner, halo, boundary): u, v, w, P and T
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

void FlowSolverRHEA::initializeFromRestart() {

    /// Read from file to restart solver: data, time and time iteration
    char char_restart_data_file[ restart_data_file.length() + 1 ]; 
    strcpy( char_restart_data_file, restart_data_file.c_str() );
    writer_reader->read( char_restart_data_file );
    current_time      = writer_reader->getAttributeDouble( "Time" );
    current_time_iter = writer_reader->getAttributeInt( "Iteration" );
    averaging_time    = writer_reader->getAttributeDouble( "AveragingTime" );
    if( reset_time_averaging ) {

	/// Reset time averaging
        averaging_time = 0.0;

	/// Reset avg and rmsf fields
        avg_rho_field   = 0.0;
        avg_rhou_field  = 0.0;
        avg_rhov_field  = 0.0;
        avg_rhow_field  = 0.0;
        avg_rhoE_field  = 0.0;
        avg_P_field     = 0.0;
        avg_T_field     = 0.0;
        avg_sos_field   = 0.0;
        rmsf_rho_field  = 0.0;
        rmsf_rhou_field = 0.0;
        rmsf_rhov_field = 0.0;
        rmsf_rhow_field = 0.0;
        rmsf_rhoE_field = 0.0;
        rmsf_P_field    = 0.0;
        rmsf_T_field    = 0.0;
        rmsf_sos_field  = 0.0;

    }

    /// Update halo values
    rho_field.update();
    u_field.update();
    v_field.update();
    w_field.update();
    E_field.update();
    P_field.update();
    T_field.update();
    sos_field.update();
    mu_field.update();
    kappa_field.update();
    avg_rho_field.update();
    avg_rhou_field.update();
    avg_rhov_field.update();
    avg_rhow_field.update();
    avg_rhoE_field.update();
    avg_P_field.update();
    avg_T_field.update();
    avg_sos_field.update();
    rmsf_rho_field.update();
    rmsf_rhou_field.update();
    rmsf_rhov_field.update();
    rmsf_rhow_field.update();
    rmsf_rhoE_field.update();
    rmsf_P_field.update();
    rmsf_T_field.update();
    rmsf_sos_field.update();
   
    /// Fill x, y and z fields
    this->fillMeshCoordinateFields();

};

void FlowSolverRHEA::initializeThermodynamics() {

    /// All (inner, halo, boundary): rho, E and sos
    double rho, e, ke;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho, e, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)] );
                rho_field[I1D(i,j,k)] = rho;
                ke                    = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) );
                E_field[I1D(i,j,k)]   = e + ke;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    rho_field.update();
    E_field.update();
    sos_field.update();

};

void FlowSolverRHEA::primitiveToConservedVariables() {

    /// All (inner, halo, boundary): rhou, rhov, rhow and rhoE
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

    /// Inner: u, v, w and E
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
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

    /// Inner points: P, T and sos
    double ke, e, P, T;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                ke                    = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) ); 
                e                     = E_field[I1D(i,j,k)] - ke;
                P = P_field[I1D(i,j,k)];	/// Initial pressure guess
                T = T_field[I1D(i,j,k)]; 	/// Initial temperature guess
                thermodynamics->calculatePressureTemperatureFromDensityInternalEnergy( P, T, rho_field[I1D(i,j,k)], e );
                P_field[I1D(i,j,k)]   = P; 
                T_field[I1D(i,j,k)]   = T; 
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    //P_field.update();
    //T_field.update();
    //sos_field.update();

};

void FlowSolverRHEA::updateBoundaries() {

    /// General form: w_g*phi_g + w_in*phi_in = phi_b
    /// phi_g is ghost cell value
    /// phi_in is inner cell value
    /// phi_b is boundary value/flux
    /// w_g is ghost cell weight
    /// w_in is inner cell weight

    /// Declare weights and ghost & inner values
    double wg_g = 0.0, wg_in = 0.0;
    double u_g, v_g, w_g, P_g, T_g, rho_g, e_g, ke_g, E_g;
    double u_in, v_in, w_in, P_in, T_in;

    /// West boundary points: rho, rhou, rhov, rhow and rhoE
    if( ( bocos_type[_WEST_] == _DIRICHLET_ ) or ( bocos_type[_WEST_] == _SUBSONIC_INFLOW_ ) ) {
        wg_g  = 1.0 - ( x_0 - mesh->getGlobx(0) )/( mesh->getGlobx(1) - mesh->getGlobx(0) );
        wg_in = 1.0 - ( mesh->getGlobx(1) - x_0 )/( mesh->getGlobx(1) - mesh->getGlobx(0) );
    }
    if( bocos_type[_WEST_] == _NEUMANN_ ) {
        wg_g  = (  1.0 )/( mesh->getGlobx(1) - mesh->getGlobx(0) );
        wg_in = ( -1.0 )/( mesh->getGlobx(1) - mesh->getGlobx(0) );
    }
    for(int i = topo->iter_bound[_WEST_][_INIX_]; i <= topo->iter_bound[_WEST_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_WEST_][_INIY_]; j <= topo->iter_bound[_WEST_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_WEST_][_INIZ_]; k <= topo->iter_bound[_WEST_][_ENDZ_]; k++) {
		/// Get/calculate inner values
                u_in   = u_field[I1D(i+1,j,k)];
                v_in   = v_field[I1D(i+1,j,k)];
                w_in   = w_field[I1D(i+1,j,k)];
                P_in   = P_field[I1D(i+1,j,k)];
                T_in   = T_field[I1D(i+1,j,k)];	
		/// Calculate ghost primitive variables
                u_g  = ( bocos_u[_WEST_] - wg_in*u_in )/wg_g;
                v_g  = ( bocos_v[_WEST_] - wg_in*v_in )/wg_g;
                w_g  = ( bocos_w[_WEST_] - wg_in*w_in )/wg_g;
                if( ( bocos_type[_WEST_] == _DIRICHLET_ ) and ( bocos_P[_WEST_] < 0.0 ) ) {
                    P_g = P_in;
                } else {
                    P_g = ( bocos_P[_WEST_] - wg_in*P_in )/wg_g;
                }
                if( ( bocos_type[_WEST_] == _DIRICHLET_ ) and ( bocos_T[_WEST_] < 0.0 ) ) {
                    T_g = T_in;
                } else {
                    T_g = ( bocos_T[_WEST_] - wg_in*T_in )/wg_g;
                }
                if( bocos_type[_WEST_] == _SUBSONIC_INFLOW_ ) {
                    double rho_in   = rho_field[I1D(i+1,j,k)]; 
                    double sos_in   = sos_field[I1D(i+1,j,k)];
                    double P_in_in  = P_field[I1D(i+2,j,k)]; 
                    double u_in_in  = u_field[I1D(i+2,j,k)]; 
                    double Delta    = mesh->x[i+2] - mesh->x[i+1];
                    double lambda_1 = u_in - sos_in;
		    double L_1      = lambda_1*( ( ( P_in_in - P_in )/Delta ) - rho_in*sos_in*( ( u_in_in - u_in )/Delta ) );
		    double L_5      = L_1;	// Steady-state velocity assumption
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_order )*0.5*( L_1 + L_5 );
                    T_g = ( bocos_T[_WEST_] - wg_in*T_in )/wg_g;
		} else if( bocos_type[_WEST_] == _SUBSONIC_OUTFLOW_ ) {
                    P_g = bocos_P[_WEST_]; 
                    double rho_in    = rho_field[I1D(i+1,j,k)]; 
                    double sos_in    = sos_field[I1D(i+1,j,k)];
                    double Ma_in     = u_in/sos_in;
                    double P_in_in   = P_field[I1D(i+2,j,k)]; 
                    //double u_in_in   = u_field[I1D(i+2,j,k)];
                    double v_in_in   = v_field[I1D(i+2,j,k)];
                    double w_in_in   = w_field[I1D(i+2,j,k)];
                    double rho_in_in = rho_field[I1D(i+2,j,k)];
                    double Delta_b   = mesh->x[i+1] - mesh->x[i];
                    double Delta     = mesh->x[i+2] - mesh->x[i+1];
                    double lambda_2  = u_in;
                    double lambda_3  = u_in;
                    double lambda_4  = u_in;
                    //double lambda_5  = u_in + sos_in;
		    double L_1       = sos_in*( 1.0 - Ma_in )*( P_in - P_g )/Delta_b;
		    double L_2       = lambda_2*( sos_in*sos_in*( ( rho_in_in - rho_in )/Delta ) - ( ( P_in_in - P_in )/Delta ) );
		    double L_3       = lambda_3*( ( v_in_in - v_in )/Delta );
		    double L_4       = lambda_4*( ( w_in_in - w_in )/Delta );
		    //double L_5       = lambda_5*( ( ( P_in_in - P_in )/Delta ) + rho_in*sos_in*( ( u_in_in - u_in )/Delta ) );
		    double L_5       = ( -1.0 )*L_1;		// Steady-state pressure assumption
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_4;
                    T_g   = thermodynamics->calculateTemperatureFromPressureDensity( P_g, rho_g );
                }
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho_g, e_g, P_g, T_g );
                ke_g = 0.5*( u_g*u_g + v_g*v_g + w_g*w_g );
                E_g  = e_g + ke_g;
		/// Update ghost conserved variables
                rho_field[I1D(i,j,k)]  = rho_g;
                rhou_field[I1D(i,j,k)] = rho_g*u_g;
                rhov_field[I1D(i,j,k)] = rho_g*v_g;
                rhow_field[I1D(i,j,k)] = rho_g*w_g;
                rhoE_field[I1D(i,j,k)] = rho_g*E_g;
		/// Update u, v, w, E, P, T and sos variables
                u_field[I1D(i,j,k)]   = u_g;
                v_field[I1D(i,j,k)]   = v_g;
                w_field[I1D(i,j,k)]   = w_g;
                E_field[I1D(i,j,k)]   = E_g;
                P_field[I1D(i,j,k)]   = P_g;
                T_field[I1D(i,j,k)]   = T_g;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
            }
        }
    }

    /// East boundary points: rho, rhou, rhov, rhow and rhoE
    if( ( bocos_type[_EAST_] == _DIRICHLET_ ) or ( bocos_type[_EAST_] == _SUBSONIC_INFLOW_ ) ) {
        wg_g  = 1.0 - ( mesh->getGlobx(mesh->getGNx()+1) - ( x_0 + L_x ) )/( mesh->getGlobx(mesh->getGNx()+1) - mesh->getGlobx(mesh->getGNx()) );
        wg_in = 1.0 - ( ( x_0 + L_x ) - mesh->getGlobx(mesh->getGNx()) )/( mesh->getGlobx(mesh->getGNx()+1) - mesh->getGlobx(mesh->getGNx()) );
    }
    if( bocos_type[_EAST_] == _NEUMANN_ ) {
        wg_g  = (  1.0 )/( mesh->getGlobx(mesh->getGNx()+1) - mesh->getGlobx(mesh->getGNx()) );
        wg_in = ( -1.0 )/( mesh->getGlobx(mesh->getGNx()+1) - mesh->getGlobx(mesh->getGNx()) );
    }
    for(int i = topo->iter_bound[_EAST_][_INIX_]; i <= topo->iter_bound[_EAST_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_EAST_][_INIY_]; j <= topo->iter_bound[_EAST_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_EAST_][_INIZ_]; k <= topo->iter_bound[_EAST_][_ENDZ_]; k++) {
		/// Get/calculate inner values
                u_in   = u_field[I1D(i-1,j,k)];
                v_in   = v_field[I1D(i-1,j,k)];
                w_in   = w_field[I1D(i-1,j,k)];
                P_in   = P_field[I1D(i-1,j,k)];
                T_in   = T_field[I1D(i-1,j,k)];	
		/// Calculate ghost primitive variables
                u_g  = ( bocos_u[_EAST_] - wg_in*u_in )/wg_g;
                v_g  = ( bocos_v[_EAST_] - wg_in*v_in )/wg_g;
                w_g  = ( bocos_w[_EAST_] - wg_in*w_in )/wg_g;
                if( ( bocos_type[_EAST_] == _DIRICHLET_ ) and ( bocos_P[_EAST_] < 0.0 ) ) {
                    P_g = P_in;
                } else {
                    P_g = ( bocos_P[_EAST_] - wg_in*P_in )/wg_g;
                }
                if( ( bocos_type[_EAST_] == _DIRICHLET_ ) and ( bocos_T[_EAST_] < 0.0 ) ) {
                    T_g = T_in;
                } else {
                    T_g = ( bocos_T[_EAST_] - wg_in*T_in )/wg_g;
                }
                if( bocos_type[_EAST_] == _SUBSONIC_INFLOW_ ) {
                    double rho_in   = rho_field[I1D(i-1,j,k)]; 
                    double sos_in   = sos_field[I1D(i-1,j,k)];
                    double P_in_in  = P_field[I1D(i-2,j,k)]; 
                    double u_in_in  = u_field[I1D(i-2,j,k)]; 
                    double Delta    = mesh->x[i-2] - mesh->x[i-1];
                    double lambda_1 = u_in - sos_in;
		    double L_1      = lambda_1*( ( ( P_in_in - P_in )/Delta ) - rho_in*sos_in*( ( u_in_in - u_in )/Delta ) );
		    double L_5      = L_1;	// Steady-state velocity assumption
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_order )*0.5*( L_1 + L_5 );
                    T_g = ( bocos_T[_EAST_] - wg_in*T_in )/wg_g;
		} else if( bocos_type[_EAST_] == _SUBSONIC_OUTFLOW_ ) {
                    P_g = bocos_P[_EAST_]; 
                    double rho_in    = rho_field[I1D(i-1,j,k)]; 
                    double sos_in    = sos_field[I1D(i-1,j,k)];
                    double Ma_in     = u_in/sos_in;
                    double P_in_in   = P_field[I1D(i-2,j,k)]; 
                    //double u_in_in   = u_field[I1D(i-2,j,k)];
                    double v_in_in   = v_field[I1D(i-2,j,k)];
                    double w_in_in   = w_field[I1D(i-2,j,k)];
                    double rho_in_in = rho_field[I1D(i-2,j,k)];
                    double Delta_b   = mesh->x[i] - mesh->x[i-1];
                    double Delta     = mesh->x[i-2] - mesh->x[i-1];
                    double lambda_2  = u_in;
                    double lambda_3  = u_in;
                    double lambda_4  = u_in;
                    //double lambda_5  = u_in + sos_in;
		    double L_1       = sos_in*( 1.0 - Ma_in )*( P_in - P_g )/Delta_b;
		    double L_2       = lambda_2*( sos_in*sos_in*( ( rho_in_in - rho_in )/Delta ) - ( ( P_in_in - P_in )/Delta ) );
		    double L_3       = lambda_3*( ( v_in_in - v_in )/Delta );
		    double L_4       = lambda_4*( ( w_in_in - w_in )/Delta );
		    //double L_5       = lambda_5*( ( ( P_in_in - P_in )/Delta ) + rho_in*sos_in*( ( u_in_in - u_in )/Delta ) );
		    double L_5       = ( -1.0 )*L_1;		// Steady-state pressure assumption
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_4;
                    T_g   = thermodynamics->calculateTemperatureFromPressureDensity( P_g, rho_g );
                }
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho_g, e_g, P_g, T_g );
                ke_g = 0.5*( u_g*u_g + v_g*v_g + w_g*w_g );
                E_g  = e_g + ke_g;
		/// Update ghost conserved variables
                rho_field[I1D(i,j,k)]  = rho_g;
                rhou_field[I1D(i,j,k)] = rho_g*u_g;
                rhov_field[I1D(i,j,k)] = rho_g*v_g;
                rhow_field[I1D(i,j,k)] = rho_g*w_g;
                rhoE_field[I1D(i,j,k)] = rho_g*E_g;
		/// Update u, v, w, E, P, T and sos variables
                u_field[I1D(i,j,k)]   = u_g;
                v_field[I1D(i,j,k)]   = v_g;
                w_field[I1D(i,j,k)]   = w_g;
                E_field[I1D(i,j,k)]   = E_g;
                P_field[I1D(i,j,k)]   = P_g;
                T_field[I1D(i,j,k)]   = T_g;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
            }
        }
    }

    /// South boundary points: rho, rhou, rhov, rhow and rhoE
    if( ( bocos_type[_SOUTH_] == _DIRICHLET_ ) or ( bocos_type[_SOUTH_] == _SUBSONIC_INFLOW_ ) ) {
        wg_g  = 1.0 - ( y_0 - mesh->getGloby(0) )/( mesh->getGloby(1) - mesh->getGloby(0) );
        wg_in = 1.0 - ( mesh->getGloby(1) - y_0 )/( mesh->getGloby(1) - mesh->getGloby(0) );
    }
    if( bocos_type[_SOUTH_] == _NEUMANN_ ) {
        wg_g  = (  1.0 )/( mesh->getGloby(1) - mesh->getGloby(0) );
        wg_in = ( -1.0 )/( mesh->getGloby(1) - mesh->getGloby(0) );
    }
    for(int i = topo->iter_bound[_SOUTH_][_INIX_]; i <= topo->iter_bound[_SOUTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_SOUTH_][_INIY_]; j <= topo->iter_bound[_SOUTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_SOUTH_][_INIZ_]; k <= topo->iter_bound[_SOUTH_][_ENDZ_]; k++) {
		/// Get/calculate inner values
                u_in   = u_field[I1D(i,j+1,k)];
                v_in   = v_field[I1D(i,j+1,k)];
                w_in   = w_field[I1D(i,j+1,k)];
                P_in   = P_field[I1D(i,j+1,k)];
                T_in   = T_field[I1D(i,j+1,k)];	
		/// Calculate ghost primitive variables
                u_g  = ( bocos_u[_SOUTH_] - wg_in*u_in )/wg_g;
                v_g  = ( bocos_v[_SOUTH_] - wg_in*v_in )/wg_g;
                w_g  = ( bocos_w[_SOUTH_] - wg_in*w_in )/wg_g;
                if( ( bocos_type[_SOUTH_] == _DIRICHLET_ ) and ( bocos_P[_SOUTH_] < 0.0 ) ) {
                    P_g = P_in;
                } else {
                    P_g = ( bocos_P[_SOUTH_] - wg_in*P_in )/wg_g;
                }
                if( ( bocos_type[_SOUTH_] == _DIRICHLET_ ) and ( bocos_T[_SOUTH_] < 0.0 ) ) {
                    T_g = T_in;
                } else {
                    T_g = ( bocos_T[_SOUTH_] - wg_in*T_in )/wg_g;
                }
                if( bocos_type[_SOUTH_] == _SUBSONIC_INFLOW_ ) {
                    double rho_in   = rho_field[I1D(i,j+1,k)]; 
                    double sos_in   = sos_field[I1D(i,j+1,k)];
                    double P_in_in  = P_field[I1D(i,j+2,k)]; 
                    double v_in_in  = v_field[I1D(i,j+2,k)]; 
                    double Delta    = mesh->y[j+2] - mesh->y[j+1];
                    double lambda_1 = v_in - sos_in;
		    double L_1      = lambda_1*( ( ( P_in_in - P_in )/Delta ) - rho_in*sos_in*( ( v_in_in - v_in )/Delta ) );
		    double L_5      = L_1;	// Steady-state velocity assumption
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_order )*0.5*( L_1 + L_5 );
                    T_g = ( bocos_T[_SOUTH_] - wg_in*T_in )/wg_g;
		} else if( bocos_type[_SOUTH_] == _SUBSONIC_OUTFLOW_ ) {
                    P_g = bocos_P[_SOUTH_]; 
                    double rho_in    = rho_field[I1D(i,j+1,k)]; 
                    double sos_in    = sos_field[I1D(i,j+1,k)];
                    double Ma_in     = u_in/sos_in;
                    double P_in_in   = P_field[I1D(i,j+2,k)]; 
                    double u_in_in   = u_field[I1D(i,j+2,k)];
                    //double v_in_in   = v_field[I1D(i,j+2,k)];
                    double w_in_in   = w_field[I1D(i,j+2,k)];
                    double rho_in_in = rho_field[I1D(i,j+2,k)];
                    double Delta_b   = mesh->y[j+1] - mesh->y[j];
                    double Delta     = mesh->y[j+2] - mesh->y[j+1];
                    double lambda_2  = v_in;
                    double lambda_3  = v_in;
                    double lambda_4  = v_in;
                    //double lambda_5  = v_in + sos_in;
		    double L_1       = sos_in*( 1.0 - Ma_in )*( P_in - P_g )/Delta_b;
		    double L_2       = lambda_2*( sos_in*sos_in*( ( rho_in_in - rho_in )/Delta ) - ( ( P_in_in - P_in )/Delta ) );
		    double L_3       = lambda_3*( ( u_in_in - u_in )/Delta );
		    double L_4       = lambda_4*( ( w_in_in - w_in )/Delta );
		    //double L_5       = lambda_5*( ( ( P_in_in - P_in )/Delta ) + rho_in*sos_in*( ( v_in_in - v_in )/Delta ) );
		    double L_5       = ( -1.0 )*L_1;		// Steady-state pressure assumption
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_4;
                    T_g   = thermodynamics->calculateTemperatureFromPressureDensity( P_g, rho_g );
                }
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho_g, e_g, P_g, T_g );
                ke_g = 0.5*( u_g*u_g + v_g*v_g + w_g*w_g );
                E_g  = e_g + ke_g;
		/// Update ghost conserved variables
                rho_field[I1D(i,j,k)]  = rho_g;
                rhou_field[I1D(i,j,k)] = rho_g*u_g;
                rhov_field[I1D(i,j,k)] = rho_g*v_g;
                rhow_field[I1D(i,j,k)] = rho_g*w_g;
                rhoE_field[I1D(i,j,k)] = rho_g*E_g;
		/// Update u, v, w, E, P, T and sos variables
                u_field[I1D(i,j,k)]   = u_g;
                v_field[I1D(i,j,k)]   = v_g;
                w_field[I1D(i,j,k)]   = w_g;
                E_field[I1D(i,j,k)]   = E_g;
                P_field[I1D(i,j,k)]   = P_g;
                T_field[I1D(i,j,k)]   = T_g;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
            }
        }
    }

    /// North boundary points: rho, rhou, rhov, rhow and rhoE
    if( ( bocos_type[_NORTH_] == _DIRICHLET_ ) or ( bocos_type[_NORTH_] == _SUBSONIC_INFLOW_ ) ) {
        wg_g  = 1.0 - ( mesh->getGloby(mesh->getGNy()+1) - ( y_0 + L_y ) )/( mesh->getGloby(mesh->getGNy()+1) - mesh->getGloby(mesh->getGNy()) );
        wg_in = 1.0 - ( ( y_0 + L_y ) - mesh->getGloby(mesh->getGNy()) )/( mesh->getGloby(mesh->getGNy()+1) - mesh->getGloby(mesh->getGNy()) );
    }
    if( bocos_type[_NORTH_] == _NEUMANN_ ) {
        wg_g  = (  1.0 )/( mesh->getGloby(mesh->getGNy()+1) - mesh->getGloby(mesh->getGNy()) );
        wg_in = ( -1.0 )/( mesh->getGloby(mesh->getGNy()+1) - mesh->getGloby(mesh->getGNy()) );
    }
    for(int i = topo->iter_bound[_NORTH_][_INIX_]; i <= topo->iter_bound[_NORTH_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_NORTH_][_INIY_]; j <= topo->iter_bound[_NORTH_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_NORTH_][_INIZ_]; k <= topo->iter_bound[_NORTH_][_ENDZ_]; k++) {
		/// Get/calculate inner values
                u_in   = u_field[I1D(i,j-1,k)];
                v_in   = v_field[I1D(i,j-1,k)];
                w_in   = w_field[I1D(i,j-1,k)];
                P_in   = P_field[I1D(i,j-1,k)];
                T_in   = T_field[I1D(i,j-1,k)];	
		/// Calculate ghost primitive variables
                u_g  = ( bocos_u[_NORTH_] - wg_in*u_in )/wg_g;
                v_g  = ( bocos_v[_NORTH_] - wg_in*v_in )/wg_g;
                w_g  = ( bocos_w[_NORTH_] - wg_in*w_in )/wg_g;
                if( ( bocos_type[_NORTH_] == _DIRICHLET_ ) and ( bocos_P[_NORTH_] < 0.0 ) ) {
                    P_g = P_in;
                } else {
                    P_g = ( bocos_P[_NORTH_] - wg_in*P_in )/wg_g;
                }
                if( ( bocos_type[_NORTH_] == _DIRICHLET_ ) and ( bocos_T[_NORTH_] < 0.0 ) ) {
                    T_g = T_in;
                } else {
                    T_g = ( bocos_T[_NORTH_] - wg_in*T_in )/wg_g;
                }
                if( bocos_type[_NORTH_] == _SUBSONIC_INFLOW_ ) {
                    double rho_in   = rho_field[I1D(i,j-1,k)]; 
                    double sos_in   = sos_field[I1D(i,j-1,k)];
                    double P_in_in  = P_field[I1D(i,j-2,k)]; 
                    double v_in_in  = v_field[I1D(i,j-2,k)]; 
                    double Delta    = mesh->y[j-2] - mesh->y[j-1];
                    double lambda_1 = v_in - sos_in;
		    double L_1      = lambda_1*( ( ( P_in_in - P_in )/Delta ) - rho_in*sos_in*( ( v_in_in - v_in )/Delta ) );
		    double L_5      = L_1;	// Steady-state velocity assumption
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_order )*0.5*( L_1 + L_5 );
                    T_g = ( bocos_T[_NORTH_] - wg_in*T_in )/wg_g;
		} else if( bocos_type[_NORTH_] == _SUBSONIC_OUTFLOW_ ) {
                    P_g = bocos_P[_NORTH_]; 
                    double rho_in    = rho_field[I1D(i,j-1,k)]; 
                    double sos_in    = sos_field[I1D(i,j-1,k)];
                    double Ma_in     = u_in/sos_in;
                    double P_in_in   = P_field[I1D(i,j-2,k)]; 
                    double u_in_in   = u_field[I1D(i,j-2,k)];
                    //double v_in_in   = v_field[I1D(i,j-2,k)];
                    double w_in_in   = w_field[I1D(i,j-2,k)];
                    double rho_in_in = rho_field[I1D(i,j-2,k)];
                    double Delta_b   = mesh->y[j] - mesh->y[j-1];
                    double Delta     = mesh->y[j-2] - mesh->y[j-1];
                    double lambda_2  = v_in;
                    double lambda_3  = v_in;
                    double lambda_4  = v_in;
                    //double lambda_5  = v_in + sos_in;
		    double L_1       = sos_in*( 1.0 - Ma_in )*( P_in - P_g )/Delta_b;
		    double L_2       = lambda_2*( sos_in*sos_in*( ( rho_in_in - rho_in )/Delta ) - ( ( P_in_in - P_in )/Delta ) );
		    double L_3       = lambda_3*( ( u_in_in - u_in )/Delta );
		    double L_4       = lambda_4*( ( w_in_in - w_in )/Delta );
		    //double L_5       = lambda_5*( ( ( P_in_in - P_in )/Delta ) + rho_in*sos_in*( ( v_in_in - v_in )/Delta ) );
		    double L_5       = ( -1.0 )*L_1;		// Steady-state pressure assumption
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_4;
                    T_g   = thermodynamics->calculateTemperatureFromPressureDensity( P_g, rho_g );
                }
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho_g, e_g, P_g, T_g );
                ke_g = 0.5*( u_g*u_g + v_g*v_g + w_g*w_g );
                E_g  = e_g + ke_g;
		/// Update ghost conserved variables
                rho_field[I1D(i,j,k)]  = rho_g;
                rhou_field[I1D(i,j,k)] = rho_g*u_g;
                rhov_field[I1D(i,j,k)] = rho_g*v_g;
                rhow_field[I1D(i,j,k)] = rho_g*w_g;
                rhoE_field[I1D(i,j,k)] = rho_g*E_g;
		/// Update u, v, w, E, P, T and sos variables
                u_field[I1D(i,j,k)]   = u_g;
                v_field[I1D(i,j,k)]   = v_g;
                w_field[I1D(i,j,k)]   = w_g;
                E_field[I1D(i,j,k)]   = E_g;
                P_field[I1D(i,j,k)]   = P_g;
                T_field[I1D(i,j,k)]   = T_g;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
            }
        }
    }

    /// Back boundary points: rho, rhou, rhov, rhow and rhoE
    if( ( bocos_type[_BACK_] == _DIRICHLET_ ) or ( bocos_type[_BACK_] == _SUBSONIC_INFLOW_ ) ) {
        wg_g  = 1.0 - ( z_0 - mesh->getGlobz(0) )/( mesh->getGlobz(1) - mesh->getGlobz(0) );
        wg_in = 1.0 - ( mesh->getGlobz(1) - z_0 )/( mesh->getGlobz(1) - mesh->getGlobz(0) );
    }
    if( bocos_type[_BACK_] == _NEUMANN_ ) {
        wg_g  = (  1.0 )/( mesh->getGlobz(1) - mesh->getGlobz(0) );
        wg_in = ( -1.0 )/( mesh->getGlobz(1) - mesh->getGlobz(0) );
    }
    for(int i = topo->iter_bound[_BACK_][_INIX_]; i <= topo->iter_bound[_BACK_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_BACK_][_INIY_]; j <= topo->iter_bound[_BACK_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_BACK_][_INIZ_]; k <= topo->iter_bound[_BACK_][_ENDZ_]; k++) {
		/// Get/calculate inner values
                u_in   = u_field[I1D(i,j,k+1)];
                v_in   = v_field[I1D(i,j,k+1)];
                w_in   = w_field[I1D(i,j,k+1)];
                P_in   = P_field[I1D(i,j,k+1)];
                T_in   = T_field[I1D(i,j,k+1)];	
		/// Calculate ghost primitive variables
                u_g  = ( bocos_u[_BACK_] - wg_in*u_in )/wg_g;
                v_g  = ( bocos_v[_BACK_] - wg_in*v_in )/wg_g;
                w_g  = ( bocos_w[_BACK_] - wg_in*w_in )/wg_g;
                if( ( bocos_type[_BACK_] == _DIRICHLET_ ) and ( bocos_P[_BACK_] < 0.0 ) ) {
                    P_g = P_in;
                } else {
                    P_g = ( bocos_P[_BACK_] - wg_in*P_in )/wg_g;
                }
                if( ( bocos_type[_BACK_] == _DIRICHLET_ ) and ( bocos_T[_BACK_] < 0.0 ) ) {
                    T_g = T_in;
                } else {
                    T_g = ( bocos_T[_BACK_] - wg_in*T_in )/wg_g;
                }
                if( bocos_type[_BACK_] == _SUBSONIC_INFLOW_ ) {
                    double rho_in   = rho_field[I1D(i,j,k+1)]; 
                    double sos_in   = sos_field[I1D(i,j,k+1)];
                    double P_in_in  = P_field[I1D(i,j,k+2)]; 
                    double w_in_in  = w_field[I1D(i,j,k+2)]; 
                    double Delta    = mesh->z[k+2] - mesh->z[k+1];
                    double lambda_1 = w_in - sos_in;
		    double L_1      = lambda_1*( ( ( P_in_in - P_in )/Delta ) - rho_in*sos_in*( ( w_in_in - w_in )/Delta ) );
		    double L_5      = L_1;	// Steady-state velocity assumption
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_order )*0.5*( L_1 + L_5 );
                    T_g = ( bocos_T[_BACK_] - wg_in*T_in )/wg_g;
		} else if( bocos_type[_BACK_] == _SUBSONIC_OUTFLOW_ ) {
                    P_g = bocos_P[_BACK_]; 
                    double rho_in    = rho_field[I1D(i,j,k+1)]; 
                    double sos_in    = sos_field[I1D(i,j,k+1)];
                    double Ma_in     = u_in/sos_in;
                    double P_in_in   = P_field[I1D(i,j,k+2)]; 
                    double u_in_in   = u_field[I1D(i,j,k+2)];
                    double v_in_in   = v_field[I1D(i,j,k+2)];
                    //double w_in_in   = w_field[I1D(i,j,k+2)];
                    double rho_in_in = rho_field[I1D(i,j,k+2)];
                    double Delta_b   = mesh->z[k+1] - mesh->z[k];
                    double Delta     = mesh->z[k+2] - mesh->z[k+1];
                    double lambda_2  = w_in;
                    double lambda_3  = w_in;
                    double lambda_4  = w_in;
                    //double lambda_5  = w_in + sos_in;
		    double L_1       = sos_in*( 1.0 - Ma_in )*( P_in - P_g )/Delta_b;
		    double L_2       = lambda_2*( sos_in*sos_in*( ( rho_in_in - rho_in )/Delta ) - ( ( P_in_in - P_in )/Delta ) );
		    double L_3       = lambda_3*( ( v_in_in - v_in )/Delta );
		    double L_4       = lambda_4*( ( u_in_in - u_in )/Delta );
		    //double L_5       = lambda_5*( ( ( P_in_in - P_in )/Delta ) + rho_in*sos_in*( ( w_in_in - w_in )/Delta ) );
		    double L_5       = ( -1.0 )*L_1;		// Steady-state pressure assumption
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_4;
                    T_g   = thermodynamics->calculateTemperatureFromPressureDensity( P_g, rho_g );
                }
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho_g, e_g, P_g, T_g );
                ke_g = 0.5*( u_g*u_g + v_g*v_g + w_g*w_g );
                E_g  = e_g + ke_g;
		/// Update ghost conserved variables
                rho_field[I1D(i,j,k)]  = rho_g;
                rhou_field[I1D(i,j,k)] = rho_g*u_g;
                rhov_field[I1D(i,j,k)] = rho_g*v_g;
                rhow_field[I1D(i,j,k)] = rho_g*w_g;
                rhoE_field[I1D(i,j,k)] = rho_g*E_g;
		/// Update u, v, w, E, P, T and sos variables
                u_field[I1D(i,j,k)]   = u_g;
                v_field[I1D(i,j,k)]   = v_g;
                w_field[I1D(i,j,k)]   = w_g;
                E_field[I1D(i,j,k)]   = E_g;
                P_field[I1D(i,j,k)]   = P_g;
                T_field[I1D(i,j,k)]   = T_g;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
            }
        }
    }

    /// Front boundary points: rho, rhou, rhov, rhow and rhoE
    if( ( bocos_type[_FRONT_] == _DIRICHLET_ ) or ( bocos_type[_FRONT_] == _SUBSONIC_INFLOW_ ) ) {
        wg_g  = 1.0 - ( mesh->getGlobz(mesh->getGNz()+1) - ( z_0 + L_z ) )/( mesh->getGlobz(mesh->getGNz()+1) - mesh->getGlobz(mesh->getGNz()) );
        wg_in = 1.0 - ( ( z_0 + L_z ) - mesh->getGlobz(mesh->getGNz()) )/( mesh->getGlobz(mesh->getGNz()+1) - mesh->getGlobz(mesh->getGNz()) );
    }
    if( bocos_type[_FRONT_] == _NEUMANN_ ) {
        wg_g  = (  1.0 )/( mesh->getGlobz(mesh->getGNz()+1) - mesh->getGlobz(mesh->getGNz()) );
        wg_in = ( -1.0 )/( mesh->getGlobz(mesh->getGNz()+1) - mesh->getGlobz(mesh->getGNz()) );
    }
    for(int i = topo->iter_bound[_FRONT_][_INIX_]; i <= topo->iter_bound[_FRONT_][_ENDX_]; i++) {
        for(int j = topo->iter_bound[_FRONT_][_INIY_]; j <= topo->iter_bound[_FRONT_][_ENDY_]; j++) {
            for(int k = topo->iter_bound[_FRONT_][_INIZ_]; k <= topo->iter_bound[_FRONT_][_ENDZ_]; k++) {
		/// Get/calculate inner values
                u_in   = u_field[I1D(i,j,k-1)];
                v_in   = v_field[I1D(i,j,k-1)];
                w_in   = w_field[I1D(i,j,k-1)];
                P_in   = P_field[I1D(i,j,k-1)];
                T_in   = T_field[I1D(i,j,k-1)];	
		/// Calculate ghost primitive variables
                u_g  = ( bocos_u[_FRONT_] - wg_in*u_in )/wg_g;
                v_g  = ( bocos_v[_FRONT_] - wg_in*v_in )/wg_g;
                w_g  = ( bocos_w[_FRONT_] - wg_in*w_in )/wg_g;
                if( ( bocos_type[_FRONT_] == _DIRICHLET_ ) and ( bocos_P[_FRONT_] < 0.0 ) ) {
                    P_g = P_in;
                } else {
                    P_g = ( bocos_P[_FRONT_] - wg_in*P_in )/wg_g;
                }
                if( ( bocos_type[_FRONT_] == _DIRICHLET_ ) and ( bocos_T[_FRONT_] < 0.0 ) ) {
                    T_g = T_in;
                } else {
                    T_g = ( bocos_T[_FRONT_] - wg_in*T_in )/wg_g;
                }
                if( bocos_type[_FRONT_] == _SUBSONIC_INFLOW_ ) {
                    double rho_in   = rho_field[I1D(i,j,k-1)]; 
                    double sos_in   = sos_field[I1D(i,j,k-1)];
                    double P_in_in  = P_field[I1D(i,j,k-2)]; 
                    double w_in_in  = w_field[I1D(i,j,k-2)]; 
                    double Delta    = mesh->z[k-2] - mesh->z[k-1];
                    double lambda_1 = w_in - sos_in;
		    double L_1      = lambda_1*( ( ( P_in_in - P_in )/Delta ) - rho_in*sos_in*( ( w_in_in - w_in )/Delta ) );
		    double L_5      = L_1;	// Steady-state velocity assumption
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_order )*0.5*( L_1 + L_5 );
                    T_g = ( bocos_T[_FRONT_] - wg_in*T_in )/wg_g;
		} else if( bocos_type[_FRONT_] == _SUBSONIC_OUTFLOW_ ) {
                    P_g = bocos_P[_FRONT_]; 
                    double rho_in    = rho_field[I1D(i,j,k-1)]; 
                    double sos_in    = sos_field[I1D(i,j,k-1)];
                    double Ma_in     = u_in/sos_in;
                    double P_in_in   = P_field[I1D(i,j,k-2)]; 
                    double u_in_in   = u_field[I1D(i,j,k-2)];
                    double v_in_in   = v_field[I1D(i,j,k-2)];
                    //double w_in_in   = w_field[I1D(i,j,k-2)];
                    double rho_in_in = rho_field[I1D(i,j,k-2)];
                    double Delta_b   = mesh->z[k] - mesh->z[k-1];
                    double Delta     = mesh->z[k-2] - mesh->z[k-1];
                    double lambda_2  = w_in;
                    double lambda_3  = w_in;
                    double lambda_4  = w_in;
                    //double lambda_5  = w_in + sos_in;
		    double L_1       = sos_in*( 1.0 - Ma_in )*( P_in - P_g )/Delta_b;
		    double L_2       = lambda_2*( sos_in*sos_in*( ( rho_in_in - rho_in )/Delta ) - ( ( P_in_in - P_in )/Delta ) );
		    double L_3       = lambda_3*( ( v_in_in - v_in )/Delta );
		    double L_4       = lambda_4*( ( u_in_in - u_in )/Delta );
		    //double L_5       = lambda_5*( ( ( P_in_in - P_in )/Delta ) + rho_in*sos_in*( ( w_in_in - w_in )/Delta ) );
		    double L_5       = ( -1.0 )*L_1;		// Steady-state pressure assumption
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_order )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_order )*L_4;
                    T_g   = thermodynamics->calculateTemperatureFromPressureDensity( P_g, rho_g );
                }
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho_g, e_g, P_g, T_g );
                ke_g = 0.5*( u_g*u_g + v_g*v_g + w_g*w_g );
                E_g  = e_g + ke_g;
		/// Update ghost conserved variables
                rho_field[I1D(i,j,k)]  = rho_g;
                rhou_field[I1D(i,j,k)] = rho_g*u_g;
                rhov_field[I1D(i,j,k)] = rho_g*v_g;
                rhow_field[I1D(i,j,k)] = rho_g*w_g;
                rhoE_field[I1D(i,j,k)] = rho_g*E_g;
		/// Update u, v, w, E, P, T and sos variables
                u_field[I1D(i,j,k)]   = u_g;
                v_field[I1D(i,j,k)]   = v_g;
                w_field[I1D(i,j,k)]   = w_g;
                E_field[I1D(i,j,k)]   = E_g;
                P_field[I1D(i,j,k)]   = P_g;
                T_field[I1D(i,j,k)]   = T_g;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
            }
        }
    }

    /// Update halo values
    rho_field.update();
    rhou_field.update();
    rhov_field.update();
    rhow_field.update();
    rhoE_field.update();
    u_field.update();
    v_field.update();
    w_field.update();
    E_field.update();
    P_field.update();
    T_field.update();
    sos_field.update();

};

void FlowSolverRHEA::updatePreviousStateConservedVariables() {

    /// All (inner, halo, boundary) points: rho_0, rhou_0 rhov_0, rhow_0 and rhoE_0
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

    /// Initialize to largest double value
    double local_delta_t = numeric_limits<double>::max();

    /// Inner points: find minimum (local) delta_t
    double sos, c_v, c_p, Pr;
    double delta_x, delta_y, delta_z;
    double S_x, S_y, S_z;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Speed of sound
		sos = sos_field[I1D(i,j,k)];
                /// Heat capacities
                thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                /// Prandtl (Pr) number
                Pr = c_p/c_v;
                if(kappa_field[I1D(i,j,k)] > epsilon) Pr = c_p*mu_field[I1D(i,j,k)]/kappa_field[I1D(i,j,k)];
                /// Geometric stuff
                delta_x = 0.5*( mesh->x[i+1] - mesh->x[i-1] ); 
                delta_y = 0.5*( mesh->y[j+1] - mesh->y[j-1] ); 
                delta_z = 0.5*( mesh->z[k+1] - mesh->z[k-1] );                
                /// x-direction inviscid & viscous terms
                S_x           = abs( u_field[I1D(i,j,k)] ) + sos;
                local_delta_t = min( local_delta_t, CFL*delta_x/S_x );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_x, 2.0 )/( mu_field[I1D(i,j,k)]*( c_p/c_v ) + epsilon ) );
                /// y-direction inviscid & viscous terms
                S_y           = abs( v_field[I1D(i,j,k)] ) + sos;
                local_delta_t = min( local_delta_t, CFL*delta_y/S_y );
                local_delta_t = min( local_delta_t, CFL*Pr*rho_field[I1D(i,j,k)]*pow( delta_y, 2.0 )/( mu_field[I1D(i,j,k)]*( c_p/c_v ) + epsilon ) );
                /// z-direction inviscid & viscous terms
                S_z           = abs( w_field[I1D(i,j,k)] ) + sos;
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

void FlowSolverRHEA::calculateTransportCoefficients() {
    
    /// All (inner, halo, boundary) points: mu and kappa
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                mu_field[I1D(i,j,k)]    = transport_coefficients->calculateDynamicViscosity( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                kappa_field[I1D(i,j,k)] = transport_coefficients->calculateThermalConductivity( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    mu_field.update();
    kappa_field.update();

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
                //f_rhoE_field[I1D(i,j,k)] = ( -1.0 )*( f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    //f_rhou_field.update();
    //f_rhov_field.update();
    //f_rhow_field.update();
    //f_rhoE_field.update();

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
                rho_F_L  = rho_L*u_L;
                rho_F_R  = rho_R*u_R;
                rho_U_L  = rho_L;
                rho_U_R  = rho_R;
                rho_F_p  = riemann_solver->calculateIntercellFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 1;
                rhou_F_L = rho_L*u_L*u_L + P_L;
                rhou_F_R = rho_R*u_R*u_R + P_R;
                rhou_U_L = rho_L*u_L;
                rhou_U_R = rho_R*u_R;
                rhou_F_p = riemann_solver->calculateIntercellFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_L*u_L*v_L;
                rhov_F_R = rho_R*u_R*v_R;
                rhov_U_L = rho_L*v_L;
                rhov_U_R = rho_R*v_R;
                rhov_F_p = riemann_solver->calculateIntercellFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_L*u_L*w_L;
                rhow_F_R = rho_R*u_R*w_R;
                rhow_U_L = rho_L*w_L;
                rhow_U_R = rho_R*w_R;
                rhow_F_p = riemann_solver->calculateIntercellFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_L*u_L*E_L + u_L*P_L;
                rhoE_F_R = rho_R*u_R*E_R + u_R*P_R;
                rhoE_U_L = rho_L*E_L;
                rhoE_U_R = rho_R*E_R;
                rhoE_F_p = riemann_solver->calculateIntercellFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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
                rho_F_L  = rho_L*u_L;
                rho_F_R  = rho_R*u_R;
                rho_U_L  = rho_L;
                rho_U_R  = rho_R;
                rho_F_m  = riemann_solver->calculateIntercellFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 1;
                rhou_F_L = rho_L*u_L*u_L + P_L;
                rhou_F_R = rho_R*u_R*u_R + P_R;
                rhou_U_L = rho_L*u_L;
                rhou_U_R = rho_R*u_R;
                rhou_F_m = riemann_solver->calculateIntercellFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_L*u_L*v_L;
                rhov_F_R = rho_R*u_R*v_R;
                rhov_U_L = rho_L*v_L;
                rhov_U_R = rho_R*v_R;
                rhov_F_m = riemann_solver->calculateIntercellFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_L*u_L*w_L;
                rhow_F_R = rho_R*u_R*w_R;
                rhow_U_L = rho_L*w_L;
                rhow_U_R = rho_R*w_R;
                rhow_F_m = riemann_solver->calculateIntercellFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_L*u_L*E_L + u_L*P_L;
                rhoE_F_R = rho_R*u_R*E_R + u_R*P_R;
                rhoE_U_L = rho_L*E_L;
                rhoE_U_R = rho_R*E_R;
                rhoE_F_m = riemann_solver->calculateIntercellFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// Fluxes x-direction
                rho_inv_flux[I1D(i,j,k)]  = ( rho_F_p - rho_F_m )/delta_x;
                rhou_inv_flux[I1D(i,j,k)] = ( rhou_F_p - rhou_F_m )/delta_x;
                rhov_inv_flux[I1D(i,j,k)] = ( rhov_F_p - rhov_F_m )/delta_x;
                rhow_inv_flux[I1D(i,j,k)] = ( rhow_F_p - rhow_F_m )/delta_x;
                rhoE_inv_flux[I1D(i,j,k)] = ( rhoE_F_p - rhoE_F_m )/delta_x;
                /// y-direction j+1/2
                index_L = j;                           index_R = j + 1;
                rho_L   = rho_field[I1D(i,index_L,k)]; rho_R   = rho_field[I1D(i,index_R,k)];
                u_L     = u_field[I1D(i,index_L,k)];   u_R     = u_field[I1D(i,index_R,k)];
                v_L     = v_field[I1D(i,index_L,k)];   v_R     = v_field[I1D(i,index_R,k)];
                w_L     = w_field[I1D(i,index_L,k)];   w_R     = w_field[I1D(i,index_R,k)];
                E_L     = E_field[I1D(i,index_L,k)];   E_R     = E_field[I1D(i,index_R,k)];
                P_L     = P_field[I1D(i,index_L,k)];   P_R     = P_field[I1D(i,index_R,k)];
                a_L     = sos_field[I1D(i,index_L,k)]; a_R     = sos_field[I1D(i,index_R,k)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_L*v_L;
                rho_F_R  = rho_R*v_R;
                rho_U_L  = rho_L;
                rho_U_R  = rho_R;
                rho_F_p  = riemann_solver->calculateIntercellFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 2;
                rhou_F_L = rho_L*v_L*u_L;
                rhou_F_R = rho_R*v_R*u_R;
                rhou_U_L = rho_L*u_L;
                rhou_U_R = rho_R*u_R;
                rhou_F_p = riemann_solver->calculateIntercellFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 1;
                rhov_F_L = rho_L*v_L*v_L + P_L;
                rhov_F_R = rho_R*v_R*v_R + P_R;
                rhov_U_L = rho_L*v_L;
                rhov_U_R = rho_R*v_R;
                rhov_F_p = riemann_solver->calculateIntercellFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_L*v_L*w_L;
                rhow_F_R = rho_R*v_R*w_R;
                rhow_U_L = rho_L*w_L;
                rhow_U_R = rho_R*w_R;
                rhow_F_p = riemann_solver->calculateIntercellFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_L*v_L*E_L + v_L*P_L;
                rhoE_F_R = rho_R*v_R*E_R + v_R*P_R;
                rhoE_U_L = rho_L*E_L;
                rhoE_U_R = rho_R*E_R;
                rhoE_F_p = riemann_solver->calculateIntercellFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// y-direction j-1/2
                index_L = j - 1;                       index_R = j;
                rho_L   = rho_field[I1D(i,index_L,k)]; rho_R   = rho_field[I1D(i,index_R,k)];
                u_L     = u_field[I1D(i,index_L,k)];   u_R     = u_field[I1D(i,index_R,k)];
                v_L     = v_field[I1D(i,index_L,k)];   v_R     = v_field[I1D(i,index_R,k)];
                w_L     = w_field[I1D(i,index_L,k)];   w_R     = w_field[I1D(i,index_R,k)];
                E_L     = E_field[I1D(i,index_L,k)];   E_R     = E_field[I1D(i,index_R,k)];
                P_L     = P_field[I1D(i,index_L,k)];   P_R     = P_field[I1D(i,index_R,k)];
                a_L     = sos_field[I1D(i,index_L,k)]; a_R     = sos_field[I1D(i,index_R,k)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_L*v_L;
                rho_F_R  = rho_R*v_R;
                rho_U_L  = rho_L;
                rho_U_R  = rho_R;
                rho_F_m  = riemann_solver->calculateIntercellFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 2;
                rhou_F_L = rho_L*v_L*u_L;
                rhou_F_R = rho_R*v_R*u_R;
                rhou_U_L = rho_L*u_L;
                rhou_U_R = rho_R*u_R;
                rhou_F_m = riemann_solver->calculateIntercellFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 1;
                rhov_F_L = rho_L*v_L*v_L + P_L;
                rhov_F_R = rho_R*v_R*v_R + P_R;
                rhov_U_L = rho_L*v_L;
                rhov_U_R = rho_R*v_R;
                rhov_F_m = riemann_solver->calculateIntercellFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                rhow_F_L = rho_L*v_L*w_L;
                rhow_F_R = rho_R*v_R*w_R;
                rhow_U_L = rho_L*w_L;
                rhow_U_R = rho_R*w_R;
                rhow_F_m = riemann_solver->calculateIntercellFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_L*v_L*E_L + v_L*P_L;
                rhoE_F_R = rho_R*v_R*E_R + v_R*P_R;
                rhoE_U_L = rho_L*E_L;
                rhoE_U_R = rho_R*E_R;
                rhoE_F_m = riemann_solver->calculateIntercellFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// Fluxes y-direction
                rho_inv_flux[I1D(i,j,k)]  += ( rho_F_p - rho_F_m )/delta_y;
                rhou_inv_flux[I1D(i,j,k)] += ( rhou_F_p - rhou_F_m )/delta_y;
                rhov_inv_flux[I1D(i,j,k)] += ( rhov_F_p - rhov_F_m )/delta_y;
                rhow_inv_flux[I1D(i,j,k)] += ( rhow_F_p - rhow_F_m )/delta_y;
                rhoE_inv_flux[I1D(i,j,k)] += ( rhoE_F_p - rhoE_F_m )/delta_y;
                /// z-direction k+1/2
                index_L = k;                           index_R = k + 1;
                rho_L   = rho_field[I1D(i,j,index_L)]; rho_R   = rho_field[I1D(i,j,index_R)];
                u_L     = u_field[I1D(i,j,index_L)];   u_R     = u_field[I1D(i,j,index_R)];
                v_L     = v_field[I1D(i,j,index_L)];   v_R     = v_field[I1D(i,j,index_R)];
                w_L     = w_field[I1D(i,j,index_L)];   w_R     = w_field[I1D(i,j,index_R)];
                E_L     = E_field[I1D(i,j,index_L)];   E_R     = E_field[I1D(i,j,index_R)];
                P_L     = P_field[I1D(i,j,index_L)];   P_R     = P_field[I1D(i,j,index_R)];
                a_L     = sos_field[I1D(i,j,index_L)]; a_R     = sos_field[I1D(i,j,index_R)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_L*w_L;
                rho_F_R  = rho_R*w_R;
                rho_U_L  = rho_L;
                rho_U_R  = rho_R;
                rho_F_p  = riemann_solver->calculateIntercellFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 3;
                rhou_F_L = rho_L*w_L*u_L;
                rhou_F_R = rho_R*w_R*u_R;
                rhou_U_L = rho_L*u_L;
                rhou_U_R = rho_R*u_R;
                rhou_F_p = riemann_solver->calculateIntercellFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_L*w_L*v_L;
                rhov_F_R = rho_R*w_R*v_R;
                rhov_U_L = rho_L*v_L;
                rhov_U_R = rho_R*v_R;
                rhov_F_p = riemann_solver->calculateIntercellFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 1;
                rhow_F_L = rho_L*w_L*w_L + P_L;
                rhow_F_R = rho_R*w_R*w_R + P_R;
                rhow_U_L = rho_L*w_L;
                rhow_U_R = rho_R*w_R;
                rhow_F_p = riemann_solver->calculateIntercellFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_L*w_L*E_L + w_L*P_L;
                rhoE_F_R = rho_R*w_R*E_R + w_R*P_R;
                rhoE_U_L = rho_L*E_L;
                rhoE_U_R = rho_R*E_R;
                rhoE_F_p = riemann_solver->calculateIntercellFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// z-direction k-1/2
                index_L = k - 1;                       index_R = k;
                rho_L   = rho_field[I1D(i,j,index_L)]; rho_R   = rho_field[I1D(i,j,index_R)];
                u_L     = u_field[I1D(i,j,index_L)];   u_R     = u_field[I1D(i,j,index_R)];
                v_L     = v_field[I1D(i,j,index_L)];   v_R     = v_field[I1D(i,j,index_R)];
                w_L     = w_field[I1D(i,j,index_L)];   w_R     = w_field[I1D(i,j,index_R)];
                E_L     = E_field[I1D(i,j,index_L)];   E_R     = E_field[I1D(i,j,index_R)];
                P_L     = P_field[I1D(i,j,index_L)];   P_R     = P_field[I1D(i,j,index_R)];
                a_L     = sos_field[I1D(i,j,index_L)]; a_R     = sos_field[I1D(i,j,index_R)];
                /// rho
                var_type = 0;
                rho_F_L  = rho_L*w_L;
                rho_F_R  = rho_R*w_R;
                rho_U_L  = rho_L;
                rho_U_R  = rho_R;
                rho_F_m  = riemann_solver->calculateIntercellFlux( rho_F_L, rho_F_R, rho_U_L, rho_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 3;
                rhou_F_L = rho_L*w_L*u_L;
                rhou_F_R = rho_R*w_R*u_R;
                rhou_U_L = rho_L*u_L;
                rhou_U_R = rho_R*u_R;
                rhou_F_m = riemann_solver->calculateIntercellFlux( rhou_F_L, rhou_F_R, rhou_U_L, rhou_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                rhov_F_L = rho_L*w_L*v_L;
                rhov_F_R = rho_R*w_R*v_R;
                rhov_U_L = rho_L*v_L;
                rhov_U_R = rho_R*v_R;
                rhov_F_m = riemann_solver->calculateIntercellFlux( rhov_F_L, rhov_F_R, rhov_U_L, rhov_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhow
                var_type = 1;
                rhow_F_L = rho_L*w_L*w_L + P_L;
                rhow_F_R = rho_R*w_R*w_R + P_R;
                rhow_U_L = rho_L*w_L;
                rhow_U_R = rho_R*w_R;
                rhow_F_m = riemann_solver->calculateIntercellFlux( rhow_F_L, rhow_F_R, rhow_U_L, rhow_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                rhoE_F_L = rho_L*w_L*E_L + w_L*P_L;
                rhoE_F_R = rho_R*w_R*E_R + w_R*P_R;
                rhoE_U_L = rho_L*E_L;
                rhoE_U_R = rho_R*E_R;
                rhoE_F_m = riemann_solver->calculateIntercellFlux( rhoE_F_L, rhoE_F_R, rhoE_U_L, rhoE_U_R, rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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
    //rho_inv_flux.update();
    //rhou_inv_flux.update();
    //rhov_inv_flux.update();
    //rhow_inv_flux.update();
    //rhoE_inv_flux.update();

};

void FlowSolverRHEA::calculateViscousFluxes() {

    /// Second-order central finite differences for derivatives:
    /// P. Moin.
    /// Fundamentals of engineering numerical analysis.
    /// Cambridge University Press, 2010.

    /// Inner points: rhou, rhov, rhow and rhoE
    double delta_x, delta_y, delta_z;
    double d_u_x, d_u_y, d_u_z, d_v_x, d_v_y, d_v_z, d_w_x, d_w_y, d_w_z;
    double div_uvw, tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;
    double div_tau_x, div_tau_y, div_tau_z;
    double div_q, div_uvw_tau; 
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = mesh->x[i+1] - mesh->x[i-1]; 
                delta_y = mesh->y[j+1] - mesh->y[j-1]; 
                delta_z = mesh->z[k+1] - mesh->z[k-1];
                /// Velocity derivatives
                d_u_x = ( u_field[I1D(i+1,j,k)] - u_field[I1D(i-1,j,k)] )/delta_x;
                d_u_y = ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j-1,k)] )/delta_y;
                d_u_z = ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k-1)] )/delta_z;
                d_v_x = ( v_field[I1D(i+1,j,k)] - v_field[I1D(i-1,j,k)] )/delta_x;
                d_v_y = ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j-1,k)] )/delta_y;
                d_v_z = ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k-1)] )/delta_z;
                d_w_x = ( w_field[I1D(i+1,j,k)] - w_field[I1D(i-1,j,k)] )/delta_x;
                d_w_y = ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j-1,k)] )/delta_y;
                d_w_z = ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k-1)] )/delta_z;
                /// Divergence of velocity
                div_uvw = d_u_x + d_v_y + d_w_z;
                /// Viscous stresses ( symmetric tensor )
                tau_xx = 2.0*mu_field[I1D(i,j,k)]*( d_u_x - ( div_uvw/3.0 ) );
                tau_xy = mu_field[I1D(i,j,k)]*( d_u_y + d_v_x );
                tau_xz = mu_field[I1D(i,j,k)]*( d_u_z + d_w_x );
                tau_yy = 2.0*mu_field[I1D(i,j,k)]*( d_v_y - ( div_uvw/3.0 ) );
                tau_yz = mu_field[I1D(i,j,k)]*( d_v_z + d_w_y );
                tau_zz = 2.0*mu_field[I1D(i,j,k)]*( d_w_z - ( div_uvw/3.0 ) );
                /// Divergence of viscous stresses
                div_tau_x = mu_field[I1D(i,j,k)]*( ( 2.0/delta_x )*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] )
                                                                   - ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) )
                                                 + ( 2.0/delta_y )*( ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] )
                                                                   - ( u_field[I1D(i,j,k)] - u_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) )
                                                 + ( 2.0/delta_z )*( ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] )
                                                                   - ( u_field[I1D(i,j,k)] - u_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) )
              + ( 1.0/3.0 )*mu_field[I1D(i,j,k)]*( ( 2.0/delta_x )*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] )
                                                                   - ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) )
                                                 + ( 1.0/delta_x )*( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i+1,j-1,k)] )/delta_y
                                                                   - ( v_field[I1D(i-1,j+1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_y )
                                                 + ( 1.0/delta_x )*( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i+1,j,k-1)] )/delta_z
                                                                   - ( w_field[I1D(i-1,j,k+1)] - w_field[I1D(i-1,j,k-1)] )/delta_z ) );
                div_tau_y = mu_field[I1D(i,j,k)]*( ( 2.0/delta_x )*( ( v_field[I1D(i+1,j,k)] - v_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] )
                                                                   - ( v_field[I1D(i,j,k)] - v_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) )
                                                 + ( 2.0/delta_y )*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] )
                                                                   - ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) )
                                                 + ( 2.0/delta_z )*( ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] )
                                                                   - ( v_field[I1D(i,j,k)] - v_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) )
              + ( 1.0/3.0 )*mu_field[I1D(i,j,k)]*( ( 1.0/delta_y )*( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i-1,j+1,k)] )/delta_x
                                                                   - ( u_field[I1D(i+1,j-1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_x )
                                                 + ( 2.0/delta_y )*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] )
                                                                   - ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) )
                                                 + ( 1.0/delta_y )*( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j+1,k-1)] )/delta_z
                                                                   - ( w_field[I1D(i,j-1,k+1)] - w_field[I1D(i,j-1,k-1)] )/delta_z ) ); 
                div_tau_z = mu_field[I1D(i,j,k)]*( ( 2.0/delta_x )*( ( w_field[I1D(i+1,j,k)] - w_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] )
                                                                   - ( w_field[I1D(i,j,k)] - w_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) )
                                                 + ( 2.0/delta_y )*( ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] )
                                                                   - ( w_field[I1D(i,j,k)] - w_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) )
                                                 + ( 2.0/delta_z )*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] )
                                                                   - ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) )
              + ( 1.0/3.0 )*mu_field[I1D(i,j,k)]*( ( 1.0/delta_z )*( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i-1,j,k+1)] )/delta_x
                                                                   - ( u_field[I1D(i+1,j,k-1)] - u_field[I1D(i-1,j,k-1)] )/delta_x )
                                                 + ( 1.0/delta_z )*( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j-1,k+1)] )/delta_y
                                                                   - ( v_field[I1D(i,j+1,k-1)] - v_field[I1D(i,j-1,k-1)] )/delta_y )
                                                 + ( 2.0/delta_z )*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] )
                                                                   - ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) );
                /// Fourier term
                div_q = kappa_field[I1D(i,j,k)]*( ( 2.0/delta_x )*( ( T_field[I1D(i+1,j,k)] - T_field[I1D(i,j,k)] )/( mesh->x[i+1] - mesh->x[i] )
                                                                  - ( T_field[I1D(i,j,k)] - T_field[I1D(i-1,j,k)] )/( mesh->x[i] - mesh->x[i-1] ) )
                                                + ( 2.0/delta_y )*( ( T_field[I1D(i,j+1,k)] - T_field[I1D(i,j,k)] )/( mesh->y[j+1] - mesh->y[j] )
                                                                  - ( T_field[I1D(i,j,k)] - T_field[I1D(i,j-1,k)] )/( mesh->y[j] - mesh->y[j-1] ) )
                                                + ( 2.0/delta_z )*( ( T_field[I1D(i,j,k+1)] - T_field[I1D(i,j,k)] )/( mesh->z[k+1] - mesh->z[k] )
                                                                  - ( T_field[I1D(i,j,k)] - T_field[I1D(i,j,k-1)] )/( mesh->z[k] - mesh->z[k-1] ) ) );
                /// Work of viscous stresses
                div_uvw_tau = u_field[I1D(i,j,k)]*div_tau_x + v_field[I1D(i,j,k)]*div_tau_y + w_field[I1D(i,j,k)]*div_tau_z
                            + tau_xx*d_u_x + tau_xy*d_u_y + tau_xz*d_u_z
                            + tau_xy*d_v_x + tau_yy*d_v_y + tau_yz*d_v_z
                            + tau_xz*d_w_x + tau_yz*d_w_y + tau_zz*d_w_z;
                /// Viscous fluxes
                rhou_vis_flux[I1D(i,j,k)] = div_tau_x;
                rhov_vis_flux[I1D(i,j,k)] = div_tau_y;
                rhow_vis_flux[I1D(i,j,k)] = div_tau_z;
                rhoE_vis_flux[I1D(i,j,k)] = div_q + div_uvw_tau;
            }
        }
    }

    /// Update halo values
    //rhou_vis_flux.update();
    //rhov_vis_flux.update();
    //rhow_vis_flux.update();
    //rhoE_vis_flux.update();

};

void FlowSolverRHEA::sumInviscidViscousFluxesSourceTerms(const int &rk_step) {

    /// First Runge-Kutta step
    if(rk_step == 1) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        double f_rhouvw = 0.0;
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    /// Work of momentum sources
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
        //rho_rk1_flux.update();
        //rhou_rk1_flux.update();
        //rhov_rk1_flux.update();
        //rhow_rk1_flux.update();
        //rhoE_rk1_flux.update();

    }

    /// Second Runge-Kutta step
    if(rk_step == 2) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        double f_rhouvw = 0.0;
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    /// Work of momentum sources
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
        //rho_rk2_flux.update();
        //rhou_rk2_flux.update();
        //rhov_rk2_flux.update();
        //rhow_rk2_flux.update();
        //rhoE_rk2_flux.update();

    }

    /// Third Runge-Kutta step
    if(rk_step == 3) {

        /// Inner points: rho, rhou, rhov, rhow and rhoE
        double f_rhouvw = 0.0;
        for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                    /// Work of momentum sources
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
        //rho_rk3_flux.update();
        //rhou_rk3_flux.update();
        //rhov_rk3_flux.update();
        //rhow_rk3_flux.update();
        //rhoE_rk3_flux.update();

    }

};

void FlowSolverRHEA::timeAdvanceConservedVariables(const int &rk_step) {

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
    //rho_field.update();
    //rhou_field.update();
    //rhov_field.update();
    //rhow_field.update();
    //rhoE_field.update();

};

void FlowSolverRHEA::outputCurrentStateData() {

    /// Write to file current solver state, time and time iteration
    writer_reader->setAttribute( "Time", current_time );
    writer_reader->setAttribute( "Iteration", current_time_iter );
    writer_reader->setAttribute( "AveragingTime", averaging_time );
    writer_reader->write( current_time_iter );

};

void FlowSolverRHEA::updateTimeAveragedQuantities() {

    /// All (inner, boundary & halo) points: time-averaged and root-mean-square-fluctuation quantities 
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                /// Time-averaged quantities
                avg_rho_field[I1D(i,j,k)]  = ( avg_rho_field[I1D(i,j,k)]*averaging_time + rho_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_rhou_field[I1D(i,j,k)] = ( avg_rhou_field[I1D(i,j,k)]*averaging_time + rhou_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_rhov_field[I1D(i,j,k)] = ( avg_rhov_field[I1D(i,j,k)]*averaging_time + rhov_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_rhow_field[I1D(i,j,k)] = ( avg_rhow_field[I1D(i,j,k)]*averaging_time + rhow_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_rhoE_field[I1D(i,j,k)] = ( avg_rhoE_field[I1D(i,j,k)]*averaging_time + rhoE_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_P_field[I1D(i,j,k)]    = ( avg_P_field[I1D(i,j,k)]*averaging_time + P_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_T_field[I1D(i,j,k)]    = ( avg_T_field[I1D(i,j,k)]*averaging_time + T_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                avg_sos_field[I1D(i,j,k)]  = ( avg_sos_field[I1D(i,j,k)]*averaging_time + sos_field[I1D(i,j,k)]*delta_t )/( averaging_time + delta_t );
                /// Root-mean-square-fluctuation quantities
                rmsf_rho_field[I1D(i,j,k)]  = sqrt( ( pow( rmsf_rho_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( rho_field[I1D(i,j,k)] - avg_rho_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_rhou_field[I1D(i,j,k)] = sqrt( ( pow( rmsf_rhou_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( rhou_field[I1D(i,j,k)] - avg_rhou_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_rhov_field[I1D(i,j,k)] = sqrt( ( pow( rmsf_rhov_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( rhov_field[I1D(i,j,k)] - avg_rhov_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_rhow_field[I1D(i,j,k)] = sqrt( ( pow( rmsf_rhow_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( rhow_field[I1D(i,j,k)] - avg_rhow_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_rhoE_field[I1D(i,j,k)] = sqrt( ( pow( rmsf_rhoE_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( rhoE_field[I1D(i,j,k)] - avg_rhoE_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_P_field[I1D(i,j,k)]    = sqrt( ( pow( rmsf_P_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( P_field[I1D(i,j,k)] - avg_P_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_T_field[I1D(i,j,k)]    = sqrt( ( pow( rmsf_T_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( T_field[I1D(i,j,k)] - avg_T_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
                rmsf_sos_field[I1D(i,j,k)]  = sqrt( ( pow( rmsf_sos_field[I1D(i,j,k)], 2.0 )*averaging_time + pow( sos_field[I1D(i,j,k)] - avg_sos_field[I1D(i,j,k)], 2.0 )*delta_t )/( averaging_time + delta_t ) );
            }
        }
    }

    /// Update averaging time
    averaging_time += delta_t;

    /// Update halo values
    avg_rho_field.update();
    avg_rhou_field.update();
    avg_rhov_field.update();
    avg_rhow_field.update();
    avg_rhoE_field.update();
    avg_P_field.update();
    avg_T_field.update();
    avg_sos_field.update();
    rmsf_rho_field.update();
    rmsf_rhou_field.update();
    rmsf_rhov_field.update();
    rmsf_rhow_field.update();
    rmsf_rhoE_field.update();
    rmsf_P_field.update();
    rmsf_T_field.update();
    rmsf_sos_field.update();

};

void FlowSolverRHEA::execute() {

    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /// Set output (cout) precision
    cout.precision( cout_presicion );

    /// Start RHEA simulation
    if( my_rank == 0 ) cout << "RHEA: START SIMULATION" << endl;

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

        /// Start timer: output_solver_state
        timers->start( "output_solver_state" );

        /// Print time iteration information (if criterion satisfied)
        if( ( current_time_iter%print_frequency_iter == 0 ) and ( my_rank == 0 ) ) {
            cout << "Time iteration " << current_time_iter << ": " 
                 << "time = " << scientific << current_time << " [s], "
                 << "time-step = " << scientific << delta_t << " [s]" << endl;
        }

        /// Output current state data to file (if criterion satisfied)
        if( current_time_iter%output_frequency_iter == 0 ) this->outputCurrentStateData();

        /// Stop timer: output_solver_state
        timers->stop( "output_solver_state" );

        /// Start timer: rk_iteration_loop
        timers->start( "rk_iteration_loop" );

        /// Runge-Kutta time-integration steps
        for(int rk_step = 1; rk_step <= rk_order; rk_step++) {

            /// Start timer: calculate_thermophysical_properties
            timers->start( "calculate_thermophysical_properties" );

            /// Calculate transport coefficients
            this->calculateTransportCoefficients();

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

            /// Start timer: sum_inviscid_viscous_fluxes_source_terms
            timers->start( "sum_inviscid_viscous_fluxes_source_terms" );

            /// Sum inviscid & viscous fluxes, and source terms (right-hand side)
            this->sumInviscidViscousFluxesSourceTerms(rk_step);

            /// Stop timer: sum_inviscid_viscous_fluxes_source_terms
            timers->stop( "sum_inviscid_viscous_fluxes_source_terms" );

            /// Start timer: time_advance_conserved_variables
            timers->start( "time_advance_conserved_variables" );

            /// Advance conserved variables in time
            this->timeAdvanceConservedVariables(rk_step);

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
    if( my_rank == 0 ) cout << "RHEA: END SIMULATION" << endl;

};


////////// BaseRiemannSolver CLASS //////////

BaseRiemannSolver::BaseRiemannSolver() {};
        
BaseRiemannSolver::~BaseRiemannSolver() {};

void BaseRiemannSolver::calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R) {

#if _PRESSURE_BASED_WAVE_SPEED_ESTIMATES_
    /// Pressure-based wave speed estimates: ... recommended by E. F. Toro, but not sure if applicable to non-ideal gas themodynamics
    /// E. F. Toro, M. Spruce, W. Speares.
    /// Restoration of the contact surface in the HLL-Riemann solver.
    /// Shock Waves, 4, 25-34, 1994.

    double P_bar   = 0.5*( P_L + P_R );
    double rho_bar = 0.5*( rho_L + rho_R );
    double gamma   = thermodynamics->calculateHeatCapacitiesRatio( P_bar, rho_bar );
    double a_bar   = 0.5*( a_L + a_R );
    double P_pvrs  = 0.5*( P_L + P_R ) - 0.5*( u_R - u_L )*rho_bar*a_bar;
    double P_star  = max( 0.0, P_pvrs );
    double q_L     = 1.0;
    if(P_star > P_L) q_L = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_L ) - 1.0 ) );
    double q_R     = 1.0;
    if(P_star > P_R) q_R = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_R ) - 1.0 ) );
    S_L = u_L - a_L*q_L;
    S_R = u_R + a_R*q_R;
#else
    /// Direct wave speed estimates:
    /// B. Einfeldt.
    /// On Godunov-type methods for gas dynamics.
    /// SIAM Journal on Numerical Analysis, 25, 294-318, 1988.

    double hat_u = ( u_L*sqrt( rho_L ) + u_R*sqrt( rho_R ) )/( sqrt( rho_L ) + sqrt( rho_R ) );
    double hat_a = sqrt( ( ( a_L*a_L*sqrt( rho_L ) + a_R*a_R*sqrt( rho_R ) )/( sqrt( rho_L ) + sqrt( rho_R ) ) ) + 0.5*( ( sqrt( rho_L )*sqrt( rho_R ) )/( ( sqrt( rho_L ) + sqrt( rho_R ) )*( sqrt( rho_L ) + sqrt( rho_R ) ) ) )*( u_R - u_L )*( u_R - u_L ) );

    S_L = min( u_L - a_L, hat_u - hat_a );
    S_R = max( u_R + a_R, hat_u + hat_a );
#endif

};


////////// CentralFluxApproximateRiemannSolver CLASS //////////

CentralFluxApproximateRiemannSolver::CentralFluxApproximateRiemannSolver() : BaseRiemannSolver() {};

CentralFluxApproximateRiemannSolver::~CentralFluxApproximateRiemannSolver() {};

double CentralFluxApproximateRiemannSolver::calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// Central scheme obtained from a central differencing of the first derivative of the flux term:

    double F = 0.5*( F_L + F_R );

    return( F );

};


////////// MurmanRoeFluxApproximateRiemannSolver CLASS //////////

MurmanRoeFluxApproximateRiemannSolver::MurmanRoeFluxApproximateRiemannSolver() : BaseRiemannSolver() {};

MurmanRoeFluxApproximateRiemannSolver::~MurmanRoeFluxApproximateRiemannSolver() {};

double MurmanRoeFluxApproximateRiemannSolver::calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// Murman-Roe Riemman solver:
    /// P. L. Roe.
    /// Approximate Riemann solvers, parameter vectors and difference schemes.
    /// Journal of Computational Physics, 43, 357-372, 1981.

    /// Wave speed
    double S = abs( ( F_L - F_R )/( U_L - U_R + epsilon ) );

    /// Conservative + dissipative flux form
    double F = 0.5*( F_L + F_R ) - 0.5*S*( U_R - U_L );

    return( F );

};


////////// KgpFluxApproximateRiemannSolver CLASS //////////

KgpFluxApproximateRiemannSolver::KgpFluxApproximateRiemannSolver() : BaseRiemannSolver() {};

KgpFluxApproximateRiemannSolver::~KgpFluxApproximateRiemannSolver() {};

double KgpFluxApproximateRiemannSolver::calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// Kennedy, Gruber & Pirozzoli (KGP) scheme:
    /// G. Coppola, F. Capuano, S. Pirozzoli, L. de Luca.
    /// Numerically stable formulations of convective terms for turbulent compressible flows.
    /// Journal of Computational Physics, 382, 86-104, 2019.

    double F = ( 1.0/8.0 )*( rho_L + rho_R )*( u_L + u_R );
    if( var_type == 0 ) {
        F *= 1.0 + 1.0;
    } else if ( var_type == 1 ) {
        F *= u_L + u_R;
        F += ( 1.0/2.0 )*( P_L + P_R );
    } else if ( var_type == 2 ) {
        F *= v_L + v_R;
    } else if ( var_type == 3 ) {
        F *= w_L + w_R;
    } else if ( var_type == 4 ) {
        F *= E_L + E_R;
        F += ( 1.0/4.0 )*( u_L + u_R )*( P_L + P_R );
    }

    return( F );

};


////////// HllApproximateRiemannSolver CLASS //////////

HllApproximateRiemannSolver::HllApproximateRiemannSolver() : BaseRiemannSolver() {};

HllApproximateRiemannSolver::~HllApproximateRiemannSolver() {};

double HllApproximateRiemannSolver::calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// Harten-Lax-van Leer (HLL) Riemman solver:
    /// A. Harten, P. D. Lax, B. van Leer.
    /// On upstream differencing and Godunov-type schemes for hyperbolic conservation laws.
    /// SIAM Review, 25, 35-61, 1983.

    double S_L, S_R;
    this->calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );

    double F = 0.0;
    if( 0.0 <= S_L ) {
        F = F_L;
    } else if( 0.0 >= S_R ) {
        F = F_R;
    } else {
        F = ( S_R*F_L - S_L*F_R + S_L*S_R*( U_R - U_L ) )/( S_R - S_L );
    }

    return( F );

};


////////// HllcApproximateRiemannSolver CLASS //////////

HllcApproximateRiemannSolver::HllcApproximateRiemannSolver() : BaseRiemannSolver() {};

HllcApproximateRiemannSolver::~HllcApproximateRiemannSolver() {};

double HllcApproximateRiemannSolver::calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// Harten-Lax-van Leer-Contact (HLLC) Riemman solver:
    /// E. F. Toro, M. Spruce, W. Speares.
    /// Restoration of the contact surface in the HLL-Riemann solver.
    /// Shock Waves, 4, 25-34, 1994.

    double S_L, S_R;
    this->calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );

    double S_star   = ( P_R - P_L + rho_L*u_L*( S_L - u_L ) - rho_R*u_R*( S_R - u_R ) )/( rho_L*( S_L - u_L ) - rho_R*( S_R - u_R ) );
    double U_star_L = rho_L*( ( S_L - u_L )/( S_L - S_star ) );
    double U_star_R = rho_R*( ( S_R - u_R )/( S_R - S_star ) );
    if( var_type == 0 ) {
        U_star_L *= 1.0;
        U_star_R *= 1.0;       
    } else if( var_type == 1 ) {
        U_star_L *= S_star;
        U_star_R *= S_star;
    } else if( var_type == 2 ) {
        U_star_L *= v_L;
        U_star_R *= v_R;
    } else if( var_type == 3 ) {
        U_star_L *= w_L;
        U_star_R *= w_R;
    } else if( var_type == 4 ) {
        U_star_L *= ( E_L + ( S_star - u_L )*( S_star + P_L/( rho_L*( S_L - u_L ) ) ) );
        U_star_R *= ( E_R + ( S_star - u_R )*( S_star + P_R/( rho_R*( S_R - u_R ) ) ) );
    }

    double F = 0.0;
    if( 0.0 <= S_L ) {
        F = F_L;
    } else if( ( S_L <= 0.0 ) && ( 0.0 <= S_star ) ) {
        F = F_L + S_L*( U_star_L - U_L );
    } else if( ( S_star <= 0.0 ) && ( 0.0 <= S_R ) ) {
        F = F_R + S_R*( U_star_R - U_R );
    } else if( 0.0 >= S_R ) {
        F = F_R;
    }

    return( F );

};


////////// HllcPlusApproximateRiemannSolver CLASS //////////

HllcPlusApproximateRiemannSolver::HllcPlusApproximateRiemannSolver() : BaseRiemannSolver() {};

HllcPlusApproximateRiemannSolver::~HllcPlusApproximateRiemannSolver() {};

double HllcPlusApproximateRiemannSolver::calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    /// HLLC-type Riemann solver for all-speed flows:
    /// S. Chen, B. Lin, Y. Li, C. Yan.
    /// HLLC+: low-Mach shock-stable HLLC-type Riemann solver for all-speed flows.
    /// SIAM Journal of Scientific Computing, 4, B921-B950, 2020.

    double S_L, S_R;
    this->calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );

    double phi_L = rho_L*( S_L - u_L );
    double phi_R = rho_R*( S_R - u_R );
    double S_star   = ( P_R - P_L + phi_L*u_L - phi_R*u_R )/( phi_L - phi_R );
    double U_star_L = rho_L*( ( S_L - u_L )/( S_L - S_star ) );
    double U_star_R = rho_R*( ( S_R - u_R )/( S_R - S_star ) );
    double M        = min( 1.0, max( ( 1.0/a_L )*sqrt( u_L*u_L + v_L*v_L + w_L*w_L ), ( 1.0/a_R )*sqrt( u_R*u_R + v_R*v_R + w_R*w_R ) ) );
    double f_M      = M*sqrt( 4.0 + pow( 1.0 - M*M, 2.0 ) )/( 1.0 + M*M );
    double h        = min( P_L/P_R, P_R/P_L );
    double g        = 1.0 - pow( h, M );
    double A_p_L    = ( phi_L*phi_R )/( phi_R - phi_L );
    double A_p_R    = ( phi_L*phi_R )/( phi_R - phi_L );
    if( var_type == 0 ) {
        U_star_L *= 1.0;
        U_star_R *= 1.0;       
        A_p_L    *= 0.0;
        A_p_R    *= 0.0;
    } else if( var_type == 1 ) {
        U_star_L *= S_star;
        U_star_R *= S_star;
        A_p_L    *= ( f_M - 1.0 )*( u_R - u_L );
        A_p_R    *= ( f_M - 1.0 )*( u_R - u_L );
    } else if( var_type == 2 ) {
        U_star_L *= v_L;
        U_star_R *= v_R;
        A_p_L    *= ( S_L/( S_L - S_star ) )*g*( v_R - v_L );
        A_p_R    *= ( S_R/( S_R - S_star ) )*g*( v_R - v_L );
    } else if( var_type == 3 ) {
        U_star_L *= w_L;
        U_star_R *= w_R;
        A_p_L    *= ( S_L/( S_L - S_star ) )*g*( w_R - w_L );
        A_p_R    *= ( S_R/( S_R - S_star ) )*g*( w_R - w_L );
    } else if( var_type == 4 ) {
        U_star_L *= ( E_L + ( S_star - u_L )*( S_star + P_L/( rho_L*( S_L - u_L ) ) ) );
        U_star_R *= ( E_R + ( S_star - u_R )*( S_star + P_R/( rho_R*( S_R - u_R ) ) ) );
        A_p_L    *= ( f_M - 1.0 )*( u_R - u_L )*S_star;
        A_p_R    *= ( f_M - 1.0 )*( u_R - u_L )*S_star;
    }
    double F_star_L = F_L + S_L*( U_star_L - U_L ) + A_p_L;
    double F_star_R = F_R + S_R*( U_star_R - U_R ) + A_p_R;

    double F = 0.0;
    if( 0.0 <= S_L ) {
        F = F_L;
    } else if( ( S_L <= 0.0 ) && ( 0.0 <= S_star ) ) {
        F = F_star_L;
    } else if( ( S_star <= 0.0 ) && ( 0.0 <= S_R ) ) {
        F = F_star_R;
    } else if( 0.0 >= S_R ) {
        F = F_R;
    }

    return( F );

};

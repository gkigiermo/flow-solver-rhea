#include "FlowSolverRHEA.hpp"

using namespace std;


////////// COMPILATION DIRECTIVES //////////
#define _PRESSURE_BASED_WAVE_SPEED_ESTIMATES_ 0		/// Select approach for estimating wave speeds
#define _OPENACC_MANUAL_DATA_MOVEMENT_ 1		/// Select manual/managed data movement for OpenACC

////////// FIXED PARAMETERS //////////
const double epsilon     = 1.0e-15;			/// Small epsilon number (fixed)
const double pi          = 2.0*asin(1.0);		/// pi number (fixed)
const int cout_precision = 5;		                /// Output precision (fixed)
int riemann_solver_scheme_index;			/// Index of Riemann solver scheme


////////// PRAGMA DIRECTIVES //////////
#pragma acc declare create( riemann_solver_scheme_index )


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
    } else if( transport_coefficients_model == "LOW_PRESSURE_GAS" ) {
        transport_coefficients = new LowPressureGasTransportCoefficients( configuration_file );
    } else if( transport_coefficients_model == "HIGH_PRESSURE" ) {
        transport_coefficients = new HighPressureTransportCoefficients( configuration_file );
    } else {
        cout << "Transport coefficients model not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /// Construct (initialize) Riemann solver
//    if( riemann_solver_scheme == "DIVERGENCE" ) {
//        riemann_solver = new DivergenceFluxApproximateRiemannSolver();
//    } else if( riemann_solver_scheme == "MURMAN-ROE" ) {
//        riemann_solver = new MurmanRoeFluxApproximateRiemannSolver();
//    } else if( riemann_solver_scheme == "KGP" ) {
//        riemann_solver = new KgpFluxApproximateRiemannSolver();
//    } else if( riemann_solver_scheme == "SHIMA" ) {
//        riemann_solver = new ShimaFluxApproximateRiemannSolver();
//    } else if( riemann_solver_scheme == "HLL" ) {
//        riemann_solver = new HllApproximateRiemannSolver();
//    } else if( riemann_solver_scheme == "HLLC" ) {
//        riemann_solver = new HllcApproximateRiemannSolver();
//    } else if( riemann_solver_scheme == "HLLC+" ) {
//        riemann_solver = new HllcPlusApproximateRiemannSolver();
//    } else {
//        cout << "Riemann solver not available!" << endl;
//        MPI_Abort( MPI_COMM_WORLD, 1 );
//    }

    /// Riemann solver selector (index)
    if( riemann_solver_scheme == "KGP" ) {
        riemann_solver_scheme_index = 0;
    } else if( riemann_solver_scheme == "SHIMA" ) {
        riemann_solver_scheme_index = 1;
    } else if( riemann_solver_scheme == "DIVERGENCE" ) {
        riemann_solver_scheme_index = 2;
    } else if( riemann_solver_scheme == "MURMAN-ROE" ) {
        riemann_solver_scheme_index = 3;
    } else if( riemann_solver_scheme == "HLL" ) {
        riemann_solver_scheme_index = 4;
    } else if( riemann_solver_scheme == "HLLC" ) {
        riemann_solver_scheme_index = 5;
    } else if( riemann_solver_scheme == "HLLC+" ) {
        riemann_solver_scheme_index = 6;
    } else { 
        cout << "Riemann solver not available!" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    } 

    /// Construct (initialize) Runge-Kutta method
    if( runge_kutta_time_scheme == "RK1" ) {
        runge_kutta_method = new RungeKutta1Method();
        rk_number_stages = 1;
    } else if( runge_kutta_time_scheme == "SSP-RK2" ) {
        runge_kutta_method = new StrongStabilityPreservingRungeKutta2Method();
        rk_number_stages = 2;
    } else if( runge_kutta_time_scheme == "SSP-RK3" ) {
        runge_kutta_method = new StrongStabilityPreservingRungeKutta3Method();
        rk_number_stages = 3;
    } else {
        cout << "Runge-Kutta time scheme not available!" << endl;
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
    c_v_field.setTopology(topo,"c_v");
    c_p_field.setTopology(topo,"c_p");

    /// Set parallel topology of time-integration variables	
    rho_0_field.setTopology(topo,"rho_0");
    rhou_0_field.setTopology(topo,"rhou_0");
    rhov_0_field.setTopology(topo,"rhov_0");
    rhow_0_field.setTopology(topo,"rhow_0");
    rhoE_0_field.setTopology(topo,"rhoE_0");    
    P_0_field.setTopology(topo,"P_0");    

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
    work_vis_rhoe_flux.setTopology(topo,"work_vis_rhoe");

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
    avg_rhoP_field.setTopology(topo,"avg_rhoP");
    avg_rhoT_field.setTopology(topo,"avg_rhoT");
    avg_u_field.setTopology(topo,"avg_u");
    avg_v_field.setTopology(topo,"avg_v");
    avg_w_field.setTopology(topo,"avg_w");
    avg_E_field.setTopology(topo,"avg_E");
    avg_P_field.setTopology(topo,"avg_P");
    avg_T_field.setTopology(topo,"avg_T");
    avg_sos_field.setTopology(topo,"avg_sos");
    avg_mu_field.setTopology(topo,"avg_mu");
    avg_kappa_field.setTopology(topo,"avg_kappa");
    avg_c_v_field.setTopology(topo,"avg_c_v");
    avg_c_p_field.setTopology(topo,"avg_c_p");
    rmsf_rho_field.setTopology(topo,"rmsf_rho");
    rmsf_rhou_field.setTopology(topo,"rmsf_rhou");
    rmsf_rhov_field.setTopology(topo,"rmsf_rhov");
    rmsf_rhow_field.setTopology(topo,"rmsf_rhow");
    rmsf_rhoE_field.setTopology(topo,"rmsf_rhoE");
    rmsf_u_field.setTopology(topo,"rmsf_u");
    rmsf_v_field.setTopology(topo,"rmsf_v");
    rmsf_w_field.setTopology(topo,"rmsf_w");
    rmsf_E_field.setTopology(topo,"rmsf_E");
    rmsf_P_field.setTopology(topo,"rmsf_P");
    rmsf_T_field.setTopology(topo,"rmsf_T");
    rmsf_sos_field.setTopology(topo,"rmsf_sos");
    rmsf_mu_field.setTopology(topo,"rmsf_mu");
    rmsf_kappa_field.setTopology(topo,"rmsf_kappa");
    rmsf_c_v_field.setTopology(topo,"rmsf_c_v");
    rmsf_c_p_field.setTopology(topo,"rmsf_c_p");
    favre_uffuff_field.setTopology(topo,"favre_uffuff");
    favre_uffvff_field.setTopology(topo,"favre_uffvff");
    favre_uffwff_field.setTopology(topo,"favre_uffwff");
    favre_vffvff_field.setTopology(topo,"favre_vffvff");
    favre_vffwff_field.setTopology(topo,"favre_vffwff");
    favre_wffwff_field.setTopology(topo,"favre_wffwff");
    favre_uffEff_field.setTopology(topo,"favre_uffEff");
    favre_vffEff_field.setTopology(topo,"favre_vffEff");
    favre_wffEff_field.setTopology(topo,"favre_wffEff");

    /// Fill mesh x, y, z, delta_x, delta_y, delta_z fields
    this->fillMeshCoordinatesSizesFields();

    /// Construct (initialize) writer/reader
    char char_array[ output_data_file_name.length() + 1 ]; 
    strcpy( char_array,output_data_file_name.c_str() );
    writer_reader = new WriteReadHDF5( topo, char_array, generate_xdmf );
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
    writer_reader->addField(&c_v_field);
    writer_reader->addField(&c_p_field);
    writer_reader->addField(&avg_rho_field);
    writer_reader->addField(&avg_rhou_field);
    writer_reader->addField(&avg_rhov_field);
    writer_reader->addField(&avg_rhow_field);
    writer_reader->addField(&avg_rhoE_field);
    writer_reader->addField(&avg_rhoP_field);
    writer_reader->addField(&avg_rhoT_field);
    writer_reader->addField(&avg_u_field);
    writer_reader->addField(&avg_v_field);
    writer_reader->addField(&avg_w_field);
    writer_reader->addField(&avg_E_field);
    writer_reader->addField(&avg_P_field);
    writer_reader->addField(&avg_T_field);
    writer_reader->addField(&avg_sos_field);
    writer_reader->addField(&avg_mu_field);
    writer_reader->addField(&avg_kappa_field);
    writer_reader->addField(&avg_c_v_field);
    writer_reader->addField(&avg_c_p_field);
    writer_reader->addField(&rmsf_rho_field);
    writer_reader->addField(&rmsf_rhou_field);
    writer_reader->addField(&rmsf_rhov_field);
    writer_reader->addField(&rmsf_rhow_field);
    writer_reader->addField(&rmsf_rhoE_field);
    writer_reader->addField(&rmsf_u_field);
    writer_reader->addField(&rmsf_v_field);
    writer_reader->addField(&rmsf_w_field);
    writer_reader->addField(&rmsf_E_field);
    writer_reader->addField(&rmsf_P_field);
    writer_reader->addField(&rmsf_T_field);
    writer_reader->addField(&rmsf_sos_field);
    writer_reader->addField(&rmsf_mu_field);
    writer_reader->addField(&rmsf_kappa_field);
    writer_reader->addField(&rmsf_c_v_field);
    writer_reader->addField(&rmsf_c_p_field);
    writer_reader->addField(&favre_uffuff_field);
    writer_reader->addField(&favre_uffvff_field);
    writer_reader->addField(&favre_uffwff_field);
    writer_reader->addField(&favre_vffvff_field);
    writer_reader->addField(&favre_vffwff_field);
    writer_reader->addField(&favre_wffwff_field);
    writer_reader->addField(&favre_uffEff_field);
    writer_reader->addField(&favre_vffEff_field);
    writer_reader->addField(&favre_wffEff_field);

    /// Construct (initialize) timers
    timers = new ParallelTimer();
    timers->createTimer( "execute" );
    timers->createTimer( "time_iteration_loop" );
    timers->createTimer( "calculate_time_step" );
    timers->createTimer( "output_solver_state" );
    timers->createTimer( "rk_iteration_loop" );
    timers->createTimer( "calculate_thermophysical_properties" );
    timers->createTimer( "calculate_inviscid_fluxes" );
    timers->createTimer( "calculate_viscous_fluxes" );
    timers->createTimer( "calculate_source_terms" );
    timers->createTimer( "time_advance_conserved_variables" );
    timers->createTimer( "update_boundaries" );
    timers->createTimer( "conserved_to_primitive_variables" );
    timers->createTimer( "calculate_thermodynamics_from_primitive_variables" );
    timers->createTimer( "update_time_averaged_quantities" );
    timers->createTimer( "update_previous_state_conserved_variables" );

    /// Construct (initialize) temporal point probes
    TemporalPointProbe temporal_point_probe(mesh, topo);
    temporal_point_probes.resize( number_temporal_point_probes );
    for(int tpp = 0; tpp < number_temporal_point_probes; ++tpp) {
        /// Set parameters of temporal point probe
	temporal_point_probe.setPositionX( tpp_x_positions[tpp] );
	temporal_point_probe.setPositionY( tpp_y_positions[tpp] );
	temporal_point_probe.setPositionZ( tpp_z_positions[tpp] );
	temporal_point_probe.setOutputFileName( tpp_output_file_names[tpp] );
        /// Insert temporal point probe to vector
        temporal_point_probes[tpp] = temporal_point_probe;
	/// Locate closest grid point to probe
	temporal_point_probes[tpp].locateClosestGridPointToProbe();
    }	   

};

FlowSolverRHEA::~FlowSolverRHEA() {

    /// Free thermodynamics, transport_coefficients, riemann_solver, runge_kutta_method, mesh, topo, writer_reader and timers
    if( thermodynamics != NULL ) free( thermodynamics );
    if( transport_coefficients != NULL ) free( transport_coefficients );
    //if( riemann_solver != NULL ) free( riemann_solver );
    if( runge_kutta_method != NULL ) free( runge_kutta_method );
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
    num_grid_x                         = computational_parameters["num_grid_x"].as<int>();
    num_grid_y                         = computational_parameters["num_grid_y"].as<int>();
    num_grid_z                         = computational_parameters["num_grid_z"].as<int>();
    A_x                                = computational_parameters["A_x"].as<double>();
    A_y                                = computational_parameters["A_y"].as<double>();
    A_z                                = computational_parameters["A_z"].as<double>();
    CFL                                = computational_parameters["CFL"].as<double>();
    riemann_solver_scheme              = computational_parameters["riemann_solver_scheme"].as<string>();
    runge_kutta_time_scheme            = computational_parameters["runge_kutta_time_scheme"].as<string>();
    transport_pressure_scheme          = computational_parameters["transport_pressure_scheme"].as<bool>();
    artificial_compressibility_method  = computational_parameters["artificial_compressibility_method"].as<bool>();
    epsilon_acm                        = computational_parameters["epsilon_acm"].as<double>();
    final_time_iter                    = computational_parameters["final_time_iter"].as<int>();

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
    } else if( dummy_type_boco == "SUPERSONIC_INFLOW" ) {
        bocos_type[_WEST_] = _SUPERSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUPERSONIC_OUTFLOW" ) {
        bocos_type[_WEST_] = _SUPERSONIC_OUTFLOW_;
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
    } else if( dummy_type_boco == "SUPERSONIC_INFLOW" ) {
        bocos_type[_EAST_] = _SUPERSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUPERSONIC_OUTFLOW" ) {
        bocos_type[_EAST_] = _SUPERSONIC_OUTFLOW_;
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
    } else if( dummy_type_boco == "SUPERSONIC_INFLOW" ) {
        bocos_type[_SOUTH_] = _SUPERSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUPERSONIC_OUTFLOW" ) {
        bocos_type[_SOUTH_] = _SUPERSONIC_OUTFLOW_;
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
    } else if( dummy_type_boco == "SUPERSONIC_INFLOW" ) {
        bocos_type[_NORTH_] = _SUPERSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUPERSONIC_OUTFLOW" ) {
        bocos_type[_NORTH_] = _SUPERSONIC_OUTFLOW_;
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
    } else if( dummy_type_boco == "SUPERSONIC_INFLOW" ) {
        bocos_type[_BACK_] = _SUPERSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUPERSONIC_OUTFLOW" ) {
        bocos_type[_BACK_] = _SUPERSONIC_OUTFLOW_;
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
    } else if( dummy_type_boco == "SUPERSONIC_INFLOW" ) {
        bocos_type[_FRONT_] = _SUPERSONIC_INFLOW_;
    } else if( dummy_type_boco == "SUPERSONIC_OUTFLOW" ) {
        bocos_type[_FRONT_] = _SUPERSONIC_OUTFLOW_;
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

    /// Temporal point probes
    const YAML::Node & temporal_point_probes = configuration["temporal_point_probes"];
    number_temporal_point_probes = temporal_point_probes["number_temporal_point_probes"].as<int>();
    tpp_x_positions.resize( number_temporal_point_probes );
    tpp_y_positions.resize( number_temporal_point_probes );
    tpp_z_positions.resize( number_temporal_point_probes );
    tpp_output_frequency_iters.resize( number_temporal_point_probes );
    tpp_output_file_names.resize( number_temporal_point_probes );
    string yaml_input_name;
    for(int tpp = 0; tpp < number_temporal_point_probes; ++tpp) {
	yaml_input_name = "probe_" + to_string( tpp + 1 ) + "_x_position";
        tpp_x_positions[tpp] = temporal_point_probes[yaml_input_name].as<double>();
	yaml_input_name = "probe_" + to_string( tpp + 1 ) + "_y_position";
        tpp_y_positions[tpp] = temporal_point_probes[yaml_input_name].as<double>();
	yaml_input_name = "probe_" + to_string( tpp + 1 ) + "_z_position";
        tpp_z_positions[tpp] = temporal_point_probes[yaml_input_name].as<double>();
	yaml_input_name = "probe_" + to_string( tpp + 1 ) + "_output_frequency_iter";
        tpp_output_frequency_iters[tpp] = temporal_point_probes[yaml_input_name].as<int>();
	yaml_input_name = "probe_" + to_string( tpp + 1 ) + "_output_data_file_name";
        tpp_output_file_names[tpp] = temporal_point_probes[yaml_input_name].as<string>();
    }	    

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

void FlowSolverRHEA::fillMeshCoordinatesSizesFields() {

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

    /// Update halo values (do not activate!)
    //x_field.update();
    //y_field.update();
    //z_field.update();

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

	/// Reset avg and fluctuating fields
        avg_rho_field       = 0.0;
        avg_rhou_field      = 0.0;
        avg_rhov_field      = 0.0;
        avg_rhow_field      = 0.0;
        avg_rhoE_field      = 0.0;
        avg_rhoP_field      = 0.0;
        avg_rhoT_field      = 0.0;
        avg_u_field         = 0.0;
        avg_v_field         = 0.0;
        avg_w_field         = 0.0;
        avg_E_field         = 0.0;
        avg_P_field         = 0.0;
        avg_T_field         = 0.0;
        avg_sos_field       = 0.0;
        avg_mu_field        = 0.0;
        avg_kappa_field     = 0.0;
        avg_c_v_field       = 0.0;
        avg_c_p_field       = 0.0;
        rmsf_rho_field      = 0.0;
        rmsf_rhou_field     = 0.0;
        rmsf_rhov_field     = 0.0;
        rmsf_rhow_field     = 0.0;
        rmsf_rhoE_field     = 0.0;
        rmsf_u_field        = 0.0;
        rmsf_v_field        = 0.0;
        rmsf_w_field        = 0.0;
        rmsf_E_field        = 0.0;
        rmsf_P_field        = 0.0;
        rmsf_T_field        = 0.0;
        rmsf_sos_field      = 0.0;
        rmsf_mu_field       = 0.0;
        rmsf_kappa_field    = 0.0;
        rmsf_c_v_field      = 0.0;
        rmsf_c_p_field      = 0.0;
        favre_uffuff_field = 0.0;
        favre_uffvff_field = 0.0;
        favre_uffwff_field = 0.0;
        favre_vffvff_field = 0.0;
        favre_vffwff_field = 0.0;
        favre_wffwff_field = 0.0;
        favre_uffEff_field = 0.0;
        favre_vffEff_field = 0.0;
        favre_wffEff_field = 0.0;

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
    c_v_field.update();
    c_p_field.update();
    avg_rho_field.update();
    avg_rhou_field.update();
    avg_rhov_field.update();
    avg_rhow_field.update();
    avg_rhoE_field.update();
    avg_rhoP_field.update();
    avg_rhoT_field.update();
    avg_u_field.update();
    avg_v_field.update();
    avg_w_field.update();
    avg_E_field.update();
    avg_P_field.update();
    avg_T_field.update();
    avg_sos_field.update();
    avg_mu_field.update();
    avg_kappa_field.update();
    avg_c_v_field.update();
    avg_c_p_field.update();
    rmsf_rho_field.update();
    rmsf_rhou_field.update();
    rmsf_rhov_field.update();
    rmsf_rhow_field.update();
    rmsf_rhoE_field.update();
    rmsf_u_field.update();
    rmsf_v_field.update();
    rmsf_w_field.update();
    rmsf_E_field.update();
    rmsf_P_field.update();
    rmsf_T_field.update();
    rmsf_sos_field.update();
    rmsf_mu_field.update();
    rmsf_kappa_field.update();
    rmsf_c_v_field.update();
    rmsf_c_p_field.update();
    favre_uffuff_field.update();
    favre_uffvff_field.update();
    favre_uffwff_field.update();
    favre_vffvff_field.update();
    favre_vffwff_field.update();
    favre_wffwff_field.update();
    favre_uffEff_field.update();
    favre_vffEff_field.update();
    favre_wffEff_field.update();

    /// Fill mesh x, y, z, delta_x, delta_y, delta_z fields
    this->fillMeshCoordinatesSizesFields();

};

void FlowSolverRHEA::initializeThermodynamics() {

    /// All (inner, halo, boundary): rho, E, sos, c_v and c_p
    double rho, e, ke, c_v, c_p;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                thermodynamics->calculateDensityInternalEnergyFromPressureTemperature( rho, e, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)] );
                rho_field[I1D(i,j,k)] = rho;
                ke                    = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) );
                E_field[I1D(i,j,k)]   = e + ke;
                sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;
            }
        }
    }

    /// Update halo values
    rho_field.update();
    E_field.update();
    sos_field.update();
    c_v_field.update();
    c_p_field.update();

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
    //rhou_field.update();
    //rhov_field.update();
    //rhow_field.update();
    //rhoE_field.update();

};

void FlowSolverRHEA::conservedToPrimitiveVariables() {

    /// All (inner, halo, boundary) points: u, v, w and E
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    const int local_size_x = _lNx_;
    const int local_size_y = _lNy_;
    const int local_size_z = _lNz_;
    const int local_size   = local_size_x*local_size_y*local_size_z;
    const int inix = topo->iter_common[_ALL_][_INIX_];
    const int iniy = topo->iter_common[_ALL_][_INIY_];
    const int iniz = topo->iter_common[_ALL_][_INIZ_];
    const int endx = topo->iter_common[_ALL_][_ENDX_];      
    const int endy = topo->iter_common[_ALL_][_ENDY_];
    const int endz = topo->iter_common[_ALL_][_ENDZ_];
    #pragma acc enter data copyin(this)
    #pragma acc parallel loop collapse (3) copyin (rho_field.vector[0:local_size],rhou_field.vector[0:local_size],rhov_field.vector[0:local_size],rhow_field.vector[0:local_size],rhoE_field.vector[0:local_size]) copyout (u_field.vector[0:local_size],v_field.vector[0:local_size],w_field.vector[0:local_size],E_field.vector[0:local_size])
    for(int i = inix; i <= endx; i++) {
        for(int j = iniy; j <= endy; j++) {
            for(int k = iniz; k <= endz; k++) {
#else
    #pragma acc parallel loop collapse (3) async
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
#endif
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
            
    if( transport_pressure_scheme ) {

        /// All (inner, halo, boundary) points: T, E, rhoE, sos, c_v and c_p
        double T, e, ke, c_v, c_p;
        for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
	            T = T_field[I1D(i,j,k)]; 
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T, P_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                    T_field[I1D(i,j,k)]    = T;
                    e  = thermodynamics->calculateInternalEnergyFromPressureTemperatureDensity( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] ); 
                    ke = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) ); 
                    E_field[I1D(i,j,k)]    = e + ke;
                    rhoE_field[I1D(i,j,k)] = rho_field[I1D(i,j,k)]*E_field[I1D(i,j,k)];
                    sos_field[I1D(i,j,k)]  = thermodynamics->calculateSoundSpeed( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                    c_v_field[I1D(i,j,k)]  = c_v;
                    c_p_field[I1D(i,j,k)]  = c_p;
                }
            }
        }

        /// Update halo values
        //T_field.update();
        //E_field.update();
        //rhoE_field.update();
        //sos_field.update();
        //c_v_field.update();
        //c_p_field.update();

    } else {

        /// All (inner, halo, boundary) points: P, T, sos, c_v and c_p
        double ke, e, P, T, c_v, c_p;
        for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
            for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
                for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                    ke = 0.5*( pow( u_field[I1D(i,j,k)], 2.0 ) + pow( v_field[I1D(i,j,k)], 2.0 ) + pow( w_field[I1D(i,j,k)], 2.0 ) ); 
                    e  = E_field[I1D(i,j,k)] - ke;
                    P = P_field[I1D(i,j,k)];	/// Initial pressure guess
                    T = T_field[I1D(i,j,k)]; 	/// Initial temperature guess
                    thermodynamics->calculatePressureTemperatureFromDensityInternalEnergy( P, T, rho_field[I1D(i,j,k)], e );
                    P_field[I1D(i,j,k)]   = P; 
                    T_field[I1D(i,j,k)]   = T; 
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                    c_v_field[I1D(i,j,k)] = c_v;
                    c_p_field[I1D(i,j,k)] = c_p;
                }
            }
        }

        /// Update halo values
        //P_field.update();
        //T_field.update();
        //sos_field.update();
        //c_v_field.update();
        //c_p_field.update();

    }

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
                u_in = u_field[I1D(i+1,j,k)];
                v_in = v_field[I1D(i+1,j,k)];
                w_in = w_field[I1D(i+1,j,k)];
                P_in = P_field[I1D(i+1,j,k)];
                T_in = T_field[I1D(i+1,j,k)];	
		/// Calculate ghost primitive variables
                u_g = ( bocos_u[_WEST_] - wg_in*u_in )/wg_g;
                v_g = ( bocos_v[_WEST_] - wg_in*v_in )/wg_g;
                w_g = ( bocos_w[_WEST_] - wg_in*w_in )/wg_g;
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
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*0.5*( L_1 + L_5 );
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
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_4;
                    T_g   = T_field[I1D(i,j,k)];
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T_g, P_g, rho_g );
		} else if( bocos_type[_WEST_] == _SUPERSONIC_INFLOW_ ) {
                    u_g = bocos_u[_WEST_];
                    v_g = bocos_v[_WEST_];
                    w_g = bocos_w[_WEST_];
                    P_g = bocos_P[_WEST_];
                    T_g = bocos_T[_WEST_];
		} else if( bocos_type[_WEST_] == _SUPERSONIC_OUTFLOW_ ) {
                    u_g = u_in;
                    v_g = v_in;
                    w_g = w_in;
                    P_g = P_in;
                    T_g = T_in;
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
		/// Update sos, c_v and c_p
		double c_v, c_p;
		if( artificial_compressibility_method ) {
                    sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_g, rho_g );
		} else {
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_g, T_g, rho_g );
		}
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;
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
                u_in = u_field[I1D(i-1,j,k)];
                v_in = v_field[I1D(i-1,j,k)];
                w_in = w_field[I1D(i-1,j,k)];
                P_in = P_field[I1D(i-1,j,k)];
                T_in = T_field[I1D(i-1,j,k)];	
		/// Calculate ghost primitive variables
                u_g = ( bocos_u[_EAST_] - wg_in*u_in )/wg_g;
                v_g = ( bocos_v[_EAST_] - wg_in*v_in )/wg_g;
                w_g = ( bocos_w[_EAST_] - wg_in*w_in )/wg_g;
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
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*0.5*( L_1 + L_5 );
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
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_4;
                    T_g   = T_field[I1D(i,j,k)];
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T_g, P_g, rho_g );
		} else if( bocos_type[_EAST_] == _SUPERSONIC_INFLOW_ ) {
                    u_g = bocos_u[_EAST_];
                    v_g = bocos_v[_EAST_];
                    w_g = bocos_w[_EAST_];
                    P_g = bocos_P[_EAST_];
                    T_g = bocos_T[_EAST_];
		} else if( bocos_type[_EAST_] == _SUPERSONIC_OUTFLOW_ ) {
                    u_g = u_in;
                    v_g = v_in;
                    w_g = w_in;
                    P_g = P_in;
                    T_g = T_in;
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
		/// Update sos, c_v and c_p
		double c_v, c_p;
		if( artificial_compressibility_method ) {
                    sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_g, rho_g );
		} else {
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_g, T_g, rho_g );
		}
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;		
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
                u_in = u_field[I1D(i,j+1,k)];
                v_in = v_field[I1D(i,j+1,k)];
                w_in = w_field[I1D(i,j+1,k)];
                P_in = P_field[I1D(i,j+1,k)];
                T_in = T_field[I1D(i,j+1,k)];	
		/// Calculate ghost primitive variables
                u_g = ( bocos_u[_SOUTH_] - wg_in*u_in )/wg_g;
                v_g = ( bocos_v[_SOUTH_] - wg_in*v_in )/wg_g;
                w_g = ( bocos_w[_SOUTH_] - wg_in*w_in )/wg_g;
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
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*0.5*( L_1 + L_5 );
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
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_4;
                    T_g   = T_field[I1D(i,j,k)];
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T_g, P_g, rho_g );
		} else if( bocos_type[_SOUTH_] == _SUPERSONIC_INFLOW_ ) {
                    u_g = bocos_u[_SOUTH_];
                    v_g = bocos_v[_SOUTH_];
                    w_g = bocos_w[_SOUTH_];
                    P_g = bocos_P[_SOUTH_];
                    T_g = bocos_T[_SOUTH_];
		} else if( bocos_type[_SOUTH_] == _SUPERSONIC_OUTFLOW_ ) {
                    u_g = u_in;
                    v_g = v_in;
                    w_g = w_in;
                    P_g = P_in;
                    T_g = T_in;
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
		/// Update sos, c_v and c_p
		double c_v, c_p;
		if( artificial_compressibility_method ) {
                    sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_g, rho_g );
		} else {
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_g, T_g, rho_g );
		}
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;
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
                u_in = u_field[I1D(i,j-1,k)];
                v_in = v_field[I1D(i,j-1,k)];
                w_in = w_field[I1D(i,j-1,k)];
                P_in = P_field[I1D(i,j-1,k)];
                T_in = T_field[I1D(i,j-1,k)];	
		/// Calculate ghost primitive variables
                u_g = ( bocos_u[_NORTH_] - wg_in*u_in )/wg_g;
                v_g = ( bocos_v[_NORTH_] - wg_in*v_in )/wg_g;
                w_g = ( bocos_w[_NORTH_] - wg_in*w_in )/wg_g;
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
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*0.5*( L_1 + L_5 );
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
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_4;
                    T_g   = T_field[I1D(i,j,k)];
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T_g, P_g, rho_g );
		} else if( bocos_type[_NORTH_] == _SUPERSONIC_INFLOW_ ) {
                    u_g = bocos_u[_NORTH_];
                    v_g = bocos_v[_NORTH_];
                    w_g = bocos_w[_NORTH_];
                    P_g = bocos_P[_NORTH_];
                    T_g = bocos_T[_NORTH_];
		} else if( bocos_type[_NORTH_] == _SUPERSONIC_OUTFLOW_ ) {
                    u_g = u_in;
                    v_g = v_in;
                    w_g = w_in;
                    P_g = P_in;
                    T_g = T_in;
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
		/// Update sos, c_v and c_p
		double c_v, c_p;
		if( artificial_compressibility_method ) {
                    sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_g, rho_g );
		} else {
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_g, T_g, rho_g );
		}
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;
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
                u_in = u_field[I1D(i,j,k+1)];
                v_in = v_field[I1D(i,j,k+1)];
                w_in = w_field[I1D(i,j,k+1)];
                P_in = P_field[I1D(i,j,k+1)];
                T_in = T_field[I1D(i,j,k+1)];	
		/// Calculate ghost primitive variables
                u_g = ( bocos_u[_BACK_] - wg_in*u_in )/wg_g;
                v_g = ( bocos_v[_BACK_] - wg_in*v_in )/wg_g;
                w_g = ( bocos_w[_BACK_] - wg_in*w_in )/wg_g;
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
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*0.5*( L_1 + L_5 );
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
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_4;
                    T_g   = T_field[I1D(i,j,k)];
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T_g, P_g, rho_g );
		} else if( bocos_type[_BACK_] == _SUPERSONIC_INFLOW_ ) {
                    u_g = bocos_u[_BACK_];
                    v_g = bocos_v[_BACK_];
                    w_g = bocos_w[_BACK_];
                    P_g = bocos_P[_BACK_];
                    T_g = bocos_T[_BACK_];
		} else if( bocos_type[_BACK_] == _SUPERSONIC_OUTFLOW_ ) {
                    u_g = u_in;
                    v_g = v_in;
                    w_g = w_in;
                    P_g = P_in;
                    T_g = T_in;
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
		/// Update sos, c_v and c_p
		double c_v, c_p;
		if( artificial_compressibility_method ) {
                    sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_g, rho_g );
		} else {
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_g, T_g, rho_g );
		}
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;
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
                u_in = u_field[I1D(i,j,k-1)];
                v_in = v_field[I1D(i,j,k-1)];
                w_in = w_field[I1D(i,j,k-1)];
                P_in = P_field[I1D(i,j,k-1)];
                T_in = T_field[I1D(i,j,k-1)];	
		/// Calculate ghost primitive variables
                u_g = ( bocos_u[_FRONT_] - wg_in*u_in )/wg_g;
                v_g = ( bocos_v[_FRONT_] - wg_in*v_in )/wg_g;
                w_g = ( bocos_w[_FRONT_] - wg_in*w_in )/wg_g;
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
                    P_g = P_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*0.5*( L_1 + L_5 );
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
                    rho_g = rho_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( sos_in*sos_in ) )*( L_2 + 0.5*( L_1 + L_5 ) );
                    u_g   = u_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*( 1.0/( 2.0*rho_in*sos_in ) )*( L_5 - L_1 );
                    v_g   = v_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_3;
                    w_g   = w_field[I1D(i,j,k)] - ( delta_t/rk_number_stages )*L_4;
                    T_g   = T_field[I1D(i,j,k)];
                    thermodynamics->calculateTemperatureFromPressureDensityWithInitialGuess( T_g, P_g, rho_g );
		} else if( bocos_type[_FRONT_] == _SUPERSONIC_INFLOW_ ) {
                    u_g = bocos_u[_FRONT_];
                    v_g = bocos_v[_FRONT_];
                    w_g = bocos_w[_FRONT_];
                    P_g = bocos_P[_FRONT_];
                    T_g = bocos_T[_FRONT_];
		} else if( bocos_type[_FRONT_] == _SUPERSONIC_OUTFLOW_ ) {
                    u_g = u_in;
                    v_g = v_in;
                    w_g = w_in;
                    P_g = P_in;
                    T_g = T_in;
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
		/// Update sos, c_v and c_p
		double c_v, c_p;
		if( artificial_compressibility_method ) {
                    sos_field[I1D(i,j,k)] = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_g, rho_g );
		} else {
                    sos_field[I1D(i,j,k)] = thermodynamics->calculateSoundSpeed( P_g, T_g, rho_g );
                    thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_g, T_g, rho_g );
		}
                c_v_field[I1D(i,j,k)] = c_v;
                c_p_field[I1D(i,j,k)] = c_p;
            }
        }
    }

    /// Update halo values
    //rho_field.update();
    //rhou_field.update();
    //rhov_field.update();
    //rhow_field.update();
    //rhoE_field.update();
    //u_field.update();
    u_field.fillEdgeCornerBoundaries();
    //v_field.update();
    v_field.fillEdgeCornerBoundaries();
    //w_field.update();
    w_field.fillEdgeCornerBoundaries();
    //E_field.update();
    //P_field.update();
    //T_field.update();
    //sos_field.update();
    //c_v_field.update();
    //c_p_field.update();

};

void FlowSolverRHEA::updatePreviousStateConservedVariables() {

    /// All (inner, halo, boundary) points: rho_0, rhou_0 rhov_0, rhow_0, rhoE_0 and P_0
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    #pragma acc kernels loop collapse(3) independent
#endif
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                rho_0_field[I1D(i,j,k)]  = rho_field[I1D(i,j,k)]; 
                rhou_0_field[I1D(i,j,k)] = rhou_field[I1D(i,j,k)]; 
                rhov_0_field[I1D(i,j,k)] = rhov_field[I1D(i,j,k)]; 
                rhow_0_field[I1D(i,j,k)] = rhow_field[I1D(i,j,k)]; 
                rhoE_0_field[I1D(i,j,k)] = rhoE_field[I1D(i,j,k)]; 
                P_0_field[I1D(i,j,k)]    = P_field[I1D(i,j,k)]; 
            }
        }
    }

    /// Update halo values
    //rho_0_field.update();
    //rhou_0_field.update();
    //rhov_0_field.update();
    //rhow_0_field.update();
    //rhoE_0_field.update();
    //P_0_field.update();

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
    double sos, c_p;
    double delta_x, delta_y, delta_z;
    double S_x, S_y, S_z;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Speed of sound
		sos = sos_field[I1D(i,j,k)];
                /// Heat capacities
                //thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_field[I1D(i,j,k)], T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
		c_p = c_p_field[I1D(i,j,k)];
                /// Geometric stuff
                delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                /// x-direction inviscid, viscous & thermal terms
                S_x           = abs( u_field[I1D(i,j,k)] ) + sos;
                local_delta_t = min( local_delta_t, CFL*delta_x/S_x );											/// acoustic scale
                local_delta_t = min( local_delta_t, CFL*rho_field[I1D(i,j,k)]*pow( delta_x, 2.0 )/( mu_field[I1D(i,j,k)] + epsilon ) );			/// viscous scale
                local_delta_t = min( local_delta_t, CFL*rho_field[I1D(i,j,k)]*c_p*pow( delta_x, 2.0 )/( kappa_field[I1D(i,j,k)] + epsilon ) );		/// thermal diffusivity scale
                /// y-direction inviscid, viscous & thermal terms
                S_y           = abs( v_field[I1D(i,j,k)] ) + sos;
                local_delta_t = min( local_delta_t, CFL*delta_y/S_y );											/// acoustic scale
                local_delta_t = min( local_delta_t, CFL*rho_field[I1D(i,j,k)]*pow( delta_y, 2.0 )/( mu_field[I1D(i,j,k)] + epsilon ) );			/// viscous scale
                local_delta_t = min( local_delta_t, CFL*rho_field[I1D(i,j,k)]*c_p*pow( delta_y, 2.0 )/( kappa_field[I1D(i,j,k)] + epsilon ) );		/// thermal diffusivity scale
                /// z-direction inviscid, viscous & thermal terms
                S_z           = abs( w_field[I1D(i,j,k)] ) + sos;
                local_delta_t = min( local_delta_t, CFL*delta_z/S_z );											/// acoustic scale
                local_delta_t = min( local_delta_t, CFL*rho_field[I1D(i,j,k)]*pow( delta_z, 2.0 )/( mu_field[I1D(i,j,k)] + epsilon ) );			/// viscous scale
                local_delta_t = min( local_delta_t, CFL*rho_field[I1D(i,j,k)]*c_p*pow( delta_z, 2.0 )/( kappa_field[I1D(i,j,k)] + epsilon ) );		/// thermal diffusivity scale
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
    //mu_field.update();
    //kappa_field.update();

};

void FlowSolverRHEA::calculateSourceTerms() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    #pragma acc kernels loop collapse(3) independent
#endif
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

    /// Unsplit method for Euler equations:
    /// E. F. Toro.
    /// Riemann solvers and numerical methods for fluid dynamics.
    /// Springer, 2009.

    /// Inner points: rho, rhou, rhov, rhow and rhoE
    int index_L, index_R, var_type;
    double delta_x, delta_y, delta_z;
    double rho_L, u_L, v_L, w_L, E_L, P_L, P_rhouvw_L, a_L;
    double rho_R, u_R, v_R, w_R, E_R, P_R, P_rhouvw_R, a_R;
    double rho_F_p, rho_F_m, rhou_F_p, rhou_F_m, rhov_F_p;
    double rhov_F_m, rhow_F_p, rhow_F_m, rhoE_F_p, rhoE_F_m;
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    const int local_size_x = _lNx_;
    const int local_size_y = _lNy_;
    const int local_size_z = _lNz_;
    const int local_size   = local_size_x*local_size_y*local_size_z;
    const int inix = topo->iter_common[_INNER_][_INIX_];
    const int iniy = topo->iter_common[_INNER_][_INIY_];
    const int iniz = topo->iter_common[_INNER_][_INIZ_];
    const int endx = topo->iter_common[_INNER_][_ENDX_];
    const int endy = topo->iter_common[_INNER_][_ENDY_];
    const int endz = topo->iter_common[_INNER_][_ENDZ_];
    #pragma acc enter data copyin(this) 
    #pragma acc parallel loop collapse (3) copyin (rho_field.vector[0:local_size],u_field.vector[0:local_size],v_field.vector[0:local_size],w_field.vector[0:local_size],E_field.vector[0:local_size],P_field.vector[0:local_size],sos_field.vector[0:local_size],x_field.vector[0:local_size],y_field.vector[0:local_size],z_field.vector[0:local_size]) copyout (rho_inv_flux.vector[0:local_size],rhou_inv_flux.vector[0:local_size],rhov_inv_flux.vector[0:local_size],rhow_inv_flux.vector[0:local_size],rhoE_inv_flux.vector[0:local_size]) 
    for(int i = inix; i <= endx; i++) {
        for(int j = iniy; j <= endy; j++) {
            for(int k = iniz; k <= endz; k++) {
#else
    #pragma acc parallel loop collapse (3)
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
#endif
                /// Geometric stuff
                delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                /// x-direction i+1/2
                index_L = i;                           index_R = i + 1;
                rho_L   = rho_field[I1D(index_L,j,k)]; rho_R   = rho_field[I1D(index_R,j,k)]; 
                u_L     = u_field[I1D(index_L,j,k)];   u_R     = u_field[I1D(index_R,j,k)];
                v_L     = v_field[I1D(index_L,j,k)];   v_R     = v_field[I1D(index_R,j,k)];
                w_L     = w_field[I1D(index_L,j,k)];   w_R     = w_field[I1D(index_R,j,k)];
                E_L     = E_field[I1D(index_L,j,k)];   E_R     = E_field[I1D(index_R,j,k)];
                P_L     = P_field[I1D(index_L,j,k)];   P_R     = P_field[I1D(index_R,j,k)];
                a_L     = sos_field[I1D(index_L,j,k)]; a_R     = sos_field[I1D(index_R,j,k)];
                P_rhouvw_L = P_L - P_thermo;           P_rhouvw_R = P_R - P_thermo;	/// P_Thermo = 0.0 when ACM is deactivated
                /// rho
                var_type = 0;
                //rho_F_p  = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rho_F_p  = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 1;
                //rhou_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhou_F_p = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                //rhov_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhov_F_p = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                //rhow_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhow_F_p = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                //rhoE_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rhoE_F_p = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// x-direction i-1/2
                index_L = i - 1;                       index_R = i;
                rho_L   = rho_field[I1D(index_L,j,k)]; rho_R   = rho_field[I1D(index_R,j,k)];
                u_L     = u_field[I1D(index_L,j,k)];   u_R     = u_field[I1D(index_R,j,k)];
                v_L     = v_field[I1D(index_L,j,k)];   v_R     = v_field[I1D(index_R,j,k)];
                w_L     = w_field[I1D(index_L,j,k)];   w_R     = w_field[I1D(index_R,j,k)];
                E_L     = E_field[I1D(index_L,j,k)];   E_R     = E_field[I1D(index_R,j,k)];
                P_L     = P_field[I1D(index_L,j,k)];   P_R     = P_field[I1D(index_R,j,k)];
                a_L     = sos_field[I1D(index_L,j,k)]; a_R     = sos_field[I1D(index_R,j,k)];
                P_rhouvw_L = P_L - P_thermo;           P_rhouvw_R = P_R - P_thermo;	/// P_Thermo = 0.0 when ACM is deactivated
                /// rho
                var_type = 0;
                //rho_F_m  = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rho_F_m  = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 1;
                //rhou_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhou_F_m = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                //rhov_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhov_F_m = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                //rhow_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhow_F_m = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                //rhoE_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rhoE_F_m = calculateIntercellFlux( rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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
                P_rhouvw_L = P_L - P_thermo;           P_rhouvw_R = P_R - P_thermo;	/// P_Thermo = 0.0 when ACM is deactivated
                /// rho
                var_type = 0;
                //rho_F_p  = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rho_F_p  = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 2;
                //rhou_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhou_F_p = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhov
                var_type = 1;
                //rhov_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhov_F_p = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                //rhow_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhow_F_p = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                //rhoE_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rhoE_F_p = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// y-direction j-1/2
                index_L = j - 1;                       index_R = j;
                rho_L   = rho_field[I1D(i,index_L,k)]; rho_R   = rho_field[I1D(i,index_R,k)];
                u_L     = u_field[I1D(i,index_L,k)];   u_R     = u_field[I1D(i,index_R,k)];
                v_L     = v_field[I1D(i,index_L,k)];   v_R     = v_field[I1D(i,index_R,k)];
                w_L     = w_field[I1D(i,index_L,k)];   w_R     = w_field[I1D(i,index_R,k)];
                E_L     = E_field[I1D(i,index_L,k)];   E_R     = E_field[I1D(i,index_R,k)];
                P_L     = P_field[I1D(i,index_L,k)];   P_R     = P_field[I1D(i,index_R,k)];
                a_L     = sos_field[I1D(i,index_L,k)]; a_R     = sos_field[I1D(i,index_R,k)];
                P_rhouvw_L = P_L - P_thermo;           P_rhouvw_R = P_R - P_thermo;	/// P_Thermo = 0.0 when ACM is deactivated
                /// rho
                var_type = 0;
                //rho_F_m  = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rho_F_m  = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 2;
                //rhou_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhou_F_m = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhov
                var_type = 1;
                //rhov_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhov_F_m = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhow
                var_type = 3;
                //rhow_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhow_F_m = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                //rhoE_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rhoE_F_m = calculateIntercellFlux( rho_L, rho_R, v_L, v_R, u_L, u_R, w_L, w_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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
                P_rhouvw_L = P_L - P_thermo;           P_rhouvw_R = P_R - P_thermo;	/// P_Thermo = 0.0 when ACM is deactivated
                /// rho
                var_type = 0;
                //rho_F_p  = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rho_F_p  = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 3;
                //rhou_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhou_F_p = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                //rhov_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhov_F_p = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhow
                var_type = 1;
                //rhow_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhow_F_p = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                //rhoE_F_p = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rhoE_F_p = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// z-direction k-1/2
                index_L = k - 1;                       index_R = k;
                rho_L   = rho_field[I1D(i,j,index_L)]; rho_R   = rho_field[I1D(i,j,index_R)];
                u_L     = u_field[I1D(i,j,index_L)];   u_R     = u_field[I1D(i,j,index_R)];
                v_L     = v_field[I1D(i,j,index_L)];   v_R     = v_field[I1D(i,j,index_R)];
                w_L     = w_field[I1D(i,j,index_L)];   w_R     = w_field[I1D(i,j,index_R)];
                E_L     = E_field[I1D(i,j,index_L)];   E_R     = E_field[I1D(i,j,index_R)];
                P_L     = P_field[I1D(i,j,index_L)];   P_R     = P_field[I1D(i,j,index_R)];
                a_L     = sos_field[I1D(i,j,index_L)]; a_R     = sos_field[I1D(i,j,index_R)];
                P_rhouvw_L = P_L - P_thermo;           P_rhouvw_R = P_R - P_thermo;	/// P_Thermo = 0.0 when ACM is deactivated
                /// rho
                var_type = 0;
                //rho_F_m  = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rho_F_m  = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                /// rhou
                var_type = 3;
                //rhou_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhou_F_m = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhov
                var_type = 2;
                //rhov_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhov_F_m = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhow
                var_type = 1;
                //rhow_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                rhow_F_m = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_rhouvw_L, P_rhouvw_R, a_L, a_R, var_type );
                /// rhoE
                var_type = 4;
                //rhoE_F_m = riemann_solver->calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
                rhoE_F_m = calculateIntercellFlux( rho_L, rho_R, w_L, w_R, v_L, v_R, u_L, u_R, E_L, E_R, P_L, P_R, a_L, a_R, var_type );
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

    /// Inner points: rho_vis_flux, rhou_vis_flux, rhov_vis_flux, rhow_vis_flux, rhoE_vis_flux and work_vis_rhoe_flux 
    double delta_x, delta_y, delta_z;
    double d_u_x, d_u_y, d_u_z, d_v_x, d_v_y, d_v_z, d_w_x, d_w_y, d_w_z;
    double d_T_x, d_T_y, d_T_z;
    double d_mu_x, d_mu_y, d_mu_z, d_kappa_x, d_kappa_y, d_kappa_z;
    double div_uvw, tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;
    double div_tau_x, div_tau_y, div_tau_z;
    double div_q, div_uvw_tau_rhoe, div_uvw_tau_rhoke, div_uvw_tau_rhoE;
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    const int local_size_x = _lNx_;
    const int local_size_y = _lNy_;
    const int local_size_z = _lNz_;
    const int local_size   = local_size_x*local_size_y*local_size_z;
    const int inix = topo->iter_common[_INNER_][_INIX_];
    const int iniy = topo->iter_common[_INNER_][_INIY_];
    const int iniz = topo->iter_common[_INNER_][_INIZ_];
    const int endx = topo->iter_common[_INNER_][_ENDX_];
    const int endy = topo->iter_common[_INNER_][_ENDY_];
    const int endz = topo->iter_common[_INNER_][_ENDZ_];
    #pragma acc enter data copyin(this)
    #pragma acc parallel loop collapse (3) copyin (u_field.vector[0:local_size],v_field.vector[0:local_size],w_field.vector[0:local_size],T_field.vector[0:local_size],mu_field.vector[0:local_size],kappa_field.vector[0:local_size],x_field.vector[0:local_size],y_field.vector[0:local_size],z_field.vector[0:local_size]) copyout (rhou_vis_flux.vector[0:local_size],rhov_vis_flux.vector[0:local_size],rhow_vis_flux.vector[0:local_size],rhoE_vis_flux.vector[0:local_size],work_vis_rhoe_flux.vector[0:local_size])
    for(int i = inix; i <= endx; i++) {
        for(int j = iniy; j <= endy; j++) {
            for(int k = iniz; k <= endz; k++) {  
#else
    #pragma acc parallel loop collapse (3) 
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
#endif
                /// Geometric stuff
                delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                /// Velocity derivatives
                d_u_x = ( u_field[I1D(i+1,j,k)] - u_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_u_y = ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_u_z = ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                d_v_x = ( v_field[I1D(i+1,j,k)] - v_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_v_y = ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_v_z = ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                d_w_x = ( w_field[I1D(i+1,j,k)] - w_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_w_y = ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_w_z = ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                /// Temperature derivatives
                d_T_x = ( T_field[I1D(i+1,j,k)] - T_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_T_y = ( T_field[I1D(i,j+1,k)] - T_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_T_z = ( T_field[I1D(i,j,k+1)] - T_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                /// Transport coefficients derivatives
                d_mu_x    = ( mu_field[I1D(i+1,j,k)] - mu_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_mu_y    = ( mu_field[I1D(i,j+1,k)] - mu_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_mu_z    = ( mu_field[I1D(i,j,k+1)] - mu_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                d_kappa_x = ( kappa_field[I1D(i+1,j,k)] - kappa_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_kappa_y = ( kappa_field[I1D(i,j+1,k)] - kappa_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_kappa_z = ( kappa_field[I1D(i,j,k+1)] - kappa_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
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
                div_tau_x = mu_field[I1D(i,j,k)]*( ( 1.00/delta_x )*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( x_field[I1D(i+1,j,k)] - x_field[I1D(i,j,k)] )
                                                                    - ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( x_field[I1D(i,j,k)] - x_field[I1D(i-1,j,k)] ) )
                                                 + ( 1.00/delta_y )*( ( u_field[I1D(i,j+1,k)] - u_field[I1D(i,j,k)] )/( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j,k)] )
                                                                    - ( u_field[I1D(i,j,k)] - u_field[I1D(i,j-1,k)] )/( y_field[I1D(i,j,k)] - y_field[I1D(i,j-1,k)] ) )
                                                 + ( 1.00/delta_z )*( ( u_field[I1D(i,j,k+1)] - u_field[I1D(i,j,k)] )/( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k)] )
                                                                    - ( u_field[I1D(i,j,k)] - u_field[I1D(i,j,k-1)] )/( z_field[I1D(i,j,k)] - z_field[I1D(i,j,k-1)] ) ) )
              + ( 1.0/3.0 )*mu_field[I1D(i,j,k)]*( ( 1.00/delta_x )*( ( u_field[I1D(i+1,j,k)] - u_field[I1D(i,j,k)] )/( x_field[I1D(i+1,j,k)] - x_field[I1D(i,j,k)] )
                                                                    - ( u_field[I1D(i,j,k)] - u_field[I1D(i-1,j,k)] )/( x_field[I1D(i,j,k)] - x_field[I1D(i-1,j,k)] ) )
                                                 + ( 0.25/delta_x )*( ( v_field[I1D(i+1,j+1,k)] - v_field[I1D(i+1,j-1,k)] )/delta_y
                                                                    - ( v_field[I1D(i-1,j+1,k)] - v_field[I1D(i-1,j-1,k)] )/delta_y )
                                                 + ( 0.25/delta_x )*( ( w_field[I1D(i+1,j,k+1)] - w_field[I1D(i+1,j,k-1)] )/delta_z
                                                                    - ( w_field[I1D(i-1,j,k+1)] - w_field[I1D(i-1,j,k-1)] )/delta_z ) )
	                  + ( d_mu_x*tau_xx + d_mu_y*tau_xy + d_mu_z*tau_xz )/( mu_field[I1D(i,j,k)] + epsilon );
                div_tau_y = mu_field[I1D(i,j,k)]*( ( 1.00/delta_x )*( ( v_field[I1D(i+1,j,k)] - v_field[I1D(i,j,k)] )/( x_field[I1D(i+1,j,k)] - x_field[I1D(i,j,k)] )
                                                                    - ( v_field[I1D(i,j,k)] - v_field[I1D(i-1,j,k)] )/( x_field[I1D(i,j,k)] - x_field[I1D(i-1,j,k)] ) )
                                                 + ( 1.00/delta_y )*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j,k)] )
                                                                    - ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( y_field[I1D(i,j,k)] - y_field[I1D(i,j-1,k)] ) )
                                                 + ( 1.00/delta_z )*( ( v_field[I1D(i,j,k+1)] - v_field[I1D(i,j,k)] )/( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k)] )
                                                                    - ( v_field[I1D(i,j,k)] - v_field[I1D(i,j,k-1)] )/( z_field[I1D(i,j,k)] - z_field[I1D(i,j,k-1)] ) ) )
              + ( 1.0/3.0 )*mu_field[I1D(i,j,k)]*( ( 0.25/delta_y )*( ( u_field[I1D(i+1,j+1,k)] - u_field[I1D(i-1,j+1,k)] )/delta_x
                                                                    - ( u_field[I1D(i+1,j-1,k)] - u_field[I1D(i-1,j-1,k)] )/delta_x )
                                                 + ( 1.00/delta_y )*( ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j,k)] )/( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j,k)] )
                                                                    - ( v_field[I1D(i,j,k)] - v_field[I1D(i,j-1,k)] )/( y_field[I1D(i,j,k)] - y_field[I1D(i,j-1,k)] ) )
                                                 + ( 0.25/delta_y )*( ( w_field[I1D(i,j+1,k+1)] - w_field[I1D(i,j+1,k-1)] )/delta_z
                                                                    - ( w_field[I1D(i,j-1,k+1)] - w_field[I1D(i,j-1,k-1)] )/delta_z ) )
	                  + ( d_mu_x*tau_xy + d_mu_y*tau_yy + d_mu_z*tau_yz )/( mu_field[I1D(i,j,k)] + epsilon );
                div_tau_z = mu_field[I1D(i,j,k)]*( ( 1.00/delta_x )*( ( w_field[I1D(i+1,j,k)] - w_field[I1D(i,j,k)] )/( x_field[I1D(i+1,j,k)] - x_field[I1D(i,j,k)] )
                                                                    - ( w_field[I1D(i,j,k)] - w_field[I1D(i-1,j,k)] )/( x_field[I1D(i,j,k)] - x_field[I1D(i-1,j,k)] ) )
                                                 + ( 1.00/delta_y )*( ( w_field[I1D(i,j+1,k)] - w_field[I1D(i,j,k)] )/( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j,k)] )
                                                                    - ( w_field[I1D(i,j,k)] - w_field[I1D(i,j-1,k)] )/( y_field[I1D(i,j,k)] - y_field[I1D(i,j-1,k)] ) )
                                                 + ( 1.00/delta_z )*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k)] )
                                                                    - ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( z_field[I1D(i,j,k)] - z_field[I1D(i,j,k-1)] ) ) )
              + ( 1.0/3.0 )*mu_field[I1D(i,j,k)]*( ( 0.25/delta_z )*( ( u_field[I1D(i+1,j,k+1)] - u_field[I1D(i-1,j,k+1)] )/delta_x
                                                                    - ( u_field[I1D(i+1,j,k-1)] - u_field[I1D(i-1,j,k-1)] )/delta_x )
                                                 + ( 0.25/delta_z )*( ( v_field[I1D(i,j+1,k+1)] - v_field[I1D(i,j-1,k+1)] )/delta_y
                                                                    - ( v_field[I1D(i,j+1,k-1)] - v_field[I1D(i,j-1,k-1)] )/delta_y )
                                                 + ( 1.00/delta_z )*( ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k)] )/( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k)] )
                                                                    - ( w_field[I1D(i,j,k)] - w_field[I1D(i,j,k-1)] )/( z_field[I1D(i,j,k)] - z_field[I1D(i,j,k-1)] ) ) )
	                  + ( d_mu_x*tau_xz + d_mu_y*tau_yz + d_mu_z*tau_zz )/( mu_field[I1D(i,j,k)] + epsilon );
                /// Fourier term
                div_q = ( -1.0 )*kappa_field[I1D(i,j,k)]*( ( 1.0/delta_x )*( ( T_field[I1D(i+1,j,k)] - T_field[I1D(i,j,k)] )/( x_field[I1D(i+1,j,k)] - x_field[I1D(i,j,k)] )
                                                                           - ( T_field[I1D(i,j,k)] - T_field[I1D(i-1,j,k)] )/( x_field[I1D(i,j,k)] - x_field[I1D(i-1,j,k)] ) )
                                                         + ( 1.0/delta_y )*( ( T_field[I1D(i,j+1,k)] - T_field[I1D(i,j,k)] )/( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j,k)] )
                                                                           - ( T_field[I1D(i,j,k)] - T_field[I1D(i,j-1,k)] )/( y_field[I1D(i,j,k)] - y_field[I1D(i,j-1,k)] ) )
                                                         + ( 1.0/delta_z )*( ( T_field[I1D(i,j,k+1)] - T_field[I1D(i,j,k)] )/( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k)] )
                                                                           - ( T_field[I1D(i,j,k)] - T_field[I1D(i,j,k-1)] )/( z_field[I1D(i,j,k)] - z_field[I1D(i,j,k-1)] ) ) )
		      - d_kappa_x*d_T_x - d_kappa_y*d_T_y - d_kappa_z*d_T_z;
                /// Work of viscous stresses for internal energy
                div_uvw_tau_rhoe = tau_xx*d_u_x + tau_xy*d_u_y + tau_xz*d_u_z
                                 + tau_xy*d_v_x + tau_yy*d_v_y + tau_yz*d_v_z
                                 + tau_xz*d_w_x + tau_yz*d_w_y + tau_zz*d_w_z;
                /// Work of viscous stresses for kinetic energy
                div_uvw_tau_rhoke = u_field[I1D(i,j,k)]*div_tau_x + v_field[I1D(i,j,k)]*div_tau_y + w_field[I1D(i,j,k)]*div_tau_z;
                /// Work of viscous stresses for total energy
                div_uvw_tau_rhoE = div_uvw_tau_rhoe + div_uvw_tau_rhoke; 
                /// Viscous fluxes
                rhou_vis_flux[I1D(i,j,k)]      = div_tau_x;
                rhov_vis_flux[I1D(i,j,k)]      = div_tau_y;
                rhow_vis_flux[I1D(i,j,k)]      = div_tau_z;
                rhoE_vis_flux[I1D(i,j,k)]      = ( -1.0 )*div_q + div_uvw_tau_rhoE;
                work_vis_rhoe_flux[I1D(i,j,k)] = ( -1.0 )*div_q + div_uvw_tau_rhoe;
            }
        }
    }

    /// Update halo values
    //rhou_vis_flux.update();
    //rhov_vis_flux.update();
    //rhow_vis_flux.update();
    //rhoE_vis_flux.update();
    //work_vis_rhoe_flux.update();

};

void FlowSolverRHEA::timeAdvanceConservedVariables(const int &rk_time_stage) {

    /// Coefficients of explicit Runge-Kutta stages
    double rk_a = 0.0, rk_b = 0.0, rk_c = 0.0;
    runge_kutta_method->setStageCoefficients(rk_a,rk_b,rk_c,rk_time_stage);    

    /// Inner points: rho, rhou, rhov, rhow and rhoE
    double f_rhouvw = 0.0;
    double rho_rhs_flux = 0.0, rhou_rhs_flux = 0.0, rhov_rhs_flux = 0.0, rhow_rhs_flux = 0.0, rhoE_rhs_flux = 0.0;
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    const int local_size_x = _lNx_;
    const int local_size_y = _lNy_;
    const int local_size_z = _lNz_;
    const int local_size   = local_size_x*local_size_y*local_size_z;
    const int inix = topo->iter_common[_INNER_][_INIX_];
    const int iniy = topo->iter_common[_INNER_][_INIY_];
    const int iniz = topo->iter_common[_INNER_][_INIZ_];
    const int endx = topo->iter_common[_INNER_][_ENDX_];
    const int endy = topo->iter_common[_INNER_][_ENDY_];
    const int endz = topo->iter_common[_INNER_][_ENDZ_];
    #pragma acc enter data copyin (this)
    #pragma acc data copyin (u_field.vector[0:local_size],v_field.vector[0:local_size],w_field.vector[0:local_size])
    #pragma acc data copyin (rho_0_field.vector[0:local_size],rhou_0_field.vector[0:local_size],rhov_0_field.vector[0:local_size],rhow_0_field.vector[0:local_size],rhoE_0_field.vector[0:local_size])
    #pragma acc data copyin (rho_inv_flux.vector[0:local_size],rhou_inv_flux.vector[0:local_size],rhov_inv_flux.vector[0:local_size],rhow_inv_flux.vector[0:local_size],rhoE_inv_flux.vector[0:local_size])
    #pragma acc data copyin (rhou_vis_flux.vector[0:local_size],rhov_vis_flux.vector[0:local_size],rhow_vis_flux.vector[0:local_size],rhoE_vis_flux.vector[0:local_size])
    #pragma acc data copyin (f_rhou_field.vector[0:local_size],f_rhov_field.vector[0:local_size],f_rhow_field.vector[0:local_size],f_rhoE_field.vector[0:local_size])
    #pragma acc enter data copyin (rho_field.vector[0:local_size],rhou_field.vector[0:local_size],rhov_field.vector[0:local_size],rhow_field.vector[0:local_size],rhoE_field.vector[0:local_size])
    #pragma acc parallel loop collapse (3)
    for(int i = inix; i <= endx; i++) {
        for(int j = iniy; j <= endy; j++) {
            for(int k = iniz; k <= endz; k++) {
#else
    #pragma acc parallel loop collapse (3)
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
#endif
                /// Work of momentum sources
                f_rhouvw = f_rhou_field[I1D(i,j,k)]*u_field[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]*v_field[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]*w_field[I1D(i,j,k)];
                /// Sum right-hand-side (RHS) fluxes
                rho_rhs_flux  = ( -1.0 )*rho_inv_flux[I1D(i,j,k)]; 
                rhou_rhs_flux = ( -1.0 )*rhou_inv_flux[I1D(i,j,k)] + rhou_vis_flux[I1D(i,j,k)] + f_rhou_field[I1D(i,j,k)]; 
                rhov_rhs_flux = ( -1.0 )*rhov_inv_flux[I1D(i,j,k)] + rhov_vis_flux[I1D(i,j,k)] + f_rhov_field[I1D(i,j,k)]; 
                rhow_rhs_flux = ( -1.0 )*rhow_inv_flux[I1D(i,j,k)] + rhow_vis_flux[I1D(i,j,k)] + f_rhow_field[I1D(i,j,k)]; 
                rhoE_rhs_flux = ( -1.0 )*rhoE_inv_flux[I1D(i,j,k)] + rhoE_vis_flux[I1D(i,j,k)] + f_rhoE_field[I1D(i,j,k)] + f_rhouvw;
                /// Runge-Kutta step
                rho_field[I1D(i,j,k)]  = rk_a*rho_0_field[I1D(i,j,k)]  + rk_b*rho_field[I1D(i,j,k)]  + rk_c*delta_t*rho_rhs_flux;
                rhou_field[I1D(i,j,k)] = rk_a*rhou_0_field[I1D(i,j,k)] + rk_b*rhou_field[I1D(i,j,k)] + rk_c*delta_t*rhou_rhs_flux;
                rhov_field[I1D(i,j,k)] = rk_a*rhov_0_field[I1D(i,j,k)] + rk_b*rhov_field[I1D(i,j,k)] + rk_c*delta_t*rhov_rhs_flux;
                rhow_field[I1D(i,j,k)] = rk_a*rhow_0_field[I1D(i,j,k)] + rk_b*rhow_field[I1D(i,j,k)] + rk_c*delta_t*rhow_rhs_flux;
                rhoE_field[I1D(i,j,k)] = rk_a*rhoE_0_field[I1D(i,j,k)] + rk_b*rhoE_field[I1D(i,j,k)] + rk_c*delta_t*rhoE_rhs_flux;
	    }
        }
    }
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    #pragma acc exit data copyout (rho_field.vector[0:local_size],rhou_field.vector[0:local_size],rhov_field.vector[0:local_size],rhow_field.vector[0:local_size],rhoE_field.vector[0:local_size])
#endif    
    
    ///// Attention! Communications performed only at the last stage of the Runge-Kutta to improve computational performance
    ///// ... temporal integration is first-order at points connecting partitions
    //if( rk_time_stage == rk_number_stages ) {

        /// Update halo values
        rho_field.update();
        rhou_field.update();
        rhov_field.update();
        rhow_field.update();
        rhoE_field.update();

    //}

    if( transport_pressure_scheme ) {
        this->timeAdvancePressure(rk_time_stage);
    }

};

void FlowSolverRHEA::timeAdvancePressure(const int &rk_time_stage) {

    /// Coefficients of explicit Runge-Kutta stages
    double rk_a = 0.0, rk_b = 0.0, rk_c = 0.0;
    runge_kutta_method->setStageCoefficients(rk_a,rk_b,rk_c,rk_time_stage);    

    /// Inner points: P
    double P_inv_flux = 0.0, P_vis_flux = 0.0, P_rhs_flux = 0.0;
    double delta_x, delta_y, delta_z;
    double d_P_x, d_P_y, d_P_z, d_u_x, d_v_y, d_w_z;
    double div_uvw, bar_v;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                /// Pressure and velocity derivatives
                d_P_x = ( P_field[I1D(i+1,j,k)] - P_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_P_y = ( P_field[I1D(i,j+1,k)] - P_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_P_z = ( P_field[I1D(i,j,k+1)] - P_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                d_u_x = ( u_field[I1D(i+1,j,k)] - u_field[I1D(i-1,j,k)] )/( 2.0*delta_x );
                d_v_y = ( v_field[I1D(i,j+1,k)] - v_field[I1D(i,j-1,k)] )/( 2.0*delta_y );
                d_w_z = ( w_field[I1D(i,j,k+1)] - w_field[I1D(i,j,k-1)] )/( 2.0*delta_z );
                /// Divergence of velocity
                div_uvw = d_u_x + d_v_y + d_w_z;
                /// Inviscid flux
		P_inv_flux = u_field[I1D(i,j,k)]*d_P_x + v_field[I1D(i,j,k)]*d_P_y + w_field[I1D(i,j,k)]*d_P_z + rho_field[I1D(i,j,k)]*pow( sos_field[I1D(i,j,k)], 2.0 )*div_uvw;
                /// Viscous flux
		bar_v = thermodynamics->getMolecularWeight()/rho_field[I1D(i,j,k)];
		P_vis_flux = ( thermodynamics->calculateVolumeExpansivity( T_field[I1D(i,j,k)], bar_v )/( rho_field[I1D(i,j,k)]*c_v_field[I1D(i,j,k)]*thermodynamics->calculateIsothermalCompressibility( T_field[I1D(i,j,k)], bar_v ) ) )*work_vis_rhoe_flux[I1D(i,j,k)];
                /// Sum right-hand-side (RHS) fluxes
                P_rhs_flux = ( -1.0 )*P_inv_flux + P_vis_flux + f_rhoE_field[I1D(i,j,k)]; 
                /// Runge-Kutta step
                P_field[I1D(i,j,k)] = rk_a*P_0_field[I1D(i,j,k)] + rk_b*P_field[I1D(i,j,k)] + rk_c*delta_t*P_rhs_flux;
	    }
        }
    }
    
    ///// Attention! Communications performed only at the last stage of the Runge-Kutta to improve computational performance
    ///// ... temporal integration is first-order at points connecting partitions
    //if( rk_time_stage == rk_number_stages ) {

        /// Update halo values
        P_field.update();

    //}

};

void FlowSolverRHEA::outputCurrentStateData() {

    /// Write to file current solver state, time and time iteration
    writer_reader->setAttribute( "Time", current_time );
    writer_reader->setAttribute( "Iteration", current_time_iter );
    writer_reader->setAttribute( "AveragingTime", averaging_time );
    writer_reader->write( current_time_iter );

};

void FlowSolverRHEA::outputTemporalPointProbesData() {

    /// Initialize MPI stuff
    int my_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /// Iterate through temporal point probes
    for(int tpp = 0; tpp < number_temporal_point_probes; ++tpp) {
        /// Write temporal point probe data to file (if criterion satisfied)
        if( current_time_iter%tpp_output_frequency_iters[tpp] == 0 ) {
            /// Owner rank writes to file
            if( temporal_point_probes[tpp].getGlobalOwnerRank() == my_rank ) {
                int i_index, j_index, k_index;
                /// Get local indices i, j, k
		i_index = temporal_point_probes[tpp].getLocalIndexI(); 
		j_index = temporal_point_probes[tpp].getLocalIndexJ(); 
		k_index = temporal_point_probes[tpp].getLocalIndexK();
                /// Generate header string
                string output_header_string; 
	        output_header_string  = "# t [s], x [m], y[m], z[m], rho [kg/m3], u [m/s], v [m/s], w [m/s], E [J/kg], P [Pa], T [K], sos [m/s], mu [Pa·s], kappa [W/(m·K)], c_v [J/(kg·K)], c_p [J/(kg·K)]";
	        output_header_string += ", avg_rho [kg/m3], avg_rhou [kg/(s·m2)], avg_rhov [kg/(s·m2)], avg_rhow [kg/(s·m2)], avg_rhoE [J/m3], avg_rhoP [kg2/(m4·s2)], avg_rhoT [(kg·K)/m3]";
	        output_header_string += ", avg_u [m/s], avg_v [m/s], avg_w [m/s], avg_E [J/kg], avg_P [Pa], avg_T [K], avg_sos [m/s], avg_mu [Pa·s], avg_kappa [W/(m·K)], avg_c_v [J/(kg·K)], avg_c_p [J/(kg·K)]";
	        output_header_string += ", rmsf_rho [kg/m3], rmsf_rhou [kg/(s·m2)], rmsf_rhov [kg/(s·m2)], rmsf_rhow [kg/(s·m2)], rmsf_rhoE [J/m3]";
	        output_header_string += ", rmsf_u [m/s], rmsf_v [m/s], rmsf_w [m/s], rmsf_E [J/kg], rmsf_P [Pa], rmsf_T [K], rmsf_sos [m/s], rmsf_mu [Pa·s], rmsf_kappa [W/(m·K)], rmsf_c_v [J/(kg·K)], rmsf_c_p [J/(kg·K)]";
	        output_header_string += ", favre_uffuff [m2/s2], favre_uffvff [m2/s2], favre_uffwff [m2/s2], favre_vffvff [m2/s2], favre_vffwff [m2/s2], favre_wffwff [m2/s2]";
	        output_header_string += ", favre_uffEff [(m·J)/(s·kg)], favre_vffEff [(m·J)/(s·kg)], favre_wffEff [(m·J)/(s·kg)]";
                /// Generate data string
                string output_data_string; 
	        output_data_string    = to_string( current_time );
	        output_data_string   += "," + to_string( x_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( y_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( z_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( rho_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( u_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( v_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( w_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( E_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( P_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( T_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( sos_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( mu_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( kappa_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( c_v_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( c_p_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( avg_rho_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_rhou_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_rhov_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_rhow_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( avg_rhoE_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_rhoP_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_rhoT_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( avg_u_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_v_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_w_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( avg_E_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_P_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_T_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( avg_sos_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_mu_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_kappa_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( avg_c_v_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( avg_c_p_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( rmsf_rho_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_rhou_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_rhov_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_rhow_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_rhoE_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( rmsf_u_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_v_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_w_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( rmsf_E_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_P_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_T_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( rmsf_sos_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_mu_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_kappa_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( rmsf_c_v_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( rmsf_c_p_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( favre_uffuff_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( favre_uffvff_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( favre_uffwff_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( favre_vffvff_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( favre_vffwff_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( favre_wffwff_field[I1D(i_index,j_index,k_index)] );
	        output_data_string   += "," + to_string( favre_uffEff_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( favre_vffEff_field[I1D(i_index,j_index,k_index)] ) + "," + to_string( favre_wffEff_field[I1D(i_index,j_index,k_index)] );
                /// Write (header string) data string to file
                temporal_point_probes[tpp].writeDataStringToOutputFile(output_header_string, output_data_string);
	    }
	}
    }	    

};

void FlowSolverRHEA::updateTimeAveragedQuantities() {

    /// All (inner, boundary & halo) points: first- and second-order time statistics of flow quantities 
#if _OPENACC_MANUAL_DATA_MOVEMENT_
    const int local_size_x = _lNx_;
    const int local_size_y = _lNy_;
    const int local_size_z = _lNz_;
    const int local_size   = local_size_x*local_size_y*local_size_z;
    const int inix = topo->iter_common[_ALL_][_INIX_];
    const int iniy = topo->iter_common[_ALL_][_INIY_];
    const int iniz = topo->iter_common[_ALL_][_INIZ_];
    const int endx = topo->iter_common[_ALL_][_ENDX_];
    const int endy = topo->iter_common[_ALL_][_ENDY_];
    const int endz = topo->iter_common[_ALL_][_ENDZ_];
    #pragma acc enter data copyin (this)
    #pragma acc data copy (rho_field.vector[0:local_size],rhou_field.vector[0:local_size],rhov_field.vector[0:local_size],rhow_field.vector[0:local_size],rhoE_field.vector[0:local_size])
    #pragma acc data copy (u_field.vector[0:local_size],v_field.vector[0:local_size],w_field.vector[0:local_size],E_field.vector[0:local_size],P_field.vector[0:local_size],T_field.vector[0:local_size])
    #pragma acc data copy (sos_field.vector[0:local_size],mu_field.vector[0:local_size],kappa_field.vector[0:local_size],c_v_field.vector[0:local_size],c_p_field.vector[0:local_size])
    #pragma acc data copy (avg_rho_field.vector[0:local_size],avg_rhou_field.vector[0:local_size],avg_rhov_field.vector[0:local_size],avg_rhow_field.vector[0:local_size],avg_rhoE_field.vector[0:local_size],avg_rhoP_field.vector[0:local_size],avg_rhoT_field.vector[0:local_size])
    #pragma acc data copy (avg_u_field.vector[0:local_size],avg_v_field.vector[0:local_size],avg_w_field.vector[0:local_size],avg_E_field.vector[0:local_size],avg_P_field.vector[0:local_size],avg_T_field.vector[0:local_size])
    #pragma acc data copy (avg_sos_field.vector[0:local_size],avg_mu_field.vector[0:local_size],avg_kappa_field.vector[0:local_size],avg_c_v_field.vector[0:local_size],avg_c_p_field.vector[0:local_size])
    #pragma acc data copy (rmsf_rho_field.vector[0:local_size],rmsf_rhou_field.vector[0:local_size],rmsf_rhov_field.vector[0:local_size],rmsf_rhow_field.vector[0:local_size],rmsf_rhoE_field.vector[0:local_size])
    #pragma acc data copy (rmsf_u_field.vector[0:local_size],rmsf_v_field.vector[0:local_size],rmsf_w_field.vector[0:local_size],rmsf_E_field.vector[0:local_size],rmsf_P_field.vector[0:local_size],rmsf_T_field.vector[0:local_size])
    #pragma acc data copy (rmsf_sos_field.vector[0:local_size],rmsf_mu_field.vector[0:local_size],rmsf_kappa_field.vector[0:local_size],rmsf_c_v_field.vector[0:local_size],rmsf_c_p_field.vector[0:local_size])
    #pragma acc data copy (favre_uffuff_field.vector[0:local_size],favre_uffvff_field.vector[0:local_size],favre_uffwff_field.vector[0:local_size],favre_vffvff_field.vector[0:local_size],favre_vffwff_field.vector[0:local_size],favre_wffwff_field.vector[0:local_size])
    #pragma acc data copy (favre_uffEff_field.vector[0:local_size],favre_vffEff_field.vector[0:local_size],favre_wffEff_field.vector[0:local_size])
    #pragma acc parallel loop collapse (3)
    for(int i = inix; i <= endx; i++) {
        for(int j = iniy; j <= endy; j++) {
            for(int k = iniz; k <= endz; k++) {
#else
    #pragma acc parallel loop collapse (3) async
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
#endif
                /// Time-averaged quantities
                avg_rho_field[I1D(i,j,k)]   = updateTimeMeanQuantity(rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_rhou_field[I1D(i,j,k)]  = updateTimeMeanQuantity(rhou_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_rhov_field[I1D(i,j,k)]  = updateTimeMeanQuantity(rhov_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_rhow_field[I1D(i,j,k)]  = updateTimeMeanQuantity(rhow_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_rhoE_field[I1D(i,j,k)]  = updateTimeMeanQuantity(rhoE_field[I1D(i,j,k)],avg_rhoE_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_rhoP_field[I1D(i,j,k)]  = updateTimeMeanQuantity(rho_field[I1D(i,j,k)]*P_field[I1D(i,j,k)],avg_rhoP_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_rhoT_field[I1D(i,j,k)]  = updateTimeMeanQuantity(rho_field[I1D(i,j,k)]*T_field[I1D(i,j,k)],avg_rhoT_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_u_field[I1D(i,j,k)]     = updateTimeMeanQuantity(u_field[I1D(i,j,k)],avg_u_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_v_field[I1D(i,j,k)]     = updateTimeMeanQuantity(v_field[I1D(i,j,k)],avg_v_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_w_field[I1D(i,j,k)]     = updateTimeMeanQuantity(w_field[I1D(i,j,k)],avg_w_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_E_field[I1D(i,j,k)]     = updateTimeMeanQuantity(E_field[I1D(i,j,k)],avg_E_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_P_field[I1D(i,j,k)]     = updateTimeMeanQuantity(P_field[I1D(i,j,k)],avg_P_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_T_field[I1D(i,j,k)]     = updateTimeMeanQuantity(T_field[I1D(i,j,k)],avg_T_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_sos_field[I1D(i,j,k)]   = updateTimeMeanQuantity(sos_field[I1D(i,j,k)],avg_sos_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_mu_field[I1D(i,j,k)]    = updateTimeMeanQuantity(mu_field[I1D(i,j,k)],avg_mu_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_kappa_field[I1D(i,j,k)] = updateTimeMeanQuantity(kappa_field[I1D(i,j,k)],avg_kappa_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_c_v_field[I1D(i,j,k)]   = updateTimeMeanQuantity(c_v_field[I1D(i,j,k)],avg_c_v_field[I1D(i,j,k)],delta_t,averaging_time);
                avg_c_p_field[I1D(i,j,k)]   = updateTimeMeanQuantity(c_p_field[I1D(i,j,k)],avg_c_p_field[I1D(i,j,k)],delta_t,averaging_time);

                /// Root-mean-square-fluctuation quantities
                rmsf_rho_field[I1D(i,j,k)]   = updateTimeRmsfQuantity(rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],rmsf_rho_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_rhou_field[I1D(i,j,k)]  = updateTimeRmsfQuantity(rhou_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],rmsf_rhou_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_rhov_field[I1D(i,j,k)]  = updateTimeRmsfQuantity(rhov_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],rmsf_rhov_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_rhow_field[I1D(i,j,k)]  = updateTimeRmsfQuantity(rhow_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],rmsf_rhow_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_rhoE_field[I1D(i,j,k)]  = updateTimeRmsfQuantity(rhoE_field[I1D(i,j,k)],avg_rhoE_field[I1D(i,j,k)],rmsf_rhoE_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_u_field[I1D(i,j,k)]     = updateTimeRmsfQuantity(u_field[I1D(i,j,k)],avg_u_field[I1D(i,j,k)],rmsf_u_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_v_field[I1D(i,j,k)]     = updateTimeRmsfQuantity(v_field[I1D(i,j,k)],avg_v_field[I1D(i,j,k)],rmsf_v_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_w_field[I1D(i,j,k)]     = updateTimeRmsfQuantity(w_field[I1D(i,j,k)],avg_w_field[I1D(i,j,k)],rmsf_w_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_E_field[I1D(i,j,k)]     = updateTimeRmsfQuantity(E_field[I1D(i,j,k)],avg_E_field[I1D(i,j,k)],rmsf_E_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_P_field[I1D(i,j,k)]     = updateTimeRmsfQuantity(P_field[I1D(i,j,k)],avg_P_field[I1D(i,j,k)],rmsf_P_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_T_field[I1D(i,j,k)]     = updateTimeRmsfQuantity(T_field[I1D(i,j,k)],avg_T_field[I1D(i,j,k)],rmsf_T_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_sos_field[I1D(i,j,k)]   = updateTimeRmsfQuantity(sos_field[I1D(i,j,k)],avg_sos_field[I1D(i,j,k)],rmsf_sos_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_mu_field[I1D(i,j,k)]    = updateTimeRmsfQuantity(mu_field[I1D(i,j,k)],avg_mu_field[I1D(i,j,k)],rmsf_mu_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_kappa_field[I1D(i,j,k)] = updateTimeRmsfQuantity(kappa_field[I1D(i,j,k)],avg_kappa_field[I1D(i,j,k)],rmsf_kappa_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_c_v_field[I1D(i,j,k)]   = updateTimeRmsfQuantity(c_v_field[I1D(i,j,k)],avg_c_v_field[I1D(i,j,k)],rmsf_c_v_field[I1D(i,j,k)],delta_t,averaging_time);
                rmsf_c_p_field[I1D(i,j,k)]   = updateTimeRmsfQuantity(c_p_field[I1D(i,j,k)],avg_c_p_field[I1D(i,j,k)],rmsf_c_p_field[I1D(i,j,k)],delta_t,averaging_time);

                /// Favre-averaged quantities
		favre_uffuff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(u_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],u_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_uffuff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_uffvff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(u_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],v_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_uffvff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_uffwff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(u_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],w_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_uffwff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_vffvff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(v_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],v_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_vffvff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_vffwff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(v_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],w_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_vffwff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_wffwff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(w_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],w_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_wffwff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_uffEff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(u_field[I1D(i,j,k)],avg_rhou_field[I1D(i,j,k)],E_field[I1D(i,j,k)],avg_rhoE_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_uffEff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_vffEff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(v_field[I1D(i,j,k)],avg_rhov_field[I1D(i,j,k)],E_field[I1D(i,j,k)],avg_rhoE_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_vffEff_field[I1D(i,j,k)],delta_t,averaging_time); 
		favre_wffEff_field[I1D(i,j,k)] = updateTimeFavreAveragedQuantity(w_field[I1D(i,j,k)],avg_rhow_field[I1D(i,j,k)],E_field[I1D(i,j,k)],avg_rhoE_field[I1D(i,j,k)],rho_field[I1D(i,j,k)],avg_rho_field[I1D(i,j,k)],favre_wffEff_field[I1D(i,j,k)],delta_t,averaging_time); 
            }
        }
    }

    /// Update averaging time
    averaging_time += delta_t;

    /// Update halo values
    //avg_rho_field.update();
    //avg_rhou_field.update();
    //avg_rhov_field.update();
    //avg_rhow_field.update();
    //avg_rhoE_field.update();
    //avg_rhoP_field.update();
    //avg_rhoT_field.update();
    //avg_u_field.update();
    //avg_v_field.update();
    //avg_w_field.update();
    //avg_E_field.update();
    //avg_P_field.update();
    //avg_T_field.update();
    //avg_sos_field.update();
    //avg_mu_field.update();
    //avg_kappa_field.update();
    //avg_c_v_field.update();
    //avg_c_p_field.update();
    //rmsf_rho_field.update();
    //rmsf_rhou_field.update();
    //rmsf_rhov_field.update();
    //rmsf_rhow_field.update();
    //rmsf_rhoE_field.update();
    //rmsf_u_field.update();
    //rmsf_v_field.update();
    //rmsf_w_field.update();
    //rmsf_E_field.update();
    //rmsf_P_field.update();
    //rmsf_T_field.update();
    //rmsf_sos_field.update();
    //rmsf_mu_field.update();
    //rmsf_kappa_field.update();
    //rmsf_c_v_field.update();
    //rmsf_c_p_field.update();
    //favre_uffuff_field.update();
    //favre_uffvff_field.update();
    //favre_uffwff_field.update();
    //favre_vffvff_field.update();
    //favre_vffwff_field.update();
    //favre_wffwff_field.update();

};

double FlowSolverRHEA::updateTimeMeanQuantity(const double &quantity, const double &mean_quantity, const double &delta_t, const double &averaging_time) {

    double updated_mean_quantity = ( mean_quantity*averaging_time + quantity*delta_t )/( averaging_time + delta_t ); 

    return( updated_mean_quantity );

};

double FlowSolverRHEA::updateTimeRmsfQuantity(const double &quantity, const double &mean_quantity, const double &rmsf_quantity, const double &delta_t, const double &averaging_time) {

    double updated_rmsf_quantity = sqrt( ( pow( rmsf_quantity, 2.0 )*averaging_time + pow( quantity - mean_quantity, 2.0 )*delta_t )/( averaging_time + delta_t ) ); 

    return( updated_rmsf_quantity );

};

//double FlowSolverRHEA::updateTimeReynoldsAveragedQuantity(const double &quantity_1, const double &mean_quantity_1, const double &quantity_2, const double &mean_quantity_2, const double &reynolds_averaged_quantity, const double &delta_t, const double &averaging_time) {
//
//    double fluctuating_quantity_1         = quantity_1 - mean_quantity_1;
//    double fluctuating_quantity_2         = quantity_2 - mean_quantity_2;
//    double product_fluctuating_quantities = fluctuating_quantity_1*fluctuating_quantity_2;
//
//    double updated_reynolds_averaged_quantity = ( reynolds_averaged_quantity*averaging_time + product_fluctuating_quantities*delta_t )/( averaging_time + delta_t ); 
//
//    return( updated_reynolds_averaged_quantity );
//
//};

double FlowSolverRHEA::updateTimeFavreAveragedQuantity(const double &quantity_1, const double &mean_rho_quantity_1, const double &quantity_2, const double &mean_rho_quantity_2, const double &rho, const double &mean_rho, const double &favre_averaged_quantity, const double &delta_t, const double &averaging_time) {

    double old_mean_rho = ( ( averaging_time + delta_t )*mean_rho - rho*delta_t )/( averaging_time + epsilon );

    double favre_mean_quantity_1              = mean_rho_quantity_1/mean_rho;
    double favre_mean_quantity_2              = mean_rho_quantity_2/mean_rho;
    double fluctuating_quantity_1             = quantity_1 - favre_mean_quantity_1;
    double fluctuating_quantity_2             = quantity_2 - favre_mean_quantity_2;
    double product_rho_fluctuating_quantities = rho*fluctuating_quantity_1*fluctuating_quantity_2;

    double updated_favre_averaged_quantity = ( ( old_mean_rho*favre_averaged_quantity*averaging_time + product_rho_fluctuating_quantities*delta_t )/( averaging_time + delta_t ) )/mean_rho; 

    return( updated_favre_averaged_quantity );

};

double FlowSolverRHEA::calculateVolumeAveragedPressure() {

    /// Inner points: P
    double local_sum_VP = 0.0;
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
                /// Sum V*P values
                local_sum_VP += volume*P_field[I1D(i,j,k)];
                /// Sum V values
                local_sum_V += volume;
	    }
        }
    }		    

    /// Communicate local values to obtain global & average values
    double global_sum_VP;
    MPI_Allreduce(&local_sum_VP, &global_sum_VP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_sum_V;
    MPI_Allreduce(&local_sum_V, &global_sum_V, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_avg_P = global_sum_VP/global_sum_V;

    return( global_avg_P );

};

double FlowSolverRHEA::calculateAlphaArtificialCompressibilityMethod() {

    //const double P_threshold = 1.0e-5*P_thermo;
    const double P_threshold = 1.0e-3*P_thermo;

#if 0	/// L1-norm
    /// Inner points: P
    double local_sum_num = 0.0;
    double local_sum_den = 0.0;
    double delta_x, delta_y, delta_z, volume = 0.0;    
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                volume  = delta_x*delta_y*delta_z;
                /// Update values
                local_sum_num += volume*abs( P_field[I1D(i,j,k)] );
                local_sum_den += volume*( max( abs( P_field[I1D(i,j,k)] - P_thermo ), P_threshold ) );
	    }
        }
    }		    

    /// Communicate local values to obtain global values
    double global_sum_num;
    MPI_Allreduce(&local_sum_num, &global_sum_num, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_sum_den;
    MPI_Allreduce(&local_sum_den, &global_sum_den, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_alpha = sqrt( 1.0 + epsilon_acm*global_sum_num/global_sum_den );
#endif

#if 1	/// L2-norm
    /// Inner points: P
    double local_sum_num = 0.0;
    double local_sum_den = 0.0;
    double delta_x, delta_y, delta_z, volume = 0.0;    
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Geometric stuff
                delta_x = 0.5*( x_field[I1D(i+1,j,k)] - x_field[I1D(i-1,j,k)] ); 
                delta_y = 0.5*( y_field[I1D(i,j+1,k)] - y_field[I1D(i,j-1,k)] ); 
                delta_z = 0.5*( z_field[I1D(i,j,k+1)] - z_field[I1D(i,j,k-1)] );
                volume  = delta_x*delta_y*delta_z;
                /// Update values
                local_sum_num += pow( volume*P_field[I1D(i,j,k)], 2.0 );
                local_sum_den += pow( volume*( max( abs( P_field[I1D(i,j,k)] - P_thermo ), P_threshold ) ), 2.0 ); 
	    }
        }
    }		    

    /// Communicate local values to obtain global values
    double global_sum_num;
    MPI_Allreduce(&local_sum_num, &global_sum_num, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_sum_den;
    MPI_Allreduce(&local_sum_den, &global_sum_den, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_alpha = sqrt( 1.0 + epsilon_acm*sqrt( global_sum_num )/sqrt( global_sum_den ) );
#endif

#if 0	/// infinity-norm
    /// Initialize to largest double value
    double local_alpha = numeric_limits<double>::max();

    /// Inner points: P
    double alpha_aux;
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                /// Update value
                alpha_aux   = sqrt( 1.0 + ( P_field[I1D(i,j,k)]*epsilon_acm )/( max( abs( P_field[I1D(i,j,k)] - P_thermo ), P_threshold ) ) );
                local_alpha = min( local_alpha, alpha_aux );
	    }
        }
    }		    

    /// Communicate local value to obtain global value
    double global_alpha;
    MPI_Allreduce(&local_alpha, &global_alpha, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    return( global_alpha );

};

void FlowSolverRHEA::calculateArtificiallyModifiedThermodynamics() {
            
    /// All (inner, halo, boundary) points: sos, c_v, c_p
    double c_v, c_p;
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                sos_field[I1D(i,j,k)]  = ( 1.0/( alpha_acm + epsilon ) )*thermodynamics->calculateSoundSpeed( P_thermo, T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                thermodynamics->calculateSpecificHeatCapacities( c_v, c_p, P_thermo, T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                c_v_field[I1D(i,j,k)]  = c_v;
                c_p_field[I1D(i,j,k)]  = c_p;
            }
        }
    }

    /// Update halo values
    //sos_field.update();
    //c_v_field.update();
    //c_p_field.update();

};

void FlowSolverRHEA::calculateArtificiallyModifiedTransportCoefficients() {
    
    /// All (inner, halo, boundary) points: mu and kappa
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                mu_field[I1D(i,j,k)]    = transport_coefficients->calculateDynamicViscosity( P_thermo, T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
                kappa_field[I1D(i,j,k)] = transport_coefficients->calculateThermalConductivity( P_thermo, T_field[I1D(i,j,k)], rho_field[I1D(i,j,k)] );
            }
        }
    }

    /// Update halo values
    //mu_field.update();
    //kappa_field.update();

};

void FlowSolverRHEA::execute() {
    
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
        
	/// Output temporal point probes data to files
	this->outputTemporalPointProbesData();

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
void FlowSolverRHEA::calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R) {
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

double FlowSolverRHEA::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {

    double F = 0.0;

    /// Select Riemann solver (needs to be coordinated with Riemann solver selector in FlowSolverRHEA constructor )
    if( riemann_solver_scheme_index == 0 ) {	/// KGP

    /// Kennedy, Gruber & Pirozzoli (KGP) scheme:
    /// G. Coppola , F. Capuano , S. Pirozzoli, L. de Luca.
    /// Numerically stable formulations of convective terms for turbulent compressible flows.
    /// Journal of Computational Physics, 382, 86-104, 2019.

    F = ( 1.0/8.0 )*( rho_L + rho_R )*( u_L + u_R );
    if( var_type == 0 ) {
        F *= 1.0 + 1.0;
    } else if ( var_type == 1 ) {
        F *= u_L + u_R; F += ( 1.0/2.0 )*( P_L + P_R );
        //F *= u_L + u_R; F += ( 1.0/4.0 )*( rho_L + rho_R )*( P_L/rho_L + P_R/rho_R );
    } else if ( var_type == 2 ) {
        F *= v_L + v_R;
    } else if ( var_type == 3 ) {
        F *= w_L + w_R;
    } else if ( var_type == 4 ) {
        F *= E_L + P_L/rho_L + E_R + P_R/rho_R;
        //F *= E_L + E_R; F += ( 1.0/4.0 )*( u_L + u_R )*( P_L + P_R );
    }

    } else if( riemann_solver_scheme_index == 1 ) {	/// SHIMA

    /// Shima, Kuya, Tamaki & Kawai (SHIMA) scheme:
    /// N. Shima, Y. Kuya, Y. Tamaki, S. Kawai.
    /// Preventing spurious pressure oscillations in split convective form discretization for compressible flows
    /// Journal of Computational Physics, 427, 110060, 2021.

    F = ( 1.0/8.0 )*( rho_L + rho_R )*( u_L + u_R );
    if( var_type == 0 ) {
        F *= 1.0 + 1.0;
    } else if ( var_type == 1 ) {
        F *= u_L + u_R; F += ( 1.0/2.0 )*( P_L + P_R );
    } else if ( var_type == 2 ) {
        F *= v_L + v_R;
    } else if ( var_type == 3 ) {
        F *= w_L + w_R;
    } else if ( var_type == 4 ) {
        double ke_L = ( 1.0/2.0 )*( u_L*u_L + v_L*v_L + w_L*w_L );
        double e_L  = E_L - ke_L;
        double ke_R = ( 1.0/2.0 )*( u_R*u_R + v_R*v_R + w_R*w_R );
        double e_R  = E_R - ke_R;
        F *= ke_L + ke_R;
        F += ( 1.0/4.0 )*( rho_L*e_L + rho_R*e_R )*( u_L + u_R );
        F += ( 1.0/2.0 )*( u_L*P_R + u_R*P_L );
    }

    } else if( riemann_solver_scheme_index == 2 ) {	/// DIVERGENCE

    /// Divergence scheme obtained from a central differencing of the first derivative of the flux term:

    double F_L = rho_L*u_L;
    double F_R = rho_R*u_R;
    if( var_type == 0 ) {
        F_L *= 1.0;
        F_R *= 1.0;
    } else if ( var_type == 1 ) {
        F_L *= u_L; F_L += P_L;
        F_R *= u_R; F_R += P_R;
    } else if ( var_type == 2 ) {
        F_L *= v_L;
        F_R *= v_R;
    } else if ( var_type == 3 ) {
        F_L *= w_L;
        F_R *= w_R;
    } else if ( var_type == 4 ) {
        F_L *= E_L; F_L += u_L*P_L;
        F_R *= E_R; F_R += u_R*P_R;
    }
    F = 0.5*( F_L + F_R );

    } else if( riemann_solver_scheme_index == 3 ) {	/// MURMAN-ROE

    /// Murman-Roe Riemman solver:
    /// P. L. Roe.
    /// Approximate Riemann solvers, parameter vectors and difference schemes.
    /// Journal of Computational Physics, 43, 357-372, 1981.

    double F_L = rho_L*u_L;
    double F_R = rho_R*u_R;
    double U_L = rho_L;
    double U_R = rho_R;
    if( var_type == 0 ) {
        F_L *= 1.0;
        F_R *= 1.0;
        U_L *= 1.0;
        U_R *= 1.0;
    } else if ( var_type == 1 ) {
        F_L *= u_L; F_L += P_L;
        F_R *= u_R; F_R += P_R;
        U_L *= u_L;
        U_R *= u_R;
    } else if ( var_type == 2 ) {
        F_L *= v_L;
        F_R *= v_R;
        U_L *= v_L;
        U_R *= v_R;
    } else if ( var_type == 3 ) {
        F_L *= w_L;
        F_R *= w_R;
        U_L *= w_L;
        U_R *= w_R;
    } else if ( var_type == 4 ) {
        F_L *= E_L; F_L += u_L*P_L;
        F_R *= E_R; F_R += u_R*P_R;
        U_L *= E_L;
        U_R *= E_R;
    }

    /// Wave speed
    double S = abs( ( F_L - F_R )/( U_L - U_R + epsilon ) );

    /// Conservative + dissipative flux form
    F = 0.5*( F_L + F_R ) - 0.5*S*( U_R - U_L );

    } else if( riemann_solver_scheme_index == 4 ) {	/// HLL

    /// Harten-Lax-van Leer (HLL) Riemman solver:
    /// A. Harten, P. D. Lax, B. van Leer.
    /// On upstream differencing and Godunov-type schemes for hyperbolic conservation laws.
    /// SIAM Review, 25, 35-61, 1983.

    double F_L = rho_L*u_L;
    double F_R = rho_R*u_R;
    double U_L = rho_L;
    double U_R = rho_R;
    if( var_type == 0 ) {
        F_L *= 1.0;
        F_R *= 1.0;
        U_L *= 1.0;
        U_R *= 1.0;
    } else if ( var_type == 1 ) {
        F_L *= u_L; F_L += P_L;
        F_R *= u_R; F_R += P_R;
        U_L *= u_L;
        U_R *= u_R;
    } else if ( var_type == 2 ) {
        F_L *= v_L;
        F_R *= v_R;
        U_L *= v_L;
        U_R *= v_R;
    } else if ( var_type == 3 ) {
        F_L *= w_L;
        F_R *= w_R;
        U_L *= w_L;
        U_R *= w_R;
    } else if ( var_type == 4 ) {
        F_L *= E_L; F_L += u_L*P_L;
        F_R *= E_R; F_R += u_R*P_R;
        U_L *= E_L;
        U_R *= E_R;
    }

    double S_L, S_R;
    calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );

    F = 0.0;
    if( 0.0 <= S_L ) {
        F = F_L;
    } else if( 0.0 >= S_R ) {
        F = F_R;
    } else {
        F = ( S_R*F_L - S_L*F_R + S_L*S_R*( U_R - U_L ) )/( S_R - S_L );
    }

    } else if( riemann_solver_scheme_index == 5 ) {	/// HLLC

    /// Harten-Lax-van Leer-Contact (HLLC) Riemman solver:
    /// E. F. Toro, M. Spruce, W. Speares.
    /// Restoration of the contact surface in the HLL-Riemann solver.
    /// Shock Waves, 4, 25-34, 1994.

    double F_L = rho_L*u_L;
    double F_R = rho_R*u_R;
    double U_L = rho_L;
    double U_R = rho_R;
    if( var_type == 0 ) {
        F_L *= 1.0;
        F_R *= 1.0;
        U_L *= 1.0;
        U_R *= 1.0;
    } else if ( var_type == 1 ) {
        F_L *= u_L; F_L += P_L;
        F_R *= u_R; F_R += P_R;
        U_L *= u_L;
        U_R *= u_R;
    } else if ( var_type == 2 ) {
        F_L *= v_L;
        F_R *= v_R;
        U_L *= v_L;
        U_R *= v_R;
    } else if ( var_type == 3 ) {
        F_L *= w_L;
        F_R *= w_R;
        U_L *= w_L;
        U_R *= w_R;
    } else if ( var_type == 4 ) {
        F_L *= E_L; F_L += u_L*P_L;
        F_R *= E_R; F_R += u_R*P_R;
        U_L *= E_L;
        U_R *= E_R;
    }

    double S_L, S_R;
    calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );

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

    F = 0.0;
    if( 0.0 <= S_L ) {
        F = F_L;
    } else if( ( S_L <= 0.0 ) && ( 0.0 <= S_star ) ) {
        F = F_L + S_L*( U_star_L - U_L );
    } else if( ( S_star <= 0.0 ) && ( 0.0 <= S_R ) ) {
        F = F_R + S_R*( U_star_R - U_R );
    } else if( 0.0 >= S_R ) {
        F = F_R;
    }

    } else if( riemann_solver_scheme_index == 6 ) {	/// HLLC+

    /// HLLC-type Riemann solver for all-speed flows:
    /// S. Chen, B. Lin, Y. Li, C. Yan.
    /// HLLC+: low-Mach shock-stable HLLC-type Riemann solver for all-speed flows.
    /// SIAM Journal of Scientific Computing, 4, B921-B950, 2020.

    double F_L = rho_L*u_L;
    double F_R = rho_R*u_R;
    double U_L = rho_L;
    double U_R = rho_R;
    if( var_type == 0 ) {
        F_L *= 1.0;
        F_R *= 1.0;
        U_L *= 1.0;
        U_R *= 1.0;
    } else if ( var_type == 1 ) {
        F_L *= u_L; F_L += P_L;
        F_R *= u_R; F_R += P_R;
        U_L *= u_L;
        U_R *= u_R;
    } else if ( var_type == 2 ) {
        F_L *= v_L;
        F_R *= v_R;
        U_L *= v_L;
        U_R *= v_R;
    } else if ( var_type == 3 ) {
        F_L *= w_L;
        F_R *= w_R;
        U_L *= w_L;
        U_R *= w_R;
    } else if ( var_type == 4 ) {
        F_L *= E_L; F_L += u_L*P_L;
        F_R *= E_R; F_R += u_R*P_R;
        U_L *= E_L;
        U_R *= E_R;
    }

    double S_L, S_R;
    calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );

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

    F = 0.0;
    if( 0.0 <= S_L ) {
        F = F_L;
    } else if( ( S_L <= 0.0 ) && ( 0.0 <= S_star ) ) {
        F = F_star_L;
    } else if( ( S_star <= 0.0 ) && ( 0.0 <= S_R ) ) {
        F = F_star_R;
    } else if( 0.0 >= S_R ) {
        F = F_R;
    }
    
    }

    return( F );

};


//////////// BaseRiemannSolver CLASS //////////
//
//BaseRiemannSolver::BaseRiemannSolver() {};
//        
//BaseRiemannSolver::~BaseRiemannSolver() {};
//
//void BaseRiemannSolver::calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R) {
//
//#if _PRESSURE_BASED_WAVE_SPEED_ESTIMATES_
//    /// Pressure-based wave speed estimates: ... recommended by E. F. Toro, but not sure if applicable to non-ideal gas themodynamics
//    /// E. F. Toro, M. Spruce, W. Speares.
//    /// Restoration of the contact surface in the HLL-Riemann solver.
//    /// Shock Waves, 4, 25-34, 1994.
//
//    double P_bar   = 0.5*( P_L + P_R );
//    double rho_bar = 0.5*( rho_L + rho_R );
//    double gamma   = thermodynamics->calculateHeatCapacitiesRatio( P_bar, rho_bar );
//    double a_bar   = 0.5*( a_L + a_R );
//    double P_pvrs  = 0.5*( P_L + P_R ) - 0.5*( u_R - u_L )*rho_bar*a_bar;
//    double P_star  = max( 0.0, P_pvrs );
//    double q_L     = 1.0;
//    if(P_star > P_L) q_L = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_L ) - 1.0 ) );
//    double q_R     = 1.0;
//    if(P_star > P_R) q_R = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_R ) - 1.0 ) );
//    S_L = u_L - a_L*q_L;
//    S_R = u_R + a_R*q_R;
//#else
//    /// Direct wave speed estimates:
//    /// B. Einfeldt.
//    /// On Godunov-type methods for gas dynamics.
//    /// SIAM Journal on Numerical Analysis, 25, 294-318, 1988.
//
//    double hat_u = ( u_L*sqrt( rho_L ) + u_R*sqrt( rho_R ) )/( sqrt( rho_L ) + sqrt( rho_R ) );
//    double hat_a = sqrt( ( ( a_L*a_L*sqrt( rho_L ) + a_R*a_R*sqrt( rho_R ) )/( sqrt( rho_L ) + sqrt( rho_R ) ) ) + 0.5*( ( sqrt( rho_L )*sqrt( rho_R ) )/( ( sqrt( rho_L ) + sqrt( rho_R ) )*( sqrt( rho_L ) + sqrt( rho_R ) ) ) )*( u_R - u_L )*( u_R - u_L ) );
//
//    S_L = min( u_L - a_L, hat_u - hat_a );
//    S_R = max( u_R + a_R, hat_u + hat_a );
//#endif
//
//};


//////////// DivergenceFluxApproximateRiemannSolver CLASS //////////
//
//DivergenceFluxApproximateRiemannSolver::DivergenceFluxApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//DivergenceFluxApproximateRiemannSolver::~DivergenceFluxApproximateRiemannSolver() {};
//
//double DivergenceFluxApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// Divergence scheme obtained from a central differencing of the first derivative of the flux term:
//
//    double F_L = rho_L*u_L;
//    double F_R = rho_R*u_R;
//    if( var_type == 0 ) {
//        F_L *= 1.0;
//        F_R *= 1.0;
//    } else if ( var_type == 1 ) {
//        F_L *= u_L; F_L += P_L;
//        F_R *= u_R; F_R += P_R;
//    } else if ( var_type == 2 ) {
//        F_L *= v_L;
//        F_R *= v_R;
//    } else if ( var_type == 3 ) {
//        F_L *= w_L;
//        F_R *= w_R;
//    } else if ( var_type == 4 ) {
//        F_L *= E_L; F_L += u_L*P_L;
//        F_R *= E_R; F_R += u_R*P_R;
//    }
//    double F = 0.5*( F_L + F_R );
//
//    return( F );
//
//};


//////////// MurmanRoeFluxApproximateRiemannSolver CLASS //////////
//
//MurmanRoeFluxApproximateRiemannSolver::MurmanRoeFluxApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//MurmanRoeFluxApproximateRiemannSolver::~MurmanRoeFluxApproximateRiemannSolver() {};
//
//double MurmanRoeFluxApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// Murman-Roe Riemman solver:
//    /// P. L. Roe.
//    /// Approximate Riemann solvers, parameter vectors and difference schemes.
//    /// Journal of Computational Physics, 43, 357-372, 1981.
//
//    double F_L = rho_L*u_L;
//    double F_R = rho_R*u_R;
//    double U_L = rho_L;
//    double U_R = rho_R;
//    if( var_type == 0 ) {
//        F_L *= 1.0;
//        F_R *= 1.0;
//        U_L *= 1.0;
//        U_R *= 1.0;
//    } else if ( var_type == 1 ) {
//        F_L *= u_L; F_L += P_L;
//        F_R *= u_R; F_R += P_R;
//        U_L *= u_L;
//        U_R *= u_R;
//    } else if ( var_type == 2 ) {
//        F_L *= v_L;
//        F_R *= v_R;
//        U_L *= v_L;
//        U_R *= v_R;
//    } else if ( var_type == 3 ) {
//        F_L *= w_L;
//        F_R *= w_R;
//        U_L *= w_L;
//        U_R *= w_R;
//    } else if ( var_type == 4 ) {
//        F_L *= E_L; F_L += u_L*P_L;
//        F_R *= E_R; F_R += u_R*P_R;
//        U_L *= E_L;
//        U_R *= E_R;
//    }
//
//    /// Wave speed
//    double S = abs( ( F_L - F_R )/( U_L - U_R + epsilon ) );
//
//    /// Conservative + dissipative flux form
//    double F = 0.5*( F_L + F_R ) - 0.5*S*( U_R - U_L );
//
//    return( F );
//
//};


//////////// KgpFluxApproximateRiemannSolver CLASS //////////
//
//KgpFluxApproximateRiemannSolver::KgpFluxApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//KgpFluxApproximateRiemannSolver::~KgpFluxApproximateRiemannSolver() {};
//
//double KgpFluxApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// Kennedy, Gruber & Pirozzoli (KGP) scheme:
//    /// G. Coppola, F. Capuano , S. Pirozzoli, L. de Luca.
//    /// Numerically stable formulations of convective terms for turbulent compressible flows.
//    /// Journal of Computational Physics, 382, 86-104, 2019.
//
//    double F = ( 1.0/8.0 )*( rho_L + rho_R )*( u_L + u_R );
//    if( var_type == 0 ) {
//        F *= 1.0 + 1.0;
//    } else if ( var_type == 1 ) {
//        F *= u_L + u_R; F += ( 1.0/2.0 )*( P_L + P_R );
//        //F *= u_L + u_R; F += ( 1.0/4.0 )*( rho_L + rho_R )*( P_L/rho_L + P_R/rho_R );
//    } else if ( var_type == 2 ) {
//        F *= v_L + v_R;
//    } else if ( var_type == 3 ) {
//        F *= w_L + w_R;
//    } else if ( var_type == 4 ) {
//        F *= E_L + P_L/rho_L + E_R + P_R/rho_R;
//        //F *= E_L + E_R; F += ( 1.0/4.0 )*( u_L + u_R )*( P_L + P_R );
//    }
//
//    return( F );
//
//};


//////////// ShimaFluxApproximateRiemannSolver CLASS //////////
//
//ShimaFluxApproximateRiemannSolver::ShimaFluxApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//ShimaFluxApproximateRiemannSolver::~ShimaFluxApproximateRiemannSolver() {};
//
//double ShimaFluxApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// Shima, Kuya, Tamaki & Kawai (SHIMA) scheme:
//    /// N. Shima, Y. Kuya, Y. Tamaki, S. Kawai.
//    /// Preventing spurious pressure oscillations in split convective form discretization for compressible flows
//    /// Journal of Computational Physics, 427, 110060, 2021.
//    
//    double F = ( 1.0/8.0 )*( rho_L + rho_R )*( u_L + u_R );
//    if( var_type == 0 ) {
//        F *= 1.0 + 1.0;
//    } else if ( var_type == 1 ) {
//        F *= u_L + u_R; F += ( 1.0/2.0 )*( P_L + P_R );
//    } else if ( var_type == 2 ) {
//        F *= v_L + v_R;
//    } else if ( var_type == 3 ) {
//        F *= w_L + w_R;
//    } else if ( var_type == 4 ) {
//        double ke_L = ( 1.0/2.0 )*( u_L*u_L + v_L*v_L + w_L*w_L );
//        double e_L  = E_L - ke_L;
//        double ke_R = ( 1.0/2.0 )*( u_R*u_R + v_R*v_R + w_R*w_R );
//        double e_R  = E_R - ke_R;
//        F *= ke_L + ke_R;
//        F += ( 1.0/4.0 )*( rho_L*e_L + rho_R*e_R )*( u_L + u_R );
//        F += ( 1.0/2.0 )*( u_L*P_R + u_R*P_L );
//    }
//
//    return( F );
//
//    return( F );
//
//};


//////////// HllApproximateRiemannSolver CLASS //////////
//
//HllApproximateRiemannSolver::HllApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//HllApproximateRiemannSolver::~HllApproximateRiemannSolver() {};
//
//double HllApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// Harten-Lax-van Leer (HLL) Riemman solver:
//    /// A. Harten, P. D. Lax, B. van Leer.
//    /// On upstream differencing and Godunov-type schemes for hyperbolic conservation laws.
//    /// SIAM Review, 25, 35-61, 1983.
//
//    double F_L = rho_L*u_L;
//    double F_R = rho_R*u_R;
//    double U_L = rho_L;
//    double U_R = rho_R;
//    if( var_type == 0 ) {
//        F_L *= 1.0;
//        F_R *= 1.0;
//        U_L *= 1.0;
//        U_R *= 1.0;
//    } else if ( var_type == 1 ) {
//        F_L *= u_L; F_L += P_L;
//        F_R *= u_R; F_R += P_R;
//        U_L *= u_L;
//        U_R *= u_R;
//    } else if ( var_type == 2 ) {
//        F_L *= v_L;
//        F_R *= v_R;
//        U_L *= v_L;
//        U_R *= v_R;
//    } else if ( var_type == 3 ) {
//        F_L *= w_L;
//        F_R *= w_R;
//        U_L *= w_L;
//        U_R *= w_R;
//    } else if ( var_type == 4 ) {
//        F_L *= E_L; F_L += u_L*P_L;
//        F_R *= E_R; F_R += u_R*P_R;
//        U_L *= E_L;
//        U_R *= E_R;
//    }
//
//    double S_L, S_R;
//    this->calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );
//
//    double F = 0.0;
//    if( 0.0 <= S_L ) {
//        F = F_L;
//    } else if( 0.0 >= S_R ) {
//        F = F_R;
//    } else {
//        F = ( S_R*F_L - S_L*F_R + S_L*S_R*( U_R - U_L ) )/( S_R - S_L );
//    }
//
//    return( F );
//
//};


//////////// HllcApproximateRiemannSolver CLASS //////////
//
//HllcApproximateRiemannSolver::HllcApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//HllcApproximateRiemannSolver::~HllcApproximateRiemannSolver() {};
//
//double HllcApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// Harten-Lax-van Leer-Contact (HLLC) Riemman solver:
//    /// E. F. Toro, M. Spruce, W. Speares.
//    /// Restoration of the contact surface in the HLL-Riemann solver.
//    /// Shock Waves, 4, 25-34, 1994.
//
//    double F_L = rho_L*u_L;
//    double F_R = rho_R*u_R;
//    double U_L = rho_L;
//    double U_R = rho_R;
//    if( var_type == 0 ) {
//        F_L *= 1.0;
//        F_R *= 1.0;
//        U_L *= 1.0;
//        U_R *= 1.0;
//    } else if ( var_type == 1 ) {
//        F_L *= u_L; F_L += P_L;
//        F_R *= u_R; F_R += P_R;
//        U_L *= u_L;
//        U_R *= u_R;
//    } else if ( var_type == 2 ) {
//        F_L *= v_L;
//        F_R *= v_R;
//        U_L *= v_L;
//        U_R *= v_R;
//    } else if ( var_type == 3 ) {
//        F_L *= w_L;
//        F_R *= w_R;
//        U_L *= w_L;
//        U_R *= w_R;
//    } else if ( var_type == 4 ) {
//        F_L *= E_L; F_L += u_L*P_L;
//        F_R *= E_R; F_R += u_R*P_R;
//        U_L *= E_L;
//        U_R *= E_R;
//    }
//
//    double S_L, S_R;
//    this->calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );
//
//    double S_star   = ( P_R - P_L + rho_L*u_L*( S_L - u_L ) - rho_R*u_R*( S_R - u_R ) )/( rho_L*( S_L - u_L ) - rho_R*( S_R - u_R ) );
//    double U_star_L = rho_L*( ( S_L - u_L )/( S_L - S_star ) );
//    double U_star_R = rho_R*( ( S_R - u_R )/( S_R - S_star ) );
//    if( var_type == 0 ) {
//        U_star_L *= 1.0;
//        U_star_R *= 1.0;       
//    } else if( var_type == 1 ) {
//        U_star_L *= S_star;
//        U_star_R *= S_star;
//    } else if( var_type == 2 ) {
//        U_star_L *= v_L;
//        U_star_R *= v_R;
//    } else if( var_type == 3 ) {
//        U_star_L *= w_L;
//        U_star_R *= w_R;
//    } else if( var_type == 4 ) {
//        U_star_L *= ( E_L + ( S_star - u_L )*( S_star + P_L/( rho_L*( S_L - u_L ) ) ) );
//        U_star_R *= ( E_R + ( S_star - u_R )*( S_star + P_R/( rho_R*( S_R - u_R ) ) ) );
//    }
//
//    double F = 0.0;
//    if( 0.0 <= S_L ) {
//        F = F_L;
//    } else if( ( S_L <= 0.0 ) && ( 0.0 <= S_star ) ) {
//        F = F_L + S_L*( U_star_L - U_L );
//    } else if( ( S_star <= 0.0 ) && ( 0.0 <= S_R ) ) {
//        F = F_R + S_R*( U_star_R - U_R );
//    } else if( 0.0 >= S_R ) {
//        F = F_R;
//    }
//
//    return( F );
//
//};


//////////// HllcPlusApproximateRiemannSolver CLASS //////////
//
//HllcPlusApproximateRiemannSolver::HllcPlusApproximateRiemannSolver() : BaseRiemannSolver() {};
//
//HllcPlusApproximateRiemannSolver::~HllcPlusApproximateRiemannSolver() {};
//
//double HllcPlusApproximateRiemannSolver::calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) {
//
//    /// HLLC-type Riemann solver for all-speed flows:
//    /// S. Chen, B. Lin, Y. Li, C. Yan.
//    /// HLLC+: low-Mach shock-stable HLLC-type Riemann solver for all-speed flows.
//    /// SIAM Journal of Scientific Computing, 4, B921-B950, 2020.
//
//    double F_L = rho_L*u_L;
//    double F_R = rho_R*u_R;
//    double U_L = rho_L;
//    double U_R = rho_R;
//    if( var_type == 0 ) {
//        F_L *= 1.0;
//        F_R *= 1.0;
//        U_L *= 1.0;
//        U_R *= 1.0;
//    } else if ( var_type == 1 ) {
//        F_L *= u_L; F_L += P_L;
//        F_R *= u_R; F_R += P_R;
//        U_L *= u_L;
//        U_R *= u_R;
//    } else if ( var_type == 2 ) {
//        F_L *= v_L;
//        F_R *= v_R;
//        U_L *= v_L;
//        U_R *= v_R;
//    } else if ( var_type == 3 ) {
//        F_L *= w_L;
//        F_R *= w_R;
//        U_L *= w_L;
//        U_R *= w_R;
//    } else if ( var_type == 4 ) {
//        F_L *= E_L; F_L += u_L*P_L;
//        F_R *= E_R; F_R += u_R*P_R;
//        U_L *= E_L;
//        U_R *= E_R;
//    }
//
//    double S_L, S_R;
//    this->calculateWavesSpeed( S_L, S_R, rho_L, rho_R, u_L, u_R, P_L, P_R, a_L, a_R );
//
//    double phi_L = rho_L*( S_L - u_L );
//    double phi_R = rho_R*( S_R - u_R );
//    double S_star   = ( P_R - P_L + phi_L*u_L - phi_R*u_R )/( phi_L - phi_R );
//    double U_star_L = rho_L*( ( S_L - u_L )/( S_L - S_star ) );
//    double U_star_R = rho_R*( ( S_R - u_R )/( S_R - S_star ) );
//    double M        = min( 1.0, max( ( 1.0/a_L )*sqrt( u_L*u_L + v_L*v_L + w_L*w_L ), ( 1.0/a_R )*sqrt( u_R*u_R + v_R*v_R + w_R*w_R ) ) );
//    double f_M      = M*sqrt( 4.0 + pow( 1.0 - M*M, 2.0 ) )/( 1.0 + M*M );
//    double h        = min( P_L/P_R, P_R/P_L );
//    double g        = 1.0 - pow( h, M );
//    double A_p_L    = ( phi_L*phi_R )/( phi_R - phi_L );
//    double A_p_R    = ( phi_L*phi_R )/( phi_R - phi_L );
//    if( var_type == 0 ) {
//        U_star_L *= 1.0;
//        U_star_R *= 1.0;       
//        A_p_L    *= 0.0;
//        A_p_R    *= 0.0;
//    } else if( var_type == 1 ) {
//        U_star_L *= S_star;
//        U_star_R *= S_star;
//        A_p_L    *= ( f_M - 1.0 )*( u_R - u_L );
//        A_p_R    *= ( f_M - 1.0 )*( u_R - u_L );
//    } else if( var_type == 2 ) {
//        U_star_L *= v_L;
//        U_star_R *= v_R;
//        A_p_L    *= ( S_L/( S_L - S_star ) )*g*( v_R - v_L );
//        A_p_R    *= ( S_R/( S_R - S_star ) )*g*( v_R - v_L );
//    } else if( var_type == 3 ) {
//        U_star_L *= w_L;
//        U_star_R *= w_R;
//        A_p_L    *= ( S_L/( S_L - S_star ) )*g*( w_R - w_L );
//        A_p_R    *= ( S_R/( S_R - S_star ) )*g*( w_R - w_L );
//    } else if( var_type == 4 ) {
//        U_star_L *= ( E_L + ( S_star - u_L )*( S_star + P_L/( rho_L*( S_L - u_L ) ) ) );
//        U_star_R *= ( E_R + ( S_star - u_R )*( S_star + P_R/( rho_R*( S_R - u_R ) ) ) );
//        A_p_L    *= ( f_M - 1.0 )*( u_R - u_L )*S_star;
//        A_p_R    *= ( f_M - 1.0 )*( u_R - u_L )*S_star;
//    }
//    double F_star_L = F_L + S_L*( U_star_L - U_L ) + A_p_L;
//    double F_star_R = F_R + S_R*( U_star_R - U_R ) + A_p_R;
//
//    double F = 0.0;
//    if( 0.0 <= S_L ) {
//        F = F_L;
//    } else if( ( S_L <= 0.0 ) && ( 0.0 <= S_star ) ) {
//        F = F_star_L;
//    } else if( ( S_star <= 0.0 ) && ( 0.0 <= S_R ) ) {
//        F = F_star_R;
//    } else if( 0.0 >= S_R ) {
//        F = F_R;
//    }
//
//    return( F );
//
//};


////////// BaseExplicitRungeKuttaMethod CLASS //////////

BaseExplicitRungeKuttaMethod::BaseExplicitRungeKuttaMethod() {};
        
BaseExplicitRungeKuttaMethod::~BaseExplicitRungeKuttaMethod() {};


////////// RungeKutta1Method CLASS //////////

RungeKutta1Method::RungeKutta1Method() : BaseExplicitRungeKuttaMethod() {};

RungeKutta1Method::~RungeKutta1Method() {};

void RungeKutta1Method::setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage) {

    /// Explicit first-order Runge-Kutta (RK1) method:
    /// S. Gottlieb, C.-W. Shu & E. Tadmor.
    /// Strong stability-preserving high-order time discretization methods.
    /// SIAM Review 43, 89-112, 2001.

    /// First Runge-Kutta stage
    rk_a = 1.0; rk_b = 0.0; rk_c = 1.0;

};


////////// StrongStabilityPreservingRungeKutta2Method CLASS //////////

StrongStabilityPreservingRungeKutta2Method::StrongStabilityPreservingRungeKutta2Method() : BaseExplicitRungeKuttaMethod() {};

StrongStabilityPreservingRungeKutta2Method::~StrongStabilityPreservingRungeKutta2Method() {};

void StrongStabilityPreservingRungeKutta2Method::setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage) {

    /// Explicit second-order strong-stability-preserving Runge-Kutta (SSP-RK2) method:
    /// S. Gottlieb, C.-W. Shu & E. Tadmor.
    /// Strong stability-preserving high-order time discretization methods.
    /// SIAM Review 43, 89-112, 2001.

    if(rk_time_stage == 1) {
        /// First Runge-Kutta stage
        rk_a = 1.0; rk_b = 0.0; rk_c = 1.0;
    } else if(rk_time_stage == 2) {
        /// Second Runge-Kutta stage
        rk_a = 1.0/2.0; rk_b = 1.0/2.0; rk_c = 1.0/2.0;
    }

};


////////// StrongStabilityPreservingRungeKutta3Method CLASS //////////

StrongStabilityPreservingRungeKutta3Method::StrongStabilityPreservingRungeKutta3Method() : BaseExplicitRungeKuttaMethod() {};

StrongStabilityPreservingRungeKutta3Method::~StrongStabilityPreservingRungeKutta3Method() {};

void StrongStabilityPreservingRungeKutta3Method::setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage) {

    /// Explicit third-order strong-stability-preserving Runge-Kutta (SSP-RK3) method:
    /// S. Gottlieb, C.-W. Shu & E. Tadmor.
    /// Strong stability-preserving high-order time discretization methods.
    /// SIAM Review 43, 89-112, 2001.

    if(rk_time_stage == 1) {
        /// First Runge-Kutta stage
        rk_a = 1.0; rk_b = 0.0; rk_c = 1.0;
    } else if(rk_time_stage == 2) {
        /// Second Runge-Kutta stage
        rk_a = 3.0/4.0; rk_b = 1.0/4.0; rk_c = 1.0/4.0;
    } else if(rk_time_stage == 3) {
        /// Third Runge-Kutta stage
        rk_a = 1.0/3.0; rk_b = 2.0/3.0; rk_c = 2.0/3.0;
    }

};

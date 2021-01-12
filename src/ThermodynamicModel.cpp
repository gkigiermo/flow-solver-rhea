#include "ThermodynamicModel.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////


////////// BaseThermodynamicModel CLASS //////////

BaseThermodynamicModel::BaseThermodynamicModel() {};
        
BaseThermodynamicModel::BaseThermodynamicModel(const string configuration_file) : configuration_file(configuration_file) {};

BaseThermodynamicModel::~BaseThermodynamicModel() {};


////////// IdealGasModel CLASS //////////

IdealGasModel::IdealGasModel() : BaseThermodynamicModel() {};
        
IdealGasModel::IdealGasModel(const string configuration_file) : BaseThermodynamicModel(configuration_file) {

    /// Read configuration (input) file
    this->readConfigurationFile();

};

IdealGasModel::~IdealGasModel() {};

void IdealGasModel::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid & flow properties
    const YAML::Node & fluid_flow_properties = configuration["fluid_flow_properties"];
    R_specific = fluid_flow_properties["R_specific"].as<double>();
    gamma      = fluid_flow_properties["gamma"].as<double>();

};

double IdealGasModel::calculateTemperatureFromPressureDensity(const double &P, const double &rho) {

    /// Ideal-gas model:
    /// P = rho*R_specific*T is pressure

    double T = P/( rho*R_specific );

    return( T );

};

void IdealGasModel::calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {

    /// Ideal-gas model:
    /// P = e*rho*(gamma - 1) is pressure
    /// T = e/c_v is temperature

    double c_v = R_specific/( gamma - 1.0 );

    P = e*rho*( gamma - 1.0 ); 
    T = e/c_v;

};

void IdealGasModel::calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {

    /// Ideal-gas model:
    /// rho = P/( e*( gamma - 1.0 ) ) is density
    /// e = c_v*T is specific internal energy

    double c_v = R_specific/( gamma - 1.0 );

    e   = c_v*T;
    rho = P/( e*( gamma - 1.0 ) );

};

void IdealGasModel::calculateSpecificHeatCapacities(double &c_v, double &c_p, const double &P, const double &T, const double &rho) {

    /// Ideal-gas model:
    /// c_v = R_specific/(gamma - 1)
    /// c_p = c_v*gamma

    c_v = R_specific/( gamma - 1.0 );
    c_p = c_v*gamma;

};

double IdealGasModel::calculateHeatCapacitiesRatio(const double &P, const double &rho) {

    /// Ideal-gas model:
    /// gamma = c_p/c_v

    return( gamma );

};

double IdealGasModel::calculateSoundSpeed(const double &P, const double &T, const double &rho) {

    /// Ideal-gas model:
    /// sos = sqrt(gamma*P/rho) is speed of sound

    double sos = sqrt( gamma*P/rho );

    return( sos );

};


////////// StiffenedGasModel CLASS //////////

StiffenedGasModel::StiffenedGasModel() : BaseThermodynamicModel() {};
        
StiffenedGasModel::StiffenedGasModel(const string configuration_file) : BaseThermodynamicModel(configuration_file) {

    /// Read configuration (input) file
    this->readConfigurationFile();

};

StiffenedGasModel::~StiffenedGasModel() {};

void StiffenedGasModel::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid & flow properties
    const YAML::Node & fluid_flow_properties = configuration["fluid_flow_properties"];
    R_specific = fluid_flow_properties["R_specific"].as<double>();
    gamma      = fluid_flow_properties["gamma"].as<double>();
    P_inf      = fluid_flow_properties["P_inf"].as<double>();
    e_0        = fluid_flow_properties["e_0"].as<double>();
    c_v        = fluid_flow_properties["c_v"].as<double>();

};

double StiffenedGasModel::calculateTemperatureFromPressureDensity(const double &P, const double &rho) {

    /// Stiffened-gas model:
    /// P = (e - e_0)*rho*(gamma - 1) - gamma*P_inf is pressure
    /// T = ((e - e_0) - (P_inf/rho))/c_v is temperature

    double T = ( ( ( P + gamma*P_inf )/( rho*( gamma - 1.0 ) ) ) - ( P_inf/rho ) )/c_v;

    return( T );

};

void StiffenedGasModel::calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {

    /// Stiffened-gas model:
    /// P = (e - e_0)*rho*(gamma - 1) - gamma*P_inf is pressure
    /// T = ((e - e_0) - (P_inf/rho))/c_v is temperature
    
    P = ( e - e_0 )*rho*( gamma - 1.0 ) - gamma*P_inf; 
    T = ( ( e - e_0 ) - ( P_inf/rho ) )/c_v;

};

void StiffenedGasModel::calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {

    /// Stiffened-gas model:
    /// rho = P_inf/(e - e_0 - c_v*T) is density
    /// e = c_v*T*((P + gamma*P_inf)/(P + P_inf)) + e_0 is specific internal energy

    e   = c_v*T*( ( P + gamma*P_inf )/( P + P_inf ) ) + e_0;
    rho = P_inf/( e - e_0 - c_v*T );

};

void StiffenedGasModel::calculateSpecificHeatCapacities(double &c_v_, double &c_p, const double &P, const double &T, const double &rho) {

    /// Stiffened-gas model:
    /// c_p = c_v*gamma

    c_v_ = c_v;
    c_p  = c_v*gamma;

};

double StiffenedGasModel::calculateHeatCapacitiesRatio(const double &P, const double &rho) {

    /// Stiffened-gas model:
    /// gamma = c_p/c_v

    return( gamma );

};

double StiffenedGasModel::calculateSoundSpeed(const double &P, const double &T, const double &rho) {

    /// Stiffened-gas model:
    /// sos = sqrt(gamma*(P+P_inf)/rho) is speed of sound

    double sos = sqrt( gamma*( P + P_inf )/rho );

    return( sos );

};


////////// PengRobinsonModel CLASS //////////

PengRobinsonModel::PengRobinsonModel() : BaseThermodynamicModel() {};
        
PengRobinsonModel::PengRobinsonModel(const string configuration_file) : BaseThermodynamicModel(configuration_file) {

    /// Read configuration (input) file
    this->readConfigurationFile();

    /// Calculate value of selected variables
    eos_b  = 0.077796*( R_universal*critical_temperature/critical_pressure );
    eos_ac = 0.457236*( pow( R_universal*critical_temperature, 2.0 )/critical_pressure );
    if( acentric_factor > 0.49 ) {
    	eos_kappa = 0.379642 + 1.48503*acentric_factor - 0.164423*pow( acentric_factor, 2.0 ) + 0.016666*pow( acentric_factor, 3.0 );
    } else {
	eos_kappa = 0.37464 + 1.54226*acentric_factor - 0.26992*pow( acentric_factor, 2.0 );
    }

    /// Construct (initialize) nr_PT_solver
    nr_PT_unknowns.resize( 2, 0.0 );
    nr_PT_r_vec.resize( nr_PT_unknowns.size(), 0.0 );
    nr_PT_solver = new NR_P_T_from_rho_e( nr_PT_r_vec, *this );

};

PengRobinsonModel::~PengRobinsonModel() {

    /// Free nr_PT_solver
    if( nr_PT_solver != NULL ) free( nr_PT_solver );

};

void PengRobinsonModel::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid & flow properties
    const YAML::Node & fluid_flow_properties = configuration["fluid_flow_properties"];
    R_specific            = fluid_flow_properties["R_specific"].as<double>();
    molecular_weight      = fluid_flow_properties["molecular_weight"].as<double>();
    acentric_factor       = fluid_flow_properties["acentric_factor"].as<double>();
    critical_temperature  = fluid_flow_properties["critical_temperature"].as<double>();
    critical_pressure     = fluid_flow_properties["critical_pressure"].as<double>();
    critical_molar_volume = fluid_flow_properties["critical_molar_volume"].as<double>();
    NASA_coefficients[0]  = fluid_flow_properties["NASA_coefficients"][0].as<double>();
    NASA_coefficients[1]  = fluid_flow_properties["NASA_coefficients"][1].as<double>();
    NASA_coefficients[2]  = fluid_flow_properties["NASA_coefficients"][2].as<double>();
    NASA_coefficients[3]  = fluid_flow_properties["NASA_coefficients"][3].as<double>();
    NASA_coefficients[4]  = fluid_flow_properties["NASA_coefficients"][4].as<double>();
    NASA_coefficients[5]  = fluid_flow_properties["NASA_coefficients"][5].as<double>();
    NASA_coefficients[6]  = fluid_flow_properties["NASA_coefficients"][6].as<double>();
    NASA_coefficients[7]  = fluid_flow_properties["NASA_coefficients"][7].as<double>();
    NASA_coefficients[8]  = fluid_flow_properties["NASA_coefficients"][8].as<double>();
    NASA_coefficients[9]  = fluid_flow_properties["NASA_coefficients"][9].as<double>();
    NASA_coefficients[10] = fluid_flow_properties["NASA_coefficients"][10].as<double>();
    NASA_coefficients[11] = fluid_flow_properties["NASA_coefficients"][11].as<double>();
    NASA_coefficients[12] = fluid_flow_properties["NASA_coefficients"][12].as<double>();
    NASA_coefficients[13] = fluid_flow_properties["NASA_coefficients"][13].as<double>();
    NASA_coefficients[14] = fluid_flow_properties["NASA_coefficients"][14].as<double>();

};

double PengRobinsonModel::calculateTemperatureFromPressureDensity(const double &P, const double &rho) {

    /// Numerical Recipes in C++, Second Edition.
    /// W.H. Press, S.A. Teulosky, W.T. Vetterling, B.P. Flannery.
    /// 5.1 Series and Their Convergence: Aitken’s delta-squared process.

    // Calculate molar volume
    double bar_v = molecular_weight/rho;

    // Initial temperature guess using ideal-gas model
    double T = P*bar_v/R_universal;

    /// Aitken’s delta-squared process:
    double x_0 = T, x_1, x_2, denominator;
    for(int iter = 0; iter < max_aitken_iter; iter++) { 
        x_1 = ( (bar_v - eos_b )/R_universal )*( P + ( this->calculate_eos_a( x_0 )/( pow( bar_v, 2.0 ) + 2.0*eos_b*bar_v - pow( eos_b, 2.0 ) ) ) );
        x_2 = ( (bar_v - eos_b )/R_universal )*( P + ( this->calculate_eos_a( x_1 )/( pow( bar_v, 2.0 ) + 2.0*eos_b*bar_v - pow( eos_b, 2.0 ) ) ) );

        denominator = x_2 - 2.0*x_1 + x_0;
        T = x_2 - ( pow( x_2 - x_1, 2.0 ) )/denominator;
    
        if( abs( ( T - x_2 )/ T ) < aitken_relative_tolerance ) break;	/// If the result is within tolerance, leave the loop!
        x_0 = T;							/// Otherwise, update x_0 to iterate again ...                 
    }

    return( T );

};

void PengRobinsonModel::calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {

    double P_norm   = critical_pressure;	/// Set pressure normalization factor
    double T_norm   = critical_temperature;	/// Set temperature normalization factor
    double nr_f     = -1.0;			/// Newton-Raphson residual value
    int nr_num_iter = 0;			/// Number of iterations required to obtain the solution

    /// Calculate P & T from rho & e by means of a Newton-Raphson solver
    nr_PT_unknowns[0] = P/P_norm;										/// Initialize unknown with previous pressure (normalized)
    nr_PT_unknowns[1] = T/T_norm;										/// Initialize unknown with previous temperature (normalized)
    nr_PT_solver->setExternalParameters( rho, e, P_norm, T_norm );						/// Set parameters of the solver
    nr_PT_solver->solve( nr_f, nr_PT_unknowns, max_nr_iter, nr_num_iter, nr_sing, nr_relative_tolerance );	/// Newton-Raphson solver 
    P = nr_PT_unknowns[0]*P_norm;										/// Update P & T from Newton-Raphson solver (unnormalized)
    T = nr_PT_unknowns[1]*T_norm;

};

void PengRobinsonModel::calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {

    /// Auxiliar parameters
    double eos_a  = this->calculate_eos_a( T );
    double eos_en = P*eos_b + R_universal*T;
    double a      = P;
    double b      = P*2.0*eos_b - eos_en;
    double c      = P*( -1.0 )*eos_b*eos_b - eos_en*2.0*eos_b + eos_a;
    double d      = ( -1.0 )*eos_b*( eos_a + eos_en*( -1.0 )*eos_b );

    /// Cubic solve to calculate rho
    complex<double> v_1, v_2, v_3;
    this->calculateRootsCubicPolynomial( v_1, v_2, v_3, a, b, c, d );
    rho = molecular_weight/real( v_1 );		/// First root is always real

    /// Calculate e
    double bar_v = molecular_weight/rho;
    e            = ( 1.0/molecular_weight )*this->calculateMolarInternalEnergyFromPressureTemperatureMolarVolume( P, T, bar_v );

};

void PengRobinsonModel::calculateSpecificHeatCapacities(double &c_v, double &c_p, const double &P, const double &T, const double &rho) {

    double bar_v       = molecular_weight/rho;
    double std_bar_c_p = this->calculateMolarStdCpFromNASApolynomials( T );
    double std_bar_c_v = std_bar_c_p - R_universal;

    c_v = ( 1.0/molecular_weight )*( std_bar_c_v + this->calculateDepartureFunctionMolarCv( P, T, bar_v ) );
    c_p = ( 1.0/molecular_weight )*( std_bar_c_p + this->calculateDepartureFunctionMolarCp( P, T, bar_v ) );

};

double PengRobinsonModel::calculateHeatCapacitiesRatio(const double &P, const double &rho) {

    double bar_v       = molecular_weight/rho;
    double T           = this->calculateTemperatureFromPressureMolarVolume( P, bar_v );
    double std_bar_c_p = this->calculateMolarStdCpFromNASApolynomials( T );
    double std_bar_c_v = std_bar_c_p - R_universal;

    double c_v = ( 1.0/molecular_weight )*( std_bar_c_v + this->calculateDepartureFunctionMolarCv( P, T, bar_v ) );
    double c_p = ( 1.0/molecular_weight )*( std_bar_c_p + this->calculateDepartureFunctionMolarCp( P, T, bar_v ) );

    double gamma = c_p/c_v;

    return( gamma );

};

double PengRobinsonModel::calculateSoundSpeed(const double &P, const double &T, const double &rho) {

    double bar_v = molecular_weight/rho;

    double sos = sqrt( 1.0/( rho*this->calculateIsentropicCompressibility( P, T, bar_v ) ) );

    return( sos );

};

double PengRobinsonModel::calculatePressureFromTemperatureDensity(const double &T, const double &rho) {

    double bar_v = molecular_weight/rho;

    double P = R_universal*T/( bar_v - eos_b ) - this->calculate_eos_a( T )/( bar_v*bar_v + 2.0*eos_b*bar_v - eos_b*eos_b );

    return( P );

};

double PengRobinsonModel::calculateMolarInternalEnergyFromPressureTemperatureMolarVolume(const double &P, const double &T, const double &bar_v) {

    double bar_e = this->calculateMolarStdEnthalpyFromNASApolynomials( T ) + this->calculateDepartureFunctionMolarEnthalpy( P, T, bar_v ) - P*bar_v;

    return( bar_e );

};

double PengRobinsonModel::calculate_eos_a(const double &T) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a = 0.457236*( pow( R_universal*critical_temperature, 2.0 )/critical_pressure )*pow( 1.0 + eos_kappa*( 1.0 - sqrt( T/critical_temperature ) ), 2.0 );

    return( eos_a );

};

double PengRobinsonModel::calculate_eos_a_first_derivative(const double &T) {
   
    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.
 
    double eos_a_first_derivative = eos_kappa*eos_ac*( ( eos_kappa/critical_temperature ) - ( ( 1.0 + eos_kappa )/sqrt( T*critical_temperature ) ) );

    return( eos_a_first_derivative );

};

double PengRobinsonModel::calculate_eos_a_second_derivative(const double &T) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a_second_derivative = ( eos_kappa*eos_ac*( 1.0 + eos_kappa ) )/( 2.0*sqrt( pow( T, 3.0 )*critical_temperature ) );

    return( eos_a_second_derivative );

};

double PengRobinsonModel::calculate_Z(const double &P, const double &T, const double &bar_v) {

    double Z = ( P*bar_v )/( R_universal*T );

    return( Z );

};

double PengRobinsonModel::calculate_A(const double &P, const double &T) {

    double eos_a = calculate_eos_a( T );

    double A = ( eos_a*P )/pow( R_universal*T, 2.0 );

    return( A );

};

double PengRobinsonModel::calculate_B(const double &P, const double &T) {

    double B = ( eos_b*P )/( R_universal*T );

    return( B );

};

double PengRobinsonModel::calculate_M(const double &Z, const double &B) {

    double M = ( Z*Z + 2.0*B*Z - B*B )/( Z - B );

    return( M );

};

double PengRobinsonModel::calculate_N(const double &eos_a_first_derivative, const double &B) {

    double N = eos_a_first_derivative*( B/( eos_b*R_universal ) );

    return( N );

};

double PengRobinsonModel::calculateMolarStdCpFromNASApolynomials(const double &T) {

    double std_bar_c_p = 0.0;

    if( ( T >= 200.0 ) && ( T < 1000.0 ) ) {
        std_bar_c_p = R_universal*( NASA_coefficients[7] + NASA_coefficients[8]*T + NASA_coefficients[9]*pow( T, 2.0 ) + NASA_coefficients[10]*pow( T, 3.0 ) + NASA_coefficients[11]*pow( T, 4.0 ) );
    } else if( ( T >= 1000.0 ) && ( T < 6000.0 ) ) {
        std_bar_c_p = R_universal*( NASA_coefficients[0] + NASA_coefficients[1]*T + NASA_coefficients[2]*pow( T, 2.0 ) + NASA_coefficients[3]*pow( T, 3.0 ) + NASA_coefficients[4]*pow( T, 4.0 ) );
    } else if ( T < 200 ) {
        // Assume constant temperature below T = 200 K	    
        double T_min = 200.0;	    
        
        std_bar_c_p = R_universal*( NASA_coefficients[7] + NASA_coefficients[8]*T_min + NASA_coefficients[9]*pow( T_min, 2.0 ) + NASA_coefficients[10]*pow( T_min, 3.0 ) + NASA_coefficients[11]*pow( T_min, 4.0 ) );
    } else {
        cout << endl << "NASA 7-coefficient polynomials for std bar c_p. T = " << T << " is above 6000 K." << endl << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    return( std_bar_c_p );

};

double PengRobinsonModel::calculateMolarStdEnthalpyFromNASApolynomials(const double &T) {

    double std_bar_h = 0.0;

    if( (T >= 200.0 ) && ( T < 1000.0 ) ) {
	//std_bar_h = R_universal*T*( NASA_coefficients[7] + NASA_coefficients[8]*T/2.0 + NASA_coefficients[9]*pow( T, 2.0 )/3.0 + NASA_coefficients[10]*pow( T, 3.0 )/4.0 + NASA_coefficients[11]*pow( T, 4.0 )/5.0 + NASA_coefficients[12]/T) - R_universal*NASA_coefficients[14];
	std_bar_h = R_universal*T*( NASA_coefficients[7] + NASA_coefficients[8]*T/2.0 + NASA_coefficients[9]*pow( T, 2.0 )/3.0 + NASA_coefficients[10]*pow( T, 3.0 )/4.0 + NASA_coefficients[11]*pow( T, 4.0 )/5.0 + NASA_coefficients[12]/T );
    } else if( ( T >= 1000.0 ) && ( T < 6000.0 ) ) {
	//std_bar_h = R_universal*T*( NASA_coefficients[0] + NASA_coefficients[1]*T/2.0 + NASA_coefficients[2]*pow( T, 2.0 )/3.0 + NASA_coefficients[3]*pow( T, 3.0 )/4.0 + NASA_coefficients[4]*pow( T, 4.0 )/5.0 + NASA_coefficients[5]/T ) - R_universal*NASA_coefficients[14];
	std_bar_h = R_universal*T*( NASA_coefficients[0] + NASA_coefficients[1]*T/2.0 + NASA_coefficients[2]*pow( T, 2.0 )/3.0 + NASA_coefficients[3]*pow( T, 3.0 )/4.0 + NASA_coefficients[4]*pow( T, 4.0 )/5.0 + NASA_coefficients[5]/T );
    } else if( T < 200.0 ) {
	// Assume linear interpolation from T = 200 K 
        double T_min = 200.0;
	    
	//double std_bar_h_min   = R_universal*T_min*( NASA_coefficients[7] + NASA_coefficients[8]*T_min/2.0 + NASA_coefficients[9]*pow( T_min, 2.0 )/3.0 + NASA_coefficients[10]*pow( T_min, 3.0 )/4.0 + NASA_coefficients[11]*pow( T_min, 4.0 )/5.0 + NASA_coefficients[12]/T_min ) - R_universal*NASA_coefficients[14];
	double std_bar_h_min   = R_universal*T_min*( NASA_coefficients[7] + NASA_coefficients[8]*T_min/2.0 + NASA_coefficients[9]*pow( T_min, 2.0 )/3.0 + NASA_coefficients[10]*pow( T_min, 3.0 )/4.0 + NASA_coefficients[11]*pow( T_min, 4.0 )/5.0 + NASA_coefficients[12]/T_min );
	double std_bar_h_slope = R_universal*( NASA_coefficients[7] + NASA_coefficients[8]*T_min + NASA_coefficients[9]*pow( T_min, 2.0 ) + NASA_coefficients[10]*pow( T_min, 3.0 ) + NASA_coefficients[11]*pow( T_min, 4.0 ) );

	std_bar_h = std_bar_h_min + std_bar_h_slope*( T - T_min );
    } else {
	cout << endl << "NASA 7-coefficient polynomials for std bar h. T = " << T << " is above 6000 K." << endl << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    return( std_bar_h );

};

double PengRobinsonModel::calculateDepartureFunctionMolarCp(const double &P, const double &T, const double &bar_v) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a_first_derivative  = this->calculate_eos_a_first_derivative( T );
    double eos_a_second_derivative = this->calculate_eos_a_second_derivative( T );
    double Z                       = this->calculate_Z( P, T, bar_v );
    double A                       = this->calculate_A( P, T );
    double B                       = this->calculate_B( P, T );
    double M                       = this->calculate_M( Z, B );
    double N                       = this->calculate_N( eos_a_first_derivative, B );

    double Delta_bar_c_p = ( ( R_universal*pow( M - N ,2.0 ) )/( pow( M, 2.0 ) - 2.0*A*( Z + B ) ) ) - ( ( T*eos_a_second_derivative )/( 2.0*sqrt(2.0)*eos_b ) )*log( ( Z + ( 1.0 - sqrt( 2.0 ) )*B )/( Z + ( 1.0 + sqrt( 2.0 ) )*B ) ) - R_universal; 

    return( Delta_bar_c_p );
  
};

double PengRobinsonModel::calculateDepartureFunctionMolarCv(const double &P, const double &T, const double &bar_v) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a_second_derivative = this->calculate_eos_a_second_derivative( T );
    double Z                       = this->calculate_Z( P, T, bar_v );
    double B                       = this->calculate_B( P, T );

    double Delta_bar_c_v = ( -1.0 )*( ( T*eos_a_second_derivative )/( 2.0*sqrt( 2.0 )*eos_b ) )*log( ( Z + ( 1.0 - sqrt( 2.0 ) )*B )/( Z + ( 1.0 + sqrt( 2.0 ) )*B ) );

    return( Delta_bar_c_v );
  
  };

double PengRobinsonModel::calculateDepartureFunctionMolarEnthalpy(const double &P, const double &T, const double &bar_v) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a                  = this->calculate_eos_a( T );
    double eos_a_first_derivative = this->calculate_eos_a_first_derivative( T );
    double Z                      = this->calculate_Z( P, T, bar_v );
    double B                      = this->calculate_B( P, T );

    double Delta_bar_h = R_universal*T*( Z - 1.0 ) + ( ( eos_a - eos_a_first_derivative*T )/( 2.0*sqrt( 2.0 )*eos_b ) )*log( ( Z + ( 1.0 - sqrt( 2.0 ) )*B )/( Z + ( 1.0 + sqrt( 2.0 ) )*B ) );

    return( Delta_bar_h );
  
};

double PengRobinsonModel::calculateTemperatureFromPressureMolarVolume(const double &P, const double &bar_v) {

    /// Numerical Recipes in C++, Second Edition.
    /// W.H. Press, S.A. Teulosky, W.T. Vetterling, B.P. Flannery.
    /// 5.1 Series and Their Convergence: Aitken’s delta-squared process.

    // Initial temperature guess using ideal-gas model
    double T = P*bar_v/R_universal;

    /// Aitken’s delta-squared process:
    double x_0 = T, x_1, x_2, denominator;
    for(int iter = 0; iter < max_aitken_iter; iter++) { 
        x_1 = ( (bar_v - eos_b )/R_universal )*( P + ( this->calculate_eos_a( x_0 )/( pow( bar_v, 2.0 ) + 2.0*eos_b*bar_v - pow( eos_b, 2.0 ) ) ) );
        x_2 = ( (bar_v - eos_b )/R_universal )*( P + ( this->calculate_eos_a( x_1 )/( pow( bar_v, 2.0 ) + 2.0*eos_b*bar_v - pow( eos_b, 2.0 ) ) ) );

        denominator = x_2 - 2.0*x_1 + x_0;
        T = x_2 - ( pow( x_2 - x_1, 2.0 ) )/denominator;
    
        if( abs( ( T - x_2 )/ T ) < aitken_relative_tolerance ) break;	/// If the result is within tolerance, leave the loop!
        x_0 = T;							/// Otherwise, update x_0 to iterate again ...                 
    }

    return( T );

};

double PengRobinsonModel::calculateDPDTConstantMolarVolume(const double &T, const double &bar_v) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a_first_derivative = this->calculate_eos_a_first_derivative( T );

    double dP_dT_const_v = ( R_universal/( bar_v - eos_b ) ) - ( eos_a_first_derivative/( bar_v*bar_v + 2.0*bar_v*eos_b - eos_b*eos_b ) );

    return( dP_dT_const_v );
  
};

double PengRobinsonModel::calculateDPDvConstantTemperature(const double &T, const double &bar_v) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double eos_a = this->calculate_eos_a( T );
    
    double dP_dv_const_T = ( -1.0 )*( ( R_universal*T )/pow( bar_v - eos_b, 2.0 ) ) + ( eos_a*( 2.0*bar_v + 2.0*eos_b ) )/pow( pow( bar_v, 2.0 ) + 2.0*bar_v*eos_b - pow( eos_b, 2.0 ), 2.0 );

    return( dP_dv_const_T );
  
};  

double PengRobinsonModel::calculateVolumeExpansivity(const double &T, const double &bar_v) {

    double dP_dT_const_v = this->calculateDPDTConstantMolarVolume( T, bar_v );
    double dP_dv_const_T = this->calculateDPDvConstantTemperature( T, bar_v );

    double expansivity = ( -1.0 )*( dP_dT_const_v/( bar_v*dP_dv_const_T ) );

    return( expansivity );
  
};

double PengRobinsonModel::calculateIsothermalCompressibility(const double &T, const double &bar_v) {

    double dP_dv_const_T = this->calculateDPDvConstantTemperature( T, bar_v );

    double isothermal_compressibility = ( -1.0 )/( bar_v*dP_dv_const_T );

    return( isothermal_compressibility );
  
};  

double PengRobinsonModel::calculateIsentropicCompressibility(const double &P, const double &T, const double &bar_v) {

    double isothermal_compressibility = this->calculateIsothermalCompressibility( T, bar_v );
    double expansivity                = this->calculateVolumeExpansivity( T, bar_v );
    double bar_c_p                    = this->calculateMolarStdCpFromNASApolynomials( T ) + this->calculateDepartureFunctionMolarCp( P, T, bar_v );
      
    double isentropic_compressibility = ( isothermal_compressibility - ( ( bar_v*T*pow( expansivity, 2.0 ) )/bar_c_p ) );
    
    return( isentropic_compressibility );
  
};
        
void PengRobinsonModel::calculateRootsCubicPolynomial(complex<double> &root_1, complex<double> &root_2, complex<double> &root_3, double &a, double &b, double &c, double &d) {

    if( a == 0 ) {	/// End if a == 0
        cout << "The coefficient of the cube of x is 0. Please use the utility for a SECOND degree quadratic. No further action taken." << endl;
        return;
    }
    if( d == 0 ) {	/// End if d == 0
        cout << "One root is 0. Now divide through by x and use the utility for a SECOND degree quadratic to solve the resulting equation for the other two roots. No further action taken." << endl;
        return;
    }

    b /= a;
    c /= a;
    d /= a;

    double disc, q, r, dum1, s, t, term1, r13;
    q = (3.0*c - (b*b))/9.0;
    r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
    r /= 54.0;
    disc = q*q*q + r*r;

    root_1.imag( 0 );	/// The first root is always real
    term1 = (b/3.0);

    if( disc > 0 ) {	/// One root real, two are complex
        s = r + sqrt(disc);
        s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
        t = r - sqrt(disc);
        t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
        root_1.real( -term1 + s + t );
        term1 += (s + t)/2.0;
        root_2.real( -term1 );
        root_3.real( -term1 );
        term1 = sqrt(3.0)*(-t + s)/2;
        root_2.imag( term1 );
        root_3.imag( -term1 );
        return;
    }	/// End if (disc > 0)

    cout << endl;
    cout << "The PengRobinsonModel::calculateRootsCubicPolynomial has found more than one real root." << endl;
    cout << "Vapor-liquid equilibrium conditions must be solved." << endl;
    cout << "Another option is to avoid calculating rho from P and T." << endl;
    cout << endl;

    /// The remaining options are all real
    root_2.imag( 0 );
    root_3.imag( 0 );
    if(disc == 0) {	/// All roots real, at least two are equal
        r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
        root_1.real( -term1 + 2.0*r13 );
        root_2.real( -(r13 + term1) );
        root_3.real( -(r13 + term1) );
        return;
    }	/// End if (disc == 0)

    /// Only option left is that all roots are real and unequal (to get here, q < 0)
    q = -q;
    dum1 = q*q*q;
    dum1 = acos(r/sqrt(dum1));
    r13 = 2.0*sqrt(q);
    root_1.real( -term1 + r13*cos(dum1/3.0) );
    root_2.real( -term1 + r13*cos((dum1 + 2.0*3.141592654)/3.0) );
    root_3.real( -term1 + r13*cos((dum1 + 4.0*3.141592654)/3.0) );

    return;

};

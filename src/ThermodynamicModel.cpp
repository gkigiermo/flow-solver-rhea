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

    /// Fluid properties
    const YAML::Node & fluid_properties = configuration["fluid_properties"];
    R_specific = fluid_properties["R_specific"].as<double>();
    gamma      = fluid_properties["gamma"].as<double>();

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

void IdealGasModel::calculateSpecificHeatCapacities(double &c_v, double &c_p) {

    /// Ideal-gas model:
    /// c_v = R_specific/(gamma - 1)
    /// c_p = c_v*gamma

    c_v = R_specific/( gamma - 1.0 );
    c_p = c_v*gamma;

};

double IdealGasModel::calculateHeatCapacitiesRatio() {

    /// Ideal-gas model:
    /// gamma = c_p/c_v

    return( gamma );

};

double IdealGasModel::calculateSoundSpeed(const double &rho, const double &P, const double &T) {

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

    /// Fluid properties
    const YAML::Node & fluid_properties = configuration["fluid_properties"];
    R_specific = fluid_properties["R_specific"].as<double>();
    gamma      = fluid_properties["gamma"].as<double>();
    P_inf      = fluid_properties["P_inf"].as<double>();
    e_0        = fluid_properties["e_0"].as<double>();
    c_v        = fluid_properties["c_v"].as<double>();

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

void StiffenedGasModel::calculateSpecificHeatCapacities(double &c_v_, double &c_p) {

    /// Stiffened-gas model:
    /// c_p = c_v*gamma

    c_v_ = c_v;
    c_p  = c_v*gamma;

};

double StiffenedGasModel::calculateHeatCapacitiesRatio() {

    /// Stiffened-gas model:
    /// gamma = c_p/c_v

    return( gamma );

};

double StiffenedGasModel::calculateSoundSpeed(const double &rho, const double &P, const double &T) {

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

};

PengRobinsonModel::~PengRobinsonModel() {};

void PengRobinsonModel::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid properties
    const YAML::Node & fluid_properties = configuration["fluid_properties"];
    R_specific            = fluid_properties["R_specific"].as<double>();
    molecular_weight      = fluid_properties["molecular_weight"].as<double>();
    acentric_factor       = fluid_properties["acentric_factor"].as<double>();
    critical_temperature  = fluid_properties["critical_temperature"].as<double>();
    critical_pressure     = fluid_properties["critical_pressure"].as<double>();
    critical_molar_volume = fluid_properties["critical_molar_volume"].as<double>();
    NASA_coefficients[0]  = fluid_properties["NASA_coefficients"][0].as<double>();
    NASA_coefficients[1]  = fluid_properties["NASA_coefficients"][1].as<double>();
    NASA_coefficients[2]  = fluid_properties["NASA_coefficients"][2].as<double>();
    NASA_coefficients[3]  = fluid_properties["NASA_coefficients"][3].as<double>();
    NASA_coefficients[4]  = fluid_properties["NASA_coefficients"][4].as<double>();
    NASA_coefficients[5]  = fluid_properties["NASA_coefficients"][5].as<double>();
    NASA_coefficients[6]  = fluid_properties["NASA_coefficients"][6].as<double>();
    NASA_coefficients[7]  = fluid_properties["NASA_coefficients"][7].as<double>();
    NASA_coefficients[8]  = fluid_properties["NASA_coefficients"][8].as<double>();
    NASA_coefficients[9]  = fluid_properties["NASA_coefficients"][9].as<double>();
    NASA_coefficients[10] = fluid_properties["NASA_coefficients"][10].as<double>();
    NASA_coefficients[11] = fluid_properties["NASA_coefficients"][11].as<double>();
    NASA_coefficients[12] = fluid_properties["NASA_coefficients"][12].as<double>();
    NASA_coefficients[13] = fluid_properties["NASA_coefficients"][13].as<double>();
    NASA_coefficients[14] = fluid_properties["NASA_coefficients"][14].as<double>();

};

void PengRobinsonModel::calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {

    P = 1.0/0.0; 
    T = 1.0/0.0;

};

void PengRobinsonModel::calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {

    e   = 1.0/0.0;
    rho = 1.0/0.0;

};

void PengRobinsonModel::calculateSpecificHeatCapacities(double &c_v, double &c_p) {

    c_v = 1.0/0.0;
    c_p = 1.0/0.0;

};

double PengRobinsonModel::calculateHeatCapacitiesRatio() {

    double gamma = 1.0/0.0;

    return( gamma );

};

double PengRobinsonModel::calculateSoundSpeed(const double &rho, const double &P, const double &T) {

    double sos = 1.0/0.0;

    return( sos );

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

double PengRobinsonModel::calculate_Z(const double &P, const double &T, const double &v) {

    double Z = ( P*v )/( R_universal*T );

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

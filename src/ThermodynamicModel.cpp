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

};

PengRobinsonModel::~PengRobinsonModel() {};

void PengRobinsonModel::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid properties
    const YAML::Node & fluid_properties = configuration["fluid_properties"];
    R_specific = fluid_properties["R_specific"].as<double>();

};

void PengRobinsonModel::calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.
    
    P = 1.0/0.0; 
    T = 1.0/0.0;

};

void PengRobinsonModel::calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    e   = 1.0/0.0;
    rho = 1.0/0.0;

};

void PengRobinsonModel::calculateSpecificHeatCapacities(double &c_v, double &c_p) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    c_v = 1.0/0.0;
    c_p = 1.0/0.0;

};

double PengRobinsonModel::calculateSoundSpeed(const double &rho, const double &P, const double &T) {

    /// Peng-Robinson model:
    /// D. Y. Peng, D. B. Robinson.
    /// A new two-constant equation of state.
    /// Industrial and Engineering Chemistry: Fundamentals, 15, 59-64, 1976.

    double sos = 1.0/0.0;

    return( sos );

};

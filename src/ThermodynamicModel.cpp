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

void IdealGasModel::calculatePointPressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {

    /// Ideal-gas model:
    /// P = e*rho*(gamma - 1) is pressure
    /// T = e/c_v is temperature

    double c_v, c_p;
    this->calculatePointSpecificHeatCapacities( c_v, c_p );

    P = e*rho*( ( c_p/c_v ) - 1.0 ); 
    T = e/c_v;

};

void IdealGasModel::calculatePointDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {

    /// Ideal-gas model:
    /// rho = P/( e*( gamma - 1.0 ) ) is density
    /// e = c_v*T is specific internal energy

    double c_v, c_p;
    this->calculatePointSpecificHeatCapacities( c_v, c_p );

    e   = c_v*T;
    rho = P/( e*( ( c_p/c_v ) - 1.0 ) );

};

void IdealGasModel::calculatePointSpecificHeatCapacities(double &c_v, double &c_p) {

    /// Ideal-gas model:
    /// c_v = R_specific/(gamma - 1)
    /// c_p = c_v*gamma

    c_v = R_specific/( gamma - 1.0 );
    c_p = c_v*gamma;

};

double IdealGasModel::calculatePointSoundSpeed(const double &rho, const double &P, const double &T) {

    /// Ideal-gas model:
    /// sos = sqrt(gamma*P/rho) is speed of sound

    double c_v, c_p;
    this->calculatePointSpecificHeatCapacities( c_v, c_p );

    double sos = sqrt( ( c_p/c_v )*P/rho );

    return( sos );

};

#include "TransportCoefficients.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////


////////// BaseTransportCoefficients CLASS //////////

BaseTransportCoefficients::BaseTransportCoefficients() {};
        
BaseTransportCoefficients::BaseTransportCoefficients(const string configuration_file) : configuration_file(configuration_file) {};

BaseTransportCoefficients::~BaseTransportCoefficients() {};


////////// ConstantTransportCoefficients CLASS //////////

ConstantTransportCoefficients::ConstantTransportCoefficients() : BaseTransportCoefficients() {};
        
ConstantTransportCoefficients::ConstantTransportCoefficients(const string configuration_file) : BaseTransportCoefficients(configuration_file) {

    /// Read configuration (input) file
    this->readConfigurationFile();

};

ConstantTransportCoefficients::~ConstantTransportCoefficients() {};

void ConstantTransportCoefficients::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid properties
    const YAML::Node & fluid_properties = configuration["fluid_properties"];
    mu    = fluid_properties["mu"].as<double>();
    kappa = fluid_properties["kappa"].as<double>();

};
        
double ConstantTransportCoefficients::calculateDynamicViscosity(const double &P, const double &T) {

  return( mu );

};
        
double ConstantTransportCoefficients::calculateThermalConductivity(const double &P, const double &T) {

  return( kappa );

};

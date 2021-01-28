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

    /// Fluid & flow properties
    const YAML::Node & fluid_flow_properties = configuration["fluid_flow_properties"];
    if( fluid_flow_properties["substance_name"] ) {

        /// Create YAML object
	string substance_name          = fluid_flow_properties["substance_name"].as<string>();
	string substances_library_file = fluid_flow_properties["substances_library_file"].as<string>();
        YAML::Node substances_library  = YAML::LoadFile( substances_library_file );
        YAML::Node substance;

	if( substance_name == "NITROGEN" ) {
            substance = substances_library["NITROGEN"];
	} else {
            cout << "Substance not available!" << endl;
            MPI_Abort( MPI_COMM_WORLD, 1 );
	}

        mu    = substance["mu"].as<double>();
        kappa = substance["kappa"].as<double>();

    } else {	

        mu    = fluid_flow_properties["mu"].as<double>();
        kappa = fluid_flow_properties["kappa"].as<double>();

    }

};
        
double ConstantTransportCoefficients::calculateDynamicViscosity(const double &P, const double &T, const double &rho) {

    return( mu );

};
        
double ConstantTransportCoefficients::calculateThermalConductivity(const double &P, const double &T, const double &rho) {

    return( kappa );

};


////////// HighPressureTransportCoefficients CLASS //////////

HighPressureTransportCoefficients::HighPressureTransportCoefficients() : BaseTransportCoefficients() {};
        
HighPressureTransportCoefficients::HighPressureTransportCoefficients(const string configuration_file) : BaseTransportCoefficients(configuration_file) {

    /// Read configuration (input) file
    this->readConfigurationFile();

    // Adimensional dipole moment -- Poling et al. The properties of gases and liquids. McGraw-Hill, 2001.
    adimensional_dipole_moment = 131.3*( dipole_moment/sqrt( ( 1.0e6*critical_molar_volume )*critical_temperature ) );
      
    // Viscosity mu -- Poling et al. The properties of gases and liquids. McGraw-Hill, 2001. (9.40, Table 9-6)
    const double a1_mu  = 6.324,    b1_mu  = 50.412,    c1_mu  = -51.680,   d1_mu  = 1189.0;
    const double a2_mu  = 1.210e-3, b2_mu  = -1.154e-3, c2_mu  = -6.257e-3, d2_mu  = 0.03728;
    const double a3_mu  = 5.283,    b3_mu  = 254.209,   c3_mu  = -168.48,   d3_mu  = 3898.0;
    const double a4_mu  = 6.623,    b4_mu  = 38.096,    c4_mu  = -8.464,    d4_mu  = 31.42;
    const double a5_mu  = 19.745,   b5_mu  = 7.630,     c5_mu  = -14.354,   d5_mu  = 31.53;
    const double a6_mu  = -1.900,   b6_mu  = -12.537,   c6_mu  = 4.985,     d6_mu  = -18.15;
    const double a7_mu  = 24.275,   b7_mu  = 3.450,     c7_mu  = -11.291,   d7_mu  = 69.35;
    const double a8_mu  = 0.7972,   b8_mu  = 1.117,     c8_mu  = 0.01235,   d8_mu  = -4.117;
    const double a9_mu  = -0.2382,  b9_mu  = 0.06770,   c9_mu  = -0.8163,   d9_mu  = 4.025;
    const double a10_mu = 0.06863,  b10_mu = 0.3479,    c10_mu = 0.5926,    d10_mu = -0.727;

    E1_mu  = a1_mu  + b1_mu*acentric_factor  + c1_mu*pow( adimensional_dipole_moment, 4.0 )  + d1_mu*association_factor;
    E2_mu  = a2_mu  + b2_mu*acentric_factor  + c2_mu*pow( adimensional_dipole_moment, 4.0 )  + d2_mu*association_factor;
    E3_mu  = a3_mu  + b3_mu*acentric_factor  + c3_mu*pow( adimensional_dipole_moment, 4.0 )  + d3_mu*association_factor;
    E4_mu  = a4_mu  + b4_mu*acentric_factor  + c4_mu*pow( adimensional_dipole_moment, 4.0 )  + d4_mu*association_factor;
    E5_mu  = a5_mu  + b5_mu*acentric_factor  + c5_mu*pow( adimensional_dipole_moment, 4.0 )  + d5_mu*association_factor;
    E6_mu  = a6_mu  + b6_mu*acentric_factor  + c6_mu*pow( adimensional_dipole_moment, 4.0 )  + d6_mu*association_factor;
    E7_mu  = a7_mu  + b7_mu*acentric_factor  + c7_mu*pow( adimensional_dipole_moment, 4.0 )  + d7_mu*association_factor;
    E8_mu  = a8_mu  + b8_mu*acentric_factor  + c8_mu*pow( adimensional_dipole_moment, 4.0 )  + d8_mu*association_factor;
    E9_mu  = a9_mu  + b9_mu*acentric_factor  + c9_mu*pow( adimensional_dipole_moment, 4.0 )  + d9_mu*association_factor;
    E10_mu = a10_mu + b10_mu*acentric_factor + c10_mu*pow( adimensional_dipole_moment, 4.0 ) + d10_mu*association_factor;

    // Thermal conductivity k -- Poling et al. The properties of gases and liquids. McGraw-Hill, 2001. (10.23, Table 10-3)
    const double a1_k = 2.4166*1.00, b1_k = 7.4824*0.100, c1_k = -9.1858*0.10, d1_k = 1.2172*100.0;
    const double a2_k = -5.0924*0.1, b2_k = -1.5094*1.00, c2_k = -4.9991*10.0, d2_k = 6.9983*10.00;
    const double a3_k = 6.6107*1.00, b3_k = 5.6207*1.000, c3_k = 6.4760*10.00, d3_k = 2.7039*10.00;
    const double a4_k = 1.4543*10.0, b4_k = -8.9139*1.00, c4_k = -5.6379*1.00, d4_k = 7.4344*10.00;
    const double a5_k = 7.9274*0.10, b5_k = 8.2019*0.100, c5_k = -6.9369*0.10, d5_k = 6.3173*1.000;
    const double a6_k = -5.8634*1.0, b6_k = 1.2801*10.00, c6_k = 9.5893*1.000, d6_k = 6.5529*10.00;
    const double a7_k = 9.1089*10.0, b7_k = 1.2811*100.0, c7_k = -5.4217*10.0, d7_k = 5.2381*100.0;

    B1_k = a1_k + b1_k*acentric_factor + c1_k*pow( adimensional_dipole_moment, 4.0 ) + d1_k*association_factor;
    B2_k = a2_k + b2_k*acentric_factor + c2_k*pow( adimensional_dipole_moment, 4.0 ) + d2_k*association_factor;
    B3_k = a3_k + b3_k*acentric_factor + c3_k*pow( adimensional_dipole_moment, 4.0 ) + d3_k*association_factor;
    B4_k = a4_k + b4_k*acentric_factor + c4_k*pow( adimensional_dipole_moment, 4.0 ) + d4_k*association_factor;
    B5_k = a5_k + b5_k*acentric_factor + c5_k*pow( adimensional_dipole_moment, 4.0 ) + d5_k*association_factor;
    B6_k = a6_k + b6_k*acentric_factor + c6_k*pow( adimensional_dipole_moment, 4.0 ) + d6_k*association_factor;
    B7_k = a7_k + b7_k*acentric_factor + c7_k*pow( adimensional_dipole_moment, 4.0 ) + d7_k*association_factor;

};

HighPressureTransportCoefficients::~HighPressureTransportCoefficients() {};

void HighPressureTransportCoefficients::readConfigurationFile() {

    /// Create YAML object
    YAML::Node configuration = YAML::LoadFile(configuration_file);

    /// Fluid & flow properties
    const YAML::Node & fluid_flow_properties = configuration["fluid_flow_properties"];
    if( fluid_flow_properties["substance_name"] ) {

        /// Create YAML object
	string substance_name          = fluid_flow_properties["substance_name"].as<string>();
	string substances_library_file = fluid_flow_properties["substances_library_file"].as<string>();
        YAML::Node substances_library  = YAML::LoadFile( substances_library_file );
        YAML::Node substance;

	if( substance_name == "NITROGEN" ) {
            substance = substances_library["NITROGEN"];
	} else {
            cout << "Substance not available!" << endl;
            MPI_Abort( MPI_COMM_WORLD, 1 );
	}

        molecular_weight      = substance["molecular_weight"].as<double>();
        critical_temperature  = substance["critical_temperature"].as<double>();
        critical_molar_volume = substance["critical_molar_volume"].as<double>();
        acentric_factor       = substance["acentric_factor"].as<double>();
        dipole_moment         = substance["dipole_moment"].as<double>();
        association_factor    = substance["association_factor"].as<double>();
        NASA_coefficients[0]  = substance["NASA_coefficients"][0].as<double>();
        NASA_coefficients[1]  = substance["NASA_coefficients"][1].as<double>();
        NASA_coefficients[2]  = substance["NASA_coefficients"][2].as<double>();
        NASA_coefficients[3]  = substance["NASA_coefficients"][3].as<double>();
        NASA_coefficients[4]  = substance["NASA_coefficients"][4].as<double>();
        NASA_coefficients[5]  = substance["NASA_coefficients"][5].as<double>();
        NASA_coefficients[6]  = substance["NASA_coefficients"][6].as<double>();
        NASA_coefficients[7]  = substance["NASA_coefficients"][7].as<double>();
        NASA_coefficients[8]  = substance["NASA_coefficients"][8].as<double>();
        NASA_coefficients[9]  = substance["NASA_coefficients"][9].as<double>();
        NASA_coefficients[10] = substance["NASA_coefficients"][10].as<double>();
        NASA_coefficients[11] = substance["NASA_coefficients"][11].as<double>();
        NASA_coefficients[12] = substance["NASA_coefficients"][12].as<double>();
        NASA_coefficients[13] = substance["NASA_coefficients"][13].as<double>();
        NASA_coefficients[14] = substance["NASA_coefficients"][14].as<double>();

    } else {

        molecular_weight      = fluid_flow_properties["molecular_weight"].as<double>();
        critical_temperature  = fluid_flow_properties["critical_temperature"].as<double>();
        critical_molar_volume = fluid_flow_properties["critical_molar_volume"].as<double>();
        acentric_factor       = fluid_flow_properties["acentric_factor"].as<double>();
        dipole_moment         = fluid_flow_properties["dipole_moment"].as<double>();
        association_factor    = fluid_flow_properties["association_factor"].as<double>();
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

    }

};
        
double HighPressureTransportCoefficients::calculateDynamicViscosity(const double &P, const double &T, const double &rho) {

    /// T. H. Chung, L. L. Lee, K. E. Starling.
    /// Applications of kinetic gas theories and multiparameter correlation for prediction of dilute gas viscosity and thermal conductivity.
    /// Industrial & Engineering Chemistry Fundamentals, 23, 8-13, 1984.

    /// T. H. Chung, M. Ajlan, L. L. Lee, K. E. Starling.
    /// Generalized multiparameter correlation for nonpolar and polar fluid transport properties.
    /// Industrial & Engineering Chemistry Fundamentals, 27, 671-679, 1988.

    /// Auxiliar coefficients
    double v     = molecular_weight/rho;
    double Y     = critical_molar_volume/( 6.0*v );
    double T_ast = 1.2593*( T/critical_temperature );
    double Omega = 1.16145*pow( T_ast, -0.14874 ) + 0.52487*exp( -0.77320*T_ast ) + 2.16178*exp( -2.43787*T_ast );
    double G1    = ( 1.0 - 0.5*Y )/pow( 1.0 - Y, 3.0 );
    double Fc    = 1.0 - 0.2756*acentric_factor + 0.059035*pow( adimensional_dipole_moment, 4.0 ) + association_factor;

    /// Additional auxiliar coefficients
    double G2_mu      = ( E1_mu*( 1.0 - exp( -E4_mu*Y ) )/Y + E2_mu*G1*exp( E5_mu*Y ) + E3_mu*G1 )/( E1_mu*E4_mu + E2_mu + E3_mu );
    double mu_ast_ast = E7_mu*pow( Y, 2.0 )*G2_mu*exp( E8_mu + E9_mu/T_ast + E10_mu*pow( T_ast, -2.0 ) );
    double mu_ast     = ( sqrt( T_ast )*Fc/Omega )*( 1.0/G2_mu + E6_mu*Y ) + mu_ast_ast;

    /// Calculate viscosity
    double mu = 1.0e-7*mu_ast*( ( 36.344*sqrt( ( 1.0e3*molecular_weight )*critical_temperature ) )/pow( 1.0e6*critical_molar_volume, 2.0/3.0 ) );

    return( mu );

};
        
double HighPressureTransportCoefficients::calculateThermalConductivity(const double &P, const double &T, const double &rho) {

    /// T. H. Chung, L. L. Lee, K. E. Starling.
    /// Applications of kinetic gas theories and multiparameter correlation for prediction of dilute gas viscosity and thermal conductivity.
    /// Industrial & Engineering Chemistry Fundamentals, 23, 8-13, 1984.

    /// T. H. Chung, M. Ajlan, L. L. Lee, K. E. Starling.
    /// Generalized multiparameter correlation for nonpolar and polar fluid transport properties.
    /// Industrial & Engineering Chemistry Fundamentals, 27, 671-679, 1988.

    /// Auxiliar coefficients
    double v           = molecular_weight/rho;
    double Y           = critical_molar_volume/( 6.0*v );
    double T_ast       = 1.2593*( T/critical_temperature );
    double Omega       = 1.16145*pow( T_ast, -0.14874 ) + 0.52487*exp( -0.77320*T_ast ) + 2.16178*exp( -2.43787*T_ast );
    double G1          = ( 1.0 - 0.5*Y )/pow( 1.0 - Y, 3.0 );
    double Fc          = 1.0 - 0.2756*acentric_factor + 0.059035*pow( adimensional_dipole_moment, 4.0 ) + association_factor;
    double std_bar_c_p = this->calculateMolarStdCpFromNASApolynomials( T );

    /// Additional auxiliar coefficients
    double mu_0_k  = 40.785e-7*Fc*sqrt( 1.0e3*molecular_weight*T )/( pow( 1.0e6*critical_molar_volume, 2.0/3.0 )*Omega );
    double alpha_k = ( std_bar_c_p/R_universal - 1.0 ) - 1.5;
    double beta_k  = 0.7862 - 0.7109*acentric_factor + 1.3168*pow( acentric_factor, 2.0 );
    double gamma_k = 2.0 + 10.5*pow( T/critical_temperature, 2.0 );
    double Psi_k   = 1.0 + alpha_k*( ( 0.215 + 0.28288*alpha_k - 1.061*beta_k + 0.26665*gamma_k )/( 0.6366 + beta_k*gamma_k + 1.061*alpha_k*beta_k) );
    double q_k     = 0.003586*( sqrt( critical_temperature/molecular_weight )/pow( ( 1.0e6*critical_molar_volume ), 2.0/3.0 ) );
    double G3_k    = ( ( ( B1_k/Y )*( 1.0 - exp( (-1.0)*B4_k*Y ) ) ) + ( B2_k*G1*exp( B5_k*Y ) ) + ( B3_k*G1 ) )/( B1_k*B4_k + B2_k + B3_k );

    /// Calculate thermal conductivity
    double kappa = ( 31.2*mu_0_k*Psi_k/molecular_weight )*( 1.0/G3_k + B6_k*Y ) + q_k*B7_k*pow( Y, 2.0 )*sqrt( T/critical_temperature)*G3_k;

    return( kappa );

};

double HighPressureTransportCoefficients::calculateMolarStdCpFromNASApolynomials(const double &T) {

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

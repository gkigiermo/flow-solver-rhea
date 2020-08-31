#ifndef _TRANSPORT_COEFFICIENTS_
#define _TRANSPORT_COEFFICIENTS_

////////// INCLUDES //////////
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include "yaml-cpp/yaml.h"

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class BaseTransportCoefficients;					/// Base transport coefficients
class ConstantTransportCoefficients;					/// Constant transport coefficients

////////// FUNCTION DECLARATION //////////


////////// BaseTransportCoefficients CLASS //////////
class BaseTransportCoefficients {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        BaseTransportCoefficients();					/// Default constructor
        BaseTransportCoefficients(const string configuration_file);	/// Parametrized constructor
        virtual ~BaseTransportCoefficients();				/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getDynamicViscosity() { return( mu ); };
        inline double getThermalConductivity() { return( kappa ); };

	////////// SET FUNCTIONS //////////
        inline void setDynamicViscosity(double mu_) { mu = mu_; };
        inline void setThermalConductivity(double kappa_) { kappa = kappa_; };

	////////// METHODS //////////
       
        /// Read configuration (input) file written in YAML language
        virtual void readConfigurationFile() {};
 
        /// Calculate dynamic viscosity
        virtual double calculateDynamicViscosity(const double &P, const double &T) = 0;

        /// Calculate thermal conductivity
        virtual double calculateThermalConductivity(const double &P, const double &T) = 0;

    protected:

        ////////// PARAMETERS //////////

        /// Transport coefficients
        double mu;						/// Dynamic viscosity [Pa·s]
        double kappa;						/// Thermal conductivity [W/(m·k)]

        /// Model parameters
        string configuration_file;				/// Configuration file name (YAML language)

    private:

};


////////// ConstantTransportCoefficients CLASS //////////
class ConstantTransportCoefficients : public BaseTransportCoefficients {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        ConstantTransportCoefficients();					/// Default constructor
        ConstantTransportCoefficients(const string configuration_file);		/// Parametrized constructor
        virtual ~ConstantTransportCoefficients();				/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate dynamic viscosity
        double calculateDynamicViscosity(const double &P, const double &T);

        /// Calculate thermal conductivity
        double calculateThermalConductivity(const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

    private:

};

#endif /*_TRANSPORT_COEFFICIENTS_*/

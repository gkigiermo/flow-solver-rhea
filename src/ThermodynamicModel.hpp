#ifndef _THERMODYNAMIC_MODELS_
#define _THERMODYNAMIC_MODELS_

////////// INCLUDES //////////
#include <stdio.h>
#include <cmath>
#include <iostream>
#include "yaml-cpp/yaml.h"
#include "RootFindingMinimization.hpp"

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class BaseThermodynamicModel;				/// Base thermodynamic model
class IdealGasModel;					/// Ideal-gas thermodynamic model
class StiffenedGasModel;				/// Stiffened-gas thermodynamic model
class PengRobinsonModel;				/// Peng-Robinson (real-gas) thermodynamic model

////////// FUNCTION DECLARATION //////////


////////// BaseThermodynamicModel CLASS //////////
class BaseThermodynamicModel {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        BaseThermodynamicModel();					/// Default constructor
        BaseThermodynamicModel(const string configuration_file);	/// Parametrized constructor
        virtual ~BaseThermodynamicModel();				/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getSpecificGasConstant() { return( R_specific ); };

	////////// SET FUNCTIONS //////////
        inline void setSpecificGasConstant(double R_specific_) { R_specific = R_specific_; };

	////////// METHODS //////////
       
        /// Read configuration (input) file written in YAML language
        virtual void readConfigurationFile() {};
 
        /// Calculate pressure and temperature from density and internal energy
        virtual void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {};

        /// Calculate density and internal energy from pressure and temperature
        virtual void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {};

        /// Calculate specific heat capacities
        virtual void calculateSpecificHeatCapacities(double &c_v, double &c_p) {};

        /// Calculate heat capacities ratio
        virtual double calculateHeatCapacitiesRatio() = 0;

        /// Calculate speed of sound
        virtual double calculateSoundSpeed(const double &rho, const double &P, const double &T) = 0;

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double R_universal = 8.31446261815324;			/// Universal (ideal-gas) gas constant [J/(mol·K)]
        double R_specific;					/// Specific gas constant [J/(kg·K)]

        /// Model parameters
        string configuration_file;				/// Configuration file name (YAML language)

    private:

};


////////// IdealGasModel CLASS //////////
class IdealGasModel : public BaseThermodynamicModel {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        IdealGasModel();							/// Default constructor
        IdealGasModel(const string configuration_file);				/// Parametrized constructor
        virtual ~IdealGasModel();						/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p);

        /// Calculate heat capacities ratio
        virtual double calculateHeatCapacitiesRatio();

        /// Calculate speed of sound
        double calculateSoundSpeed(const double &rho, const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double gamma;						/// Heat capacities ratio [-]

    private:

};


////////// StiffenedGasModel CLASS //////////
class StiffenedGasModel : public BaseThermodynamicModel {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        StiffenedGasModel();							/// Default constructor
        StiffenedGasModel(const string configuration_file);			/// Parametrized constructor
        virtual ~StiffenedGasModel();						/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v_, double &c_p);

        /// Calculate heat capacities ratio
        virtual double calculateHeatCapacitiesRatio();

        /// Calculate speed of sound
        double calculateSoundSpeed(const double &rho, const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double gamma;						/// Heat capacities ratio [-]
        double P_inf;						/// Pressure infinity (liquid stiffnes) [Pa]
        double e_0;						/// Internal energy zero reference [J/kg]
        double c_v;						/// Specific isochoric heat capacity [J/(kg·K)]

    private:

};


////////// PengRobinsonModel CLASS //////////
class PengRobinsonModel : public BaseThermodynamicModel {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        PengRobinsonModel();							/// Default constructor
        PengRobinsonModel(const string configuration_file);			/// Parametrized constructor
        virtual ~PengRobinsonModel();						/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p);

        /// Calculate heat capacities ratio
        virtual double calculateHeatCapacitiesRatio();

        /// Calculate speed of sound
        double calculateSoundSpeed(const double &rho, const double &P, const double &T);

        /// Calculate attractive-forces a coefficient
        double calculate_eos_a(const double &T);

        /// Calculate first derivative of attractive-forces a coefficient
        double calculate_eos_a_first_derivative(const double &T);

        /// Calculate second derivative of attractive-forces a coefficient
        double calculate_eos_a_second_derivative(const double &T);

        /// Calculate compressibility factor
        double calculate_Z(const double &P, const double &T, const double &v);

        /// Calculate auxiliar parameters
        double calculate_A(const double &P, const double &T);
        double calculate_B(const double &P, const double &T);
        double calculate_M(const double &Z, const double &B);
        double calculate_N(const double &eos_a_first_derivative, const double &B);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double atmospheric_pressure = 101325.0;			/// Atmospheric pressure [Pa]
        double molecular_weight;				/// Molecular weight [kg/mol]
        double acentric_factor;					/// Acentric factor [-]
        double critical_temperature;				/// Critical temperature [K]
        double critical_pressure;				/// Critical pressure [Pa]
        double critical_molar_volume;				/// Critical molar volume [m3/mol]
        double NASA_coefficients[15];				/// NASA 7-coefficient polynomial

        /// Equation of state (EoS) properties
        double eos_d1 = 2.0;					/// EoS constant coefficient (fixed)
        double eos_d2 = -1.0;					/// EoS constant coefficient (fixed)
        //double eos_a;						/// EoS attractive-forces coefficient
        double eos_b;						/// EoS finite-pack-volume coefficient
        double eos_ac, eos_kappa;				/// EoS dimensionless parameters


    private:

};

#endif /*_THERMODYNAMIC_MODELS_*/

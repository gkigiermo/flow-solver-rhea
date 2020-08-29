#ifndef _THERMODYNAMIC_MODELS_
#define _THERMODYNAMIC_MODELS_

////////// INCLUDES //////////
#include <stdio.h>
#include <cmath>
#include <iostream>
#include "yaml-cpp/yaml.h"

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class BaseThermodynamicModel;				/// Base thermodynamic model
class IdealGasModel;					/// Ideal-gas thermodynamic model
class StiffenedGasModel;				/// Stiffened-gas thermodynamic model

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
        inline double getHeatCapacitiesRatio() { return( gamma ); };

	////////// SET FUNCTIONS //////////
        inline void setSpecificGasConstant(double R_specific_) { R_specific = R_specific_; };
        inline void setHeatCapacitiesRatio(double gamma_) { gamma = gamma_; };

	////////// METHODS //////////
       
        /// Read configuration (input) file written in YAML language
        virtual void readConfigurationFile() {};
 
        /// Calculate point pressure and temperature from density and internal energy
        virtual void calculatePointPressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {};

        /// Calculate point density and internal energy from pressure and temperature
        virtual void calculatePointDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {};

        /// Calculate point specific heat capacities
        virtual void calculatePointSpecificHeatCapacities(double &c_v, double &c_p) {};

        /// Calculate point speed of sound
        virtual double calculatePointSoundSpeed(const double &rho, const double &P, const double &T) = 0;

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double R_specific;					/// Specific gas constant [J/(kg·K)]
        double gamma;						/// Heat capacities ratio [-]

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

        /// Calculate point pressure and temperature from density and internal energy
        void calculatePointPressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate point density and internal energy from pressure and temperature
        void calculatePointDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate point specific heat capacities
        void calculatePointSpecificHeatCapacities(double &c_v, double &c_p);

        /// Calculate point speed of sound
        double calculatePointSoundSpeed(const double &rho, const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

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

        /// Calculate point pressure and temperature from density and internal energy
        void calculatePointPressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate point density and internal energy from pressure and temperature
        void calculatePointDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate point specific heat capacities
        void calculatePointSpecificHeatCapacities(double &c_v_, double &c_p);

        /// Calculate point speed of sound
        double calculatePointSoundSpeed(const double &rho, const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double P_inf;						/// Pressure infinity (liquid stiffnes) [Pa]
        double e_0;						/// Internal energy zero point [J/kg]
        double c_v;						/// Specific isochoric heat capacity [J/(kg·K)]

    private:

};

#endif /*_THERMODYNAMIC_MODELS_*/

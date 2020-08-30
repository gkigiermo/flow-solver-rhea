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
 
        /// Calculate point pressure and temperature from density and internal energy
        virtual void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {};

        /// Calculate point density and internal energy from pressure and temperature
        virtual void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {};

        /// Calculate point specific heat capacities
        virtual void calculateSpecificHeatCapacities(double &c_v, double &c_p) {};

        /// Calculate point speed of sound
        virtual double calculateSoundSpeed(const double &rho, const double &P, const double &T) = 0;

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
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

        /// Calculate point pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate point density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate point specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p);

        /// Calculate point speed of sound
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

        /// Calculate point pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate point density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate point specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v_, double &c_p);

        /// Calculate point speed of sound
        double calculateSoundSpeed(const double &rho, const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double gamma;						/// Heat capacities ratio [-]
        double P_inf;						/// Pressure infinity (liquid stiffnes) [Pa]
        double e_0;						/// Internal energy zero point [J/kg]
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

        /// Calculate point pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate point density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate point specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p);

        /// Calculate point speed of sound
        double calculateSoundSpeed(const double &rho, const double &P, const double &T);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties

    private:

};

#endif /*_THERMODYNAMIC_MODELS_*/

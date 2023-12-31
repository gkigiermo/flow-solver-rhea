#ifndef _THERMODYNAMIC_MODELS_
#define _THERMODYNAMIC_MODELS_

////////// INCLUDES //////////
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <complex>
#include <mpi.h>
#include "yaml-cpp/yaml.h"
#include "RootFindingMinimization.hpp"

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
        BaseThermodynamicModel(const std::string configuration_file);	/// Parametrized constructor
        virtual ~BaseThermodynamicModel();				/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getSpecificGasConstant() { return( R_specific ); };
        inline double getMolecularWeight() { return( molecular_weight ); };

	////////// SET FUNCTIONS //////////
        inline void setSpecificGasConstant(double R_specific_) { R_specific = R_specific_; };
        inline void setMolecularWeight(double molecular_weight_) { molecular_weight = molecular_weight_; };

	////////// METHODS //////////
       
        /// Read configuration (input) file written in YAML language
        virtual void readConfigurationFile() {};

        /// Calculate temperature from pressure and density
        virtual double calculateTemperatureFromPressureDensity(const double &P, const double &rho) = 0;

        /// Calculate temperature from pressure and density with initial guess
        virtual void calculateTemperatureFromPressureDensityWithInitialGuess(double &T, const double &P, const double &rho) {};

        /// Calculate internal energy from pressure, temperature and density
        virtual double calculateInternalEnergyFromPressureTemperatureDensity(const double &P, const double &T, const double &rho) = 0;

        /// Calculate pressure and temperature from density and internal energy
        virtual void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e) {};
        
        /// Calculate density from pressure and temperature
        virtual double calculateDensityFromPressureTemperature(const double &P, const double &T) = 0;

        /// Calculate density and internal energy from pressure and temperature
        virtual void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T) {};

        /// Calculate specific heat capacities
        virtual void calculateSpecificHeatCapacities(double &c_v, double &c_p, const double &P, const double &T, const double &rho) {};

        /// Calculate heat capacities ratio
        virtual double calculateHeatCapacitiesRatio(const double &P, const double &rho) = 0;

        /// Calculate speed of sound
        virtual double calculateSoundSpeed(const double &P, const double &T, const double &rho) = 0;

        /// Calculate expansivity & compressibility
        virtual double calculateVolumeExpansivity(const double &T, const double &bar_v) = 0;
        virtual double calculateIsothermalCompressibility(const double &T, const double &bar_v) = 0;
        virtual double calculateIsentropicCompressibility(const double &P, const double &T, const double &bar_v) = 0;

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double R_universal = 8.31446261815324;					/// Universal (ideal-gas) gas constant [J/(mol·K)]
        double R_specific;							/// Specific gas constant [J/(kg·K)]
        double molecular_weight;						/// Molecular weight [kg/mol]

        /// Model parameters
	std::string configuration_file;						/// Configuration file name (YAML language)
        
    private:

};


////////// IdealGasModel CLASS //////////
class IdealGasModel : public BaseThermodynamicModel {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        IdealGasModel();							/// Default constructor
        IdealGasModel(const std::string configuration_file);			/// Parametrized constructor
        virtual ~IdealGasModel();						/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate temperature from pressure and density
        double calculateTemperatureFromPressureDensity(const double &P, const double &rho);

        /// Calculate temperature from pressure and density with initial guess
        void calculateTemperatureFromPressureDensityWithInitialGuess(double &T, const double &P, const double &rho);

        /// Calculate internal energy from pressure, temperature and density
        double calculateInternalEnergyFromPressureTemperatureDensity(const double &P, const double &T, const double &rho);

        /// Calculate pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate density from pressure and temperature
        double calculateDensityFromPressureTemperature(const double &P, const double &T);

        /// Calculate density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p, const double &P, const double &T, const double &rho);

        /// Calculate heat capacities ratio
        double calculateHeatCapacitiesRatio(const double &P, const double &rho);

        /// Calculate speed of sound
        double calculateSoundSpeed(const double &P, const double &T, const double &rho);

        /// Calculate expansivity & compressibility
        double calculateVolumeExpansivity(const double &T, const double &bar_v);
        double calculateIsothermalCompressibility(const double &T, const double &bar_v);
        double calculateIsentropicCompressibility(const double &P, const double &T, const double &bar_v);

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double gamma;						/// Heat capacities ratio (ideal-gas) [-]

    private:

};


////////// StiffenedGasModel CLASS //////////
class StiffenedGasModel : public BaseThermodynamicModel {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        StiffenedGasModel();							/// Default constructor
        StiffenedGasModel(const std::string configuration_file);		/// Parametrized constructor
        virtual ~StiffenedGasModel();						/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate temperature from pressure and density
        double calculateTemperatureFromPressureDensity(const double &P, const double &rho);

        /// Calculate temperature from pressure and density with initial guess
        void calculateTemperatureFromPressureDensityWithInitialGuess(double &T, const double &P, const double &rho);

        /// Calculate internal energy from pressure, temperature and density
        double calculateInternalEnergyFromPressureTemperatureDensity(const double &P, const double &T, const double &rho);

        /// Calculate pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate density from pressure and temperature
        double calculateDensityFromPressureTemperature(const double &P, const double &T);

        /// Calculate density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v_, double &c_p, const double &P, const double &T, const double &rho);

        /// Calculate heat capacities ratio
        double calculateHeatCapacitiesRatio(const double &P, const double &rho);

        /// Calculate speed of sound
        double calculateSoundSpeed(const double &P, const double &T, const double &rho);

        /// Calculate expansivity & compressibility
        double calculateVolumeExpansivity(const double &T, const double &bar_v);
        double calculateIsothermalCompressibility(const double &T, const double &bar_v);
        double calculateIsentropicCompressibility(const double &P, const double &T, const double &bar_v);

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
        PengRobinsonModel(const std::string configuration_file);		/// Parametrized constructor
        virtual ~PengRobinsonModel();						/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getCriticalPressure() { return( critical_pressure ); };
        inline double getCriticalTemperature() { return( critical_temperature ); };

	////////// SET FUNCTIONS //////////
        inline void setCriticalPressure(double critical_pressure_) { critical_pressure = critical_pressure_; };
        inline void setCriticalTemperature(double critical_temperature_) { critical_temperature = critical_temperature_; };

	////////// METHODS //////////
        
        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Calculate temperature from pressure and density
        double calculateTemperatureFromPressureDensity(const double &P, const double &rho);

        /// Calculate temperature from pressure and density with initial guess
        void calculateTemperatureFromPressureDensityWithInitialGuess(double &T, const double &P, const double &rho);

        /// Calculate internal energy from pressure, temperature and density
        double calculateInternalEnergyFromPressureTemperatureDensity(const double &P, const double &T, const double &rho);

        /// Calculate pressure and temperature from density and internal energy
        void calculatePressureTemperatureFromDensityInternalEnergy(double &P, double &T, const double &rho, const double &e);

        /// Calculate density from pressure and temperature
        double calculateDensityFromPressureTemperature(const double &P, const double &T);

        /// Calculate density and internal energy from pressure and temperature
        void calculateDensityInternalEnergyFromPressureTemperature(double &rho, double &e, const double &P, const double &T);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p, const double &P, const double &T, const double &rho);

        /// Calculate heat capacities ratio
        double calculateHeatCapacitiesRatio(const double &P, const double &rho);

        /// Calculate speed of sound
        double calculateSoundSpeed(const double &P, const double &T, const double &rho);

        /// Calculate expansivity & compressibility
        double calculateVolumeExpansivity(const double &T, const double &bar_v);
        double calculateIsothermalCompressibility(const double &T, const double &bar_v);
        double calculateIsentropicCompressibility(const double &P, const double &T, const double &bar_v);

        /// Calculate pressure from temperature and density
        double calculatePressureFromTemperatureDensity(const double &T, const double &rho);

        /// Calculate molar internal energy from pressure, temperature and molar volume
        double calculateMolarInternalEnergyFromPressureTemperatureMolarVolume(const double &P, const double &T, const double &bar_v);

        /// Calculate attractive-forces a coefficient
        double calculate_eos_a(const double &T);

        /// Calculate first derivative of attractive-forces a coefficient
        double calculate_eos_a_first_derivative(const double &T);

        /// Calculate second derivative of attractive-forces a coefficient
        double calculate_eos_a_second_derivative(const double &T);

        /// Calculate compressibility factor
        double calculate_Z(const double &P, const double &T, const double &bar_v);

        /// Calculate auxiliar parameters
        double calculate_A(const double &P, const double &T);
        double calculate_B(const double &P, const double &T);
        double calculate_M(const double &Z, const double &B);
        double calculate_N(const double &eos_a_first_derivative, const double &B);

        /// Calculate standard thermodynamic variables from NASA 7-coefficient polynomial
        double calculateMolarStdCpFromNASApolynomials(const double &T);
        double calculateMolarStdEnthalpyFromNASApolynomials(const double &T);

        /// Calculate high-pressure departure functions
        double calculateDepartureFunctionMolarCp(const double &P, const double &T, const double &bar_v);
        double calculateDepartureFunctionMolarCv(const double &P, const double &T, const double &bar_v);
        double calculateDepartureFunctionMolarEnthalpy(const double &P, const double &T, const double &bar_v);

        /// Calculate temperature from pressure and molar volume
        double calculateTemperatureFromPressureMolarVolume(const double &P, const double &bar_v);

        /// Calculate thermodynamic derivatives
        double calculateDPDTConstantMolarVolume(const double &T, const double &bar_v);
        double calculateDPDvConstantTemperature(const double &T, const double &bar_v);

        /// Calculate roots of cubic polynomial
        void calculateRootsCubicPolynomial(std::complex<double> &root_1, std::complex<double> &root_2, std::complex<double> &root_3, double &a, double &b, double &c, double &d);

    private:

        /// Nonlinear solver nested class used to obtain P & T from rho & e
        class NLS_P_T_from_rho_e : public NewtonRaphson { 

            public:

            ////////// CONSTRUCTORS & DESTRUCTOR //////////
            NLS_P_T_from_rho_e(std::vector<double> &fvec_, PengRobinsonModel &pr_model_) : NewtonRaphson(fvec_), pr_model(pr_model_) {};	/// Parametrized constructor
            virtual ~NLS_P_T_from_rho_e() {};													/// Destructor

	    ////////// METHODS //////////

            /// Set external (enclosing class) parameters
            void setExternalParameters(const double &target_rho_, const double &target_e_, const double &P_norm_, const double &T_norm_) {

                target_rho = target_rho_;
                target_e   = target_e_;
                P_norm     = P_norm_;
                T_norm     = T_norm_;

            };

            /// Evaluates the functions (residuals) in position x
            void function_vector(std::vector<double> &x, std::vector<double> &fx) {

                /// Set state to x
                double P = x[0]*P_norm;		// Unnormalize pressure
                double T = x[1]*T_norm;		// Unnormalize temperature
 
                /// For single-component systems:
                /// P will not oscillate for supercritical thermodynamic states
                /// P will oscillate for subcritical thermodynamic states
                /// ... small oscillations close to the critical point (slightly subcritical)
                /// ... large oscillations far from the critical point (notably subcritical) -> in that case, use two-phase solver
                double molecular_weight = pr_model.getMolecularWeight(); 
                double target_molar_v   = molecular_weight/target_rho;
                double guess_rho = pr_model.calculateDensityFromPressureTemperature( P, T );
                double guess_e   = pr_model.calculateMolarInternalEnergyFromPressureTemperatureMolarVolume( P, T, target_molar_v )/molecular_weight;
 
                /// Compute fx (residuals)
                fx[0] = ( guess_rho - target_rho )/( fabs( target_rho ) + 1.0e-14 );	/// function normalized
                fx[1] = ( guess_e - target_e )/( fabs( target_e ) + 1.0e-14 );		/// function normalized

            };
      
            ////////// PARAMETERS //////////

	    /// External (enclosing) parameters
            double target_rho;			/// External (enclosing class) target density
            double target_e;			/// External (enclosing class) target specific internal energy
            double P_norm;			/// External (enclosing class) pressure normalization
            double T_norm;			/// External (enclosing class) temperature normalization

	    /// External (enclosing) class
            PengRobinsonModel &pr_model;	/// Reference to PengRobinsonModel (enclosing) class
 
        };

        /// Nonlinear solver nested class used to obtain T from P & rho
        class NLS_T_from_P_rho : public Brent { 

            public:

            ////////// CONSTRUCTORS & DESTRUCTOR //////////
            NLS_T_from_P_rho(std::vector<double> &fvec_, PengRobinsonModel &pr_model_) : Brent(fvec_), pr_model(pr_model_) {};	/// Parametrized constructor
            virtual ~NLS_T_from_P_rho() {};											/// Destructor

	    ////////// METHODS //////////

            /// Set external (enclosing class) parameters
            void setExternalParameters(const double &target_P_, const double &target_rho_, const double &T_norm_) {

                target_P   = target_P_;
                target_rho = target_rho_;
                T_norm     = T_norm_;

            };

            /// Evaluates the functions (residuals) in position x
            void function_vector(std::vector<double> &x, std::vector<double> &fx) {

                /// Set state to x
                double T = x[0]*T_norm;		// Unnormalize temperature

                /// For single-component systems:
                /// P will not oscillate for supercritical thermodynamic states
                /// P will oscillate for subcritical thermodynamic states
                /// ... small oscillations close to the critical point (slightly subcritical)
                /// ... large oscillations far from the critical point (notably subcritical) -> in that case, use two-phase solver
                double guess_P = pr_model.calculatePressureFromTemperatureDensity( T, target_rho );

                /// Compute fx (residuals)
                fx[0] = ( ( guess_P - target_P )/( target_P + 1.0e-14 ) );		/// function normalized
        
            };
      
            ////////// PARAMETERS //////////

	    /// External (enclosing) parameters
            double target_P;			/// External (enclosing class) target pressure
            double target_rho;			/// External (enclosing class) target density
            double T_norm;			/// External (enclosing class) temperature normalization

	    /// External (enclosing) class
            PengRobinsonModel &pr_model;	/// Reference to PengRobinsonModel (enclosing) class
 
        };

    protected:

        ////////// PARAMETERS //////////

        /// Thermodynamic properties
        double acentric_factor;					/// Acentric factor [-]
        double critical_temperature;				/// Critical temperature [K]
        double critical_pressure;				/// Critical pressure [Pa]
        double critical_molar_volume;				/// Critical molar volume [m3/mol]
        double NASA_coefficients[15];				/// NASA 7-coefficient polynomial

        /// Equation of state (EoS) parameters
        //double eos_a;						/// EoS attractive-forces coefficient
        double eos_b;						/// EoS finite-pack-volume coefficient
        double eos_ac, eos_kappa;				/// EoS dimensionless parameters

        /// Aitken's delta-squared process parameters
        int aitken_max_iter              = 1000;		/// Maximum number of iterations
        double aitken_relative_tolerance = 1.0e-5;		/// Relative tolerance

        /// Nonlinear P-T solver parameters
        int nls_PT_max_iter              = 1000;		/// Maximum number of iterations
        double nls_PT_relative_tolerance = 1.0e-5;		/// Relative tolerance
        NLS_P_T_from_rho_e *nls_PT_solver;			/// Pointer to NLS_P_T_from_rho_e
	std::vector<double> nls_PT_unknowns;			/// NLS_P_T_from_rho_e unknowns: P & T
	std::vector<double> nls_PT_r_vec;			/// NLS_P_T_from_rho_e vector of functions residuals

        /// Nonlinear T solver parameters
        int nls_T_max_iter              = 1000;			/// Maximum number of iterations
        double nls_T_relative_tolerance = 1.0e-5;		/// Relative tolerance
        NLS_T_from_P_rho *nls_T_solver;				/// Pointer to NLS_T_from_P_rho
	std::vector<double> nls_T_unknowns;			/// NLS_T_from_P_rho unknowns: T
	std::vector<double> nls_T_r_vec;			/// NLS_T_from_P_rho vector of functions residuals

    private:

};

#endif /*_THERMODYNAMIC_MODELS_*/

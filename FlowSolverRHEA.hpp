#ifndef _RHEA_FLOW_SOLVER_
#define _RHEA_FLOW_SOLVER_

////////// INCLUDES //////////
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <limits>
#include <iomanip>
#include <mpi.h>
#include "yaml-cpp/yaml.h"
#include "src/MacroParameters.hpp"
#include "src/ComputationalDomain.hpp"
#include "src/ParallelTopology.hpp"
#include "src/DistributedArray.hpp"
#include "src/ManagerHDF5.hpp"

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class FlowSolverRHEA;						/// Flow solver RHEA

////////// FUNCTION DECLARATION //////////


////////// FlowSolverRHEA CLASS //////////
class FlowSolverRHEA {
   
    ////////// VARIABLES & PARAMETERS DESCRIPTION //////////

    /// Primitive variables:
    ///   - Density: rho
    ///   - Velocities: u, v, w
    ///   - Specific total energy: E = e + ke
    ///     ... sum of specific internal energy e, and specific kinetic energy ke = (u*u + v*v + w*w)/2

    /// Conserved variables:
    ///   - Mass: rho
    ///   - Momentum: rho*u, rho*v, rho*w
    ///   - Total energy: rho*E

    /// Thermodynamic state:
    ///   - Pressure: P
    ///   - Temperature: T
    ///   - Speed of sound: sos
    ///   - Specific gas constant: R_specific
    ///   - Ratio of heat capacities: gamma

    /// Thermophysical properties:
    ///   - Dynamic viscosity: mu
    ///   - Thermal conductivity: kappa
 
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        FlowSolverRHEA();					/// Default constructor
        FlowSolverRHEA(const string configuration_file);	/// Parametrized constructor
        virtual ~FlowSolverRHEA();				/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getCurrentTime() { return( current_time ); };
        inline double getFinalTime() { return( final_time ); };
        inline double getTimeStep() { return( delta_t ); };
        inline int getCurrentTimeIteration() { return( current_time_iter ); };
        inline int getFinalTimeIteration() { return( final_time_iter ); };
        inline int getCurrentRungeKuttaStep() { return( rk_step ); };
        inline int getRungeKuttaOrder() { return( rk_order ); };
        inline int getOutputIterationFrequency() { return( output_frequency_iter ); };
        inline bool getUseRestart() { return( use_restart ); };

	////////// SET FUNCTIONS //////////
        inline void setCurrentTime(double current_time_) { current_time = current_time_; };
        inline void setFinalTime(double final_time_) { final_time = final_time_; };
        inline void setTimeStep(double delta_t_) { delta_t = delta_t_; };
        inline void setCurrentTimeIteration(int current_time_iter_) { current_time_iter = current_time_iter_; };
        inline void setCurrentRungeKuttaStep(int rk_step_) { rk_step = rk_step_; };
        inline void setOutputIterationFrequency(int output_frequency_iter_) { output_frequency_iter = output_frequency_iter_; };
        inline void setUseRestart(bool use_restart_) { use_restart = use_restart_; };

	////////// SOLVER METHODS //////////

        /// Read configuration (input) file written in YAML language
        void readConfigurationFile();

        /// Fill x, y and z mesh coordinate fields (input/output data)
        void fillMeshCoordinateFields();

        /// Set initial conditions: u, v, w, P and T ... needs to be modified/overwritten according to the problem under consideration
        void setInitialConditions();

        /// Initialize from restart file: x, y, z, rho, u, v, w, E, P, T, sos, mu and kappa
        void initializeFromRestart();

        /// Initialize thermodynamic state: rho, E and sos
        void initializeThermodynamics();

        /// Calculate conserved variables from primitive variables: variable -> rho*variable
        void primitiveToConservedVariables();

        /// Calculate primitive variables from conserved variables: rho*variable -> variable
        void conservedToPrimitiveVariables();

        /// Calculate thermodynamics from primitive variables
        void calculateThermodynamicsFromPrimitiveVariables();

        /// Calculate point density and internal energy from pressure and temperature
        void calculatePointDensityInternalEnergyFromPressureTemperature(const double &P, const double &T, double &rho, double &e);

        /// Calculate specific heat capacities
        void calculateSpecificHeatCapacities(double &c_v, double &c_p);

        /// Update boundary values: rho, rhou, rhov, rhow and rhoE
        void updateBoundaries();

        /// Update previous state of conserved variables: rho*variable -> rho_0*variable_0
        void updatePreviousStateConservedVariables();

        /// Calculate time step satisfying CFL constraint
        void calculateTimeStep();

        /// Calculate thermophysical properties
        void calculateThermophysicalProperties();

        /// Calculate source terms ... needs to be modified/overwritten according to the problem under consideration
        void calculateSourceTerms();

        /// Calculate waves speed
        void calculateWavesSpeed(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, double &S_L, double &S_R);

        /// Calculate HLLC flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateHllcFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

        /// Calculate inviscid fluxes
        void calculateInviscidFluxes();

        /// Calculate viscous fluxes
        void calculateViscousFluxes();

        /// Sum inviscid & viscous fluxes, and source terms
        void sumInviscidViscousFluxesSourceTerms();

        /// Advance conserved variables in time
        void timeAdvanceConservedVariables();

        /// Output current solver state data
        void outputCurrentStateData();

    protected:

        ////////// SOLVER PARAMETERS //////////
	
        /// Fluid properties 
        double R_specific;					/// Specific gas constant [J/(kg·K)]
        double gamma;						/// Heat capacity ratio [-]
        double mu;						/// Dynamic viscosity [Pa·s]
        double kappa;						/// Thermal conductivity [W/(m·k)]       

        /// Problem parameters
        double x_0;           	        	 	 	/// Domain origin in x-direction [m]
        double y_0;             	        		/// Domain origin in y-direction [m]
        double z_0;                     			/// Domain origin in z-direction [m]
        double L_x;   						/// Domain size in x-direction [m]
        double L_y;      					/// Domain size in y-direction [m]
        double L_z;						/// Domain size in z-direction [m]
        double current_time;   					/// Current time [s]
        double final_time;		      			/// Final time [s]
        string configuration_file;				/// Configuration file name (YAML language)	

        /// Computational parameters
        int num_grid_x;						/// Number of inner grid points in x-direction
        int num_grid_y;						/// Number of inner grid points in y-direction
        int num_grid_z;						/// Number of inner grid points in z-direction
        double A_x;						/// Stretching factor in x-direction
        double A_y;						/// Stretching factor in y-direction
        double A_z;						/// Stretching factor in z-direction
        double CFL;						/// CFL coefficient
        double delta_t;		      				/// Time step [s]
        int current_time_iter;					/// Current time iteration
        int final_time_iter;					/// Final time iteration
        int rk_step;						/// Current Runge-Kutta step: 1, 2, 3
        int rk_order;						/// Order of Runge-Kutta method (fixed)

        /// Local mesh values for I1D macro
        int _lNx_;
        int _lNy_;
        int _lNz_;

        /// Boundary conditions
        int bocos_type[6];					/// Array of boundary conditions type
        double bocos_u[6];					/// Array of boundary conditions u
        double bocos_v[6];					/// Array of boundary conditions v
        double bocos_w[6];					/// Array of boundary conditions w
        double bocos_P[6];					/// Array of boundary conditions P
        double bocos_T[6];					/// Array of boundary conditions T

        /// Write/read file parameters
        string output_data_file;				/// Output data file name (HDF5 format)	
        int output_frequency_iter;				/// Data output iteration frequency
        bool generate_xdmf;					/// Generate xdmf file reader
        bool use_restart;					/// Use restart file for initialization
        int restart_data_file_iter;				/// Restart data file iteration

        /// Parallelization scheme
        int np_x;						/// Number of processes in x-direction
        int np_y;						/// Number of processes in y-direction
        int np_z;						/// Number of processes in z-direction

	////////// SOLVER (PARALLEL) VARIABLES //////////
	
        /// Mesh coordinates (input/output data)
        DistributedArray x_field;				/// 3-D field of x-coordinate
        DistributedArray y_field;				/// 3-D field of y-coordinate
        DistributedArray z_field;				/// 3-D field of z-coordinate

        /// Primitive, conserved, thermodynamic and thermophysical variables
        DistributedArray rho_field;				/// 3-D field of rho
        DistributedArray u_field;				/// 3-D field of u
        DistributedArray v_field;				/// 3-D field of v
        DistributedArray w_field;				/// 3-D field of w
        DistributedArray E_field;				/// 3-D field of E
        DistributedArray rhou_field;				/// 3-D field of rhou
        DistributedArray rhov_field;				/// 3-D field of rhov
        DistributedArray rhow_field;				/// 3-D field of rhow
        DistributedArray rhoE_field;				/// 3-D field of rhoE
        DistributedArray P_field;				/// 3-D field of P
        DistributedArray T_field;				/// 3-D field of T
        DistributedArray sos_field;				/// 3-D field of sos
        DistributedArray mu_field;				/// 3-D field of mu
        DistributedArray kappa_field;				/// 3-D field of kappa

        /// Time-integration variables
        DistributedArray rho_0_field;				/// 3-D previous field of rho
        DistributedArray rhou_0_field;				/// 3-D previous field of rhou
        DistributedArray rhov_0_field;				/// 3-D previous field of rhov
        DistributedArray rhow_0_field;				/// 3-D previous field of rhow
        DistributedArray rhoE_0_field;				/// 3-D previous field of rhoE

        /// Time-integration fluxes
        DistributedArray rho_rk1_flux;				/// 3-D Runge-Kutta flux 1 of rho
        DistributedArray rho_rk2_flux;				/// 3-D Runge-Kutta flux 2 of rho
        DistributedArray rho_rk3_flux;				/// 3-D Runge-Kutta flux 3 of rho
        DistributedArray rhou_rk1_flux;				/// 3-D Runge-Kutta flux 1 of rhou
        DistributedArray rhou_rk2_flux;				/// 3-D Runge-Kutta flux 2 of rhou
        DistributedArray rhou_rk3_flux;				/// 3-D Runge-Kutta flux 3 of rhou
        DistributedArray rhov_rk1_flux;				/// 3-D Runge-Kutta flux 1 of rhov
        DistributedArray rhov_rk2_flux;				/// 3-D Runge-Kutta flux 2 of rhov
        DistributedArray rhov_rk3_flux;				/// 3-D Runge-Kutta flux 3 of rhov
        DistributedArray rhow_rk1_flux;				/// 3-D Runge-Kutta flux 1 of rhow
        DistributedArray rhow_rk2_flux;				/// 3-D Runge-Kutta flux 2 of rhow
        DistributedArray rhow_rk3_flux;				/// 3-D Runge-Kutta flux 3 of rhow
        DistributedArray rhoE_rk1_flux;				/// 3-D Runge-Kutta flux 1 of rhoE
        DistributedArray rhoE_rk2_flux;				/// 3-D Runge-Kutta flux 2 of rhoE
        DistributedArray rhoE_rk3_flux;				/// 3-D Runge-Kutta flux 3 of rhoE

        /// Inviscid fluxes
        DistributedArray rho_inv_flux;				/// 3-D inviscid fluxes of rho
        DistributedArray rhou_inv_flux;				/// 3-D inviscid fluxes of rhou
        DistributedArray rhov_inv_flux;				/// 3-D inviscid fluxes of rhov
        DistributedArray rhow_inv_flux;				/// 3-D inviscid fluxes of rhow
        DistributedArray rhoE_inv_flux;				/// 3-D inviscid fluxes of rhoE

        /// Viscous fluxes
        DistributedArray rhou_vis_flux;				/// 3-D viscous fluxes of rhou
        DistributedArray rhov_vis_flux;				/// 3-D viscous fluxes of rhov
        DistributedArray rhow_vis_flux;				/// 3-D viscous fluxes of rhow
        DistributedArray rhoE_vis_flux;				/// 3-D viscous fluxes of rhoE

        /// Source terms
        DistributedArray f_rhou_field;				/// 3-D field of rhou
        DistributedArray f_rhov_field;				/// 3-D field of rhov
        DistributedArray f_rhow_field;				/// 3-D field of rhow
        DistributedArray f_rhoE_field;				/// 3-D field of rhoE

	////////// COMPUTATIONAL DOMAIN, PARALLEL TOPOLOGY & WRITER/READER //////////
        ComputationalDomain *mesh;				/// Computational domain
        ParallelTopology *topo;					/// Parallel topology
        ManagerHDF5 *writer_reader;				/// HDF5 data writer/reader

    private:

};

#endif /*_RHEA_FLOW_SOLVER_*/

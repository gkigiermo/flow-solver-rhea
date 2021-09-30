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
#include "MacroParameters.hpp"
#include "ThermodynamicModel.hpp"
#include "TransportCoefficients.hpp"
#include "ComputationalDomain.hpp"
#include "ParallelTopology.hpp"
#include "DistributedArray.hpp"
#include "ManagerHDF5.hpp"
#include "ParallelTimer.hpp"

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class FlowSolverRHEA;					/// Flow solver RHEA
class BaseRiemannSolver;				/// Base Riemann solver
class CentralFluxApproximateRiemannSolver;		/// Central scheme approximate Riemann solver
class MurmanRoeFluxApproximateRiemannSolver;		/// Murman-Roe scheme approximate Riemann solver
class KgpFluxApproximateRiemannSolver;			/// KGP scheme approximate Riemann solver
class HllApproximateRiemannSolver;			/// HLL approximate Riemann solver
class HllcApproximateRiemannSolver;			/// HLLC approximate Riemann solver
class HllcPlusApproximateRiemannSolver;			/// HLLC+ approximate Riemann solver

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

    /// Transport coefficients:
    ///   - Dynamic viscosity: mu
    ///   - Thermal conductivity: kappa
 
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        FlowSolverRHEA();					/// Default constructor
        FlowSolverRHEA(const string configuration_file);	/// Parametrized constructor
        virtual ~FlowSolverRHEA();				/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getCurrentTime() { return( current_time ); };
        inline int getCurrentTimeIteration() { return( current_time_iter ); };

	////////// SET FUNCTIONS //////////
        inline void setCurrentTime(double current_time_) { current_time = current_time_; };
        inline void setCurrentTimeIteration(int current_time_iter_) { current_time_iter = current_time_iter_; };

	////////// SOLVER METHODS //////////
        
	/// Execute (aggregated method) RHEA
        virtual void execute();

        /// Read configuration (input) file written in YAML language
        virtual void readConfigurationFile();

        /// Fill x, y and z mesh coordinate fields (input/output data)
        virtual void fillMeshCoordinateFields();

        /// Set initial conditions: u, v, w, P and T ... needs to be modified/overwritten according to the problem under consideration
        virtual void setInitialConditions();

        /// Initialize from restart file: x, y, z, rho, u, v, w, E, P, T, sos, mu and kappa
        virtual void initializeFromRestart();

        /// Initialize thermodynamic state: rho, E and sos
        virtual void initializeThermodynamics();

        /// Calculate conserved variables from primitive variables: variable -> rho*variable
        virtual void primitiveToConservedVariables();

        /// Calculate primitive variables from conserved variables: rho*variable -> variable
        virtual void conservedToPrimitiveVariables();

        /// Calculate thermodynamics from primitive variables
        virtual void calculateThermodynamicsFromPrimitiveVariables();

        /// Update boundary values: rho, rhou, rhov, rhow and rhoE
        virtual void updateBoundaries();

        /// Update previous state of conserved variables: rho*variable -> rho_0*variable_0
        virtual void updatePreviousStateConservedVariables();

        /// Calculate time step satisfying CFL constraint
        virtual void calculateTimeStep();

        /// Calculate transport coefficients
        virtual void calculateTransportCoefficients();

        /// Calculate rhou, rhov, rhow and rhoE source terms ... needs to be modified/overwritten according to the problem under consideration
        virtual void calculateSourceTerms();

        /// Calculate inviscid fluxes
        virtual void calculateInviscidFluxes();

        /// Calculate viscous fluxes
        virtual void calculateViscousFluxes();

        /// Sum inviscid & viscous fluxes, and source terms
        virtual void sumInviscidViscousFluxesSourceTerms(const int &rk_step);

        /// Advance conserved variables in time
        virtual void timeAdvanceConservedVariables(const int &rk_step);

        /// Output current solver state data
        virtual void outputCurrentStateData();

        /// Update time-averaged quantities
        virtual void updateTimeAveragedQuantities();

    protected:

        ////////// SOLVER PARAMETERS //////////
	
        /// Fluid & flow properties 
        string thermodynamic_model;				/// Thermodynamic model
        string transport_coefficients_model;			/// Transport coefficients model

        /// Problem parameters
        double x_0;           	        	 	 	/// Domain origin in x-direction [m]
        double y_0;             	        		/// Domain origin in y-direction [m]
        double z_0;                     			/// Domain origin in z-direction [m]
        double L_x;   						/// Domain size in x-direction [m]
        double L_y;      					/// Domain size in y-direction [m]
        double L_z;						/// Domain size in z-direction [m]
        double current_time;   					/// Current time [s]
        double final_time;		      			/// Final time [s]
        double averaging_time;		      			/// Averaging time [s]
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
        int rk_order = 3;					/// Order of Runge-Kutta method (fixed)
        string riemann_solver_scheme;				/// Riemann solver scheme

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

        /// Print/Write/Read file parameters
        int print_frequency_iter;				/// Print information iteration frequency
        string output_data_file_name;				/// Output data file name (HDF5 format)	
        int output_frequency_iter;				/// Data output iteration frequency
        bool generate_xdmf;					/// Generate xdmf file reader
        bool use_restart;					/// Use restart file for initialization
        string restart_data_file;				/// Restart data file
        bool time_averaging_active;				/// Activate time averaging
        bool reset_time_averaging;				/// Reset time averaging

        /// Timers information
        bool print_timers;					/// Print timers information
        string timers_information_file;				/// Timers information file

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
        DistributedArray c_v_field;				/// 3-D field of c_v
        DistributedArray c_p_field;				/// 3-D field of c_p

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

        /// Time averaging
        DistributedArray avg_rho_field;				/// 3-D field of time-averaged rho
        DistributedArray avg_rhou_field;			/// 3-D field of time-averaged rhou
        DistributedArray avg_rhov_field;			/// 3-D field of time-averaged rhov
        DistributedArray avg_rhow_field;			/// 3-D field of time-averaged rhow
        DistributedArray avg_rhoE_field;			/// 3-D field of time-averaged rhoE
        DistributedArray avg_u_field;				/// 3-D field of time-averaged u
        DistributedArray avg_v_field;				/// 3-D field of time-averaged v
        DistributedArray avg_w_field;				/// 3-D field of time-averaged w
        DistributedArray avg_E_field;				/// 3-D field of time-averaged E
        DistributedArray avg_P_field;				/// 3-D field of time-averaged P
        DistributedArray avg_T_field;				/// 3-D field of time-averaged T
        DistributedArray avg_sos_field;				/// 3-D field of time-averaged sos
        DistributedArray avg_c_v_field;				/// 3-D field of time-averaged c_v
        DistributedArray avg_c_p_field;				/// 3-D field of time-averaged c_p
        DistributedArray rmsf_rho_field;			/// 3-D field of root-mean-square-fluctuation rho
        DistributedArray rmsf_rhou_field;			/// 3-D field of root-mean-square-fluctuation rhou
        DistributedArray rmsf_rhov_field;			/// 3-D field of root-mean-square-fluctuation rhov
        DistributedArray rmsf_rhow_field;			/// 3-D field of root-mean-square-fluctuation rhow
        DistributedArray rmsf_rhoE_field;			/// 3-D field of root-mean-square-fluctuation rhoE
        DistributedArray rmsf_u_field;				/// 3-D field of root-mean-square-fluctuation u
        DistributedArray rmsf_v_field;				/// 3-D field of root-mean-square-fluctuation v
        DistributedArray rmsf_w_field;				/// 3-D field of root-mean-square-fluctuation w
        DistributedArray rmsf_E_field;				/// 3-D field of root-mean-square-fluctuation E
        DistributedArray rmsf_P_field;				/// 3-D field of root-mean-square-fluctuation P
        DistributedArray rmsf_T_field;				/// 3-D field of root-mean-square-fluctuation T
        DistributedArray rmsf_sos_field;			/// 3-D field of root-mean-square-fluctuation sos
        DistributedArray rmsf_c_v_field;			/// 3-D field of root-mean-square-fluctuation c_v
        DistributedArray rmsf_c_p_field;			/// 3-D field of root-mean-square-fluctuation c_p

	////////// THERMODYNAMIC MODEL, TRANSPORT COEFFICIENTS, RIEMANN SOLVER //////////
	////////// COMPUTATIONAL DOMAIN, PARALLEL TOPOLOGY, WRITER/READER, PARALLEL TIMER //////////
        BaseThermodynamicModel *thermodynamics;			/// Thermodynamic model
        BaseTransportCoefficients *transport_coefficients;	/// Transport coefficients
        BaseRiemannSolver *riemann_solver;			/// Riemann solver
        ComputationalDomain *mesh;				/// Computational domain
        ParallelTopology *topo;					/// Parallel topology
        ManagerHDF5 *writer_reader;				/// HDF5 data writer/reader
        ParallelTimer *timers;					/// Parallel timer

    private:

};

////////// BaseRiemannSolver CLASS //////////
class BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        BaseRiemannSolver();					/// Default constructor
        virtual ~BaseRiemannSolver();				/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
       
        /// Calculate waves speed
        virtual void calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R);

        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        virtual double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) = 0;

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// CentralFluxApproximateRiemannSolver CLASS //////////
class CentralFluxApproximateRiemannSolver : public BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        CentralFluxApproximateRiemannSolver();						/// Default constructor
        virtual ~CentralFluxApproximateRiemannSolver();					/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
       
        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// MurmanRoeFluxApproximateRiemannSolver CLASS //////////
class MurmanRoeFluxApproximateRiemannSolver : public BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        MurmanRoeFluxApproximateRiemannSolver();					/// Default constructor
        virtual ~MurmanRoeFluxApproximateRiemannSolver();				/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
       
        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// KgpFluxApproximateRiemannSolver CLASS //////////
class KgpFluxApproximateRiemannSolver : public BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        KgpFluxApproximateRiemannSolver();						/// Default constructor
        virtual ~KgpFluxApproximateRiemannSolver();					/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
       
        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// HllApproximateRiemannSolver CLASS //////////
class HllApproximateRiemannSolver : public BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        HllApproximateRiemannSolver();							/// Default constructor
        virtual ~HllApproximateRiemannSolver();						/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// HllcApproximateRiemannSolver CLASS //////////
class HllcApproximateRiemannSolver : public BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        HllcApproximateRiemannSolver();							/// Default constructor
        virtual ~HllcApproximateRiemannSolver();					/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// HllcPlusApproximateRiemannSolver CLASS //////////
class HllcPlusApproximateRiemannSolver : public BaseRiemannSolver {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        HllcPlusApproximateRiemannSolver();					/// Default constructor
        virtual ~HllcPlusApproximateRiemannSolver();				/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        double calculateIntercellFlux(const double &F_L, const double &F_R, const double &U_L, const double &U_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

    protected:

        ////////// PARAMETERS //////////

    private:

};

#endif /*_RHEA_FLOW_SOLVER_*/

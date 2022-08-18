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

////////// CLASS DECLARATION //////////
class FlowSolverRHEA;					/// Flow solver RHEA

//class BaseRiemannSolver;				/// Base Riemann solver
//class DivergenceFluxApproximateRiemannSolver;		/// Divergence scheme approximate Riemann solver
//class MurmanRoeFluxApproximateRiemannSolver;		/// Murman-Roe scheme approximate Riemann solver
//class KgpFluxApproximateRiemannSolver;		/// KGP scheme approximate Riemann solver
//class ShimaFluxApproximateRiemannSolver;		/// SHIMA scheme approximate Riemann solver
//class HllApproximateRiemannSolver;			/// HLL approximate Riemann solver
//class HllcApproximateRiemannSolver;			/// HLLC approximate Riemann solver
//class HllcPlusApproximateRiemannSolver;		/// HLLC+ approximate Riemann solver

class BaseExplicitRungeKuttaMethod;			/// Base explicit Runge-Kutta method
class RungeKutta1Method;				/// Runge-Kutta 1 (RK1) method
class StrongStabilityPreservingRungeKutta2Method;	/// Strong stability preserving Runge-Kutta 2 (SSP-RK2) method
class StrongStabilityPreservingRungeKutta3Method;	/// Strong stability preserving Runge-Kutta 3 (SSP-RK3) method

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
        FlowSolverRHEA(const std::string configuration_file);	/// Parametrized constructor
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

        /// Fill mesh x, y, z, delta_x, delta_y, delta_z fields
        virtual void fillMeshCoordinatesSizesFields();

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

        /// Calculate inviscid fluxes in x-, y- and z-direction
        virtual void calculateInviscidFluxes();

        /// Calculate waves speed
        //virtual void calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R);
        static void calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R);

	/// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
        //virtual double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
        static double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);

        /// Calculate viscous fluxes
        virtual void calculateViscousFluxes();

        /// Advance conserved variables in time
        virtual void timeAdvanceConservedVariables(const int &rk_step);

        /// Advance pressure in time
        virtual void timeAdvancePressure(const int &rk_step);

        /// Output current solver state data
        virtual void outputCurrentStateData();

        /// Update time-averaged quantities
        virtual void updateTimeAveragedQuantities();

        /// Update time mean quantity
        //virtual double updateTimeMeanQuantity(const double &quantity, const double &mean_quantity, const double &delta_t, const double &averaging_time);
        static double updateTimeMeanQuantity(const double &quantity, const double &mean_quantity, const double &delta_t, const double &averaging_time);

        /// Update time rmsf quantity
        //virtual double updateTimeRmsfQuantity(const double &quantity, const double &mean_quantity, const double &rmsf_quantity, const double &delta_t, const double &averaging_time);
        static double updateTimeRmsfQuantity(const double &quantity, const double &mean_quantity, const double &rmsf_quantity, const double &delta_t, const double &averaging_time);

        /// Update time Reynolds-averaged quantity
        //virtual double updateTimeReynoldsAveragedQuantity(const double &quantity_1, const double &mean_quantity_1, const double &quantity_2, const double &mean_quantity_2, const double &reynolds_averaged_quantity, const double &delta_t, const double &averaging_time);
        static double updateTimeReynoldsAveragedQuantity(const double &quantity_1, const double &mean_quantity_1, const double &quantity_2, const double &mean_quantity_2, const double &reynolds_averaged_quantity, const double &delta_t, const double &averaging_time);

        /// Update time Favre-averaged quantity
        //virtual double updateTimeFavreAveragedQuantity(const double &quantity_1, const double &mean_rho_quantity_1, const double &quantity_2, const double &mean_rho_quantity_2, const double &rho, const double &mean_rho, const double &favre_averaged_quantity, const double &delta_t, const double &averaging_time);
        static double updateTimeFavreAveragedQuantity(const double &quantity_1, const double &mean_rho_quantity_1, const double &quantity_2, const double &mean_rho_quantity_2, const double &rho, const double &mean_rho, const double &favre_averaged_quantity, const double &delta_t, const double &averaging_time);

    protected:

        ////////// SOLVER PARAMETERS //////////
	
        /// Fluid & flow properties 
	std::string thermodynamic_model;			/// Thermodynamic model
	std::string transport_coefficients_model;		/// Transport coefficients model

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
	std::string configuration_file;				/// Configuration file name (YAML language)	

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
        int rk_number_stages;					/// Number of stages of Runge-Kutta time discretization method
	std::string riemann_solver_scheme;			/// Riemann solver scheme
	std::string runge_kutta_time_scheme;			/// Runge-Kutta time scheme
        bool transport_pressure_scheme;				/// Activate transport P instead of rhoE

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
	std::string output_data_file_name;			/// Output data file name (HDF5 format)	
        int output_frequency_iter;				/// Data output iteration frequency
        bool generate_xdmf;					/// Generate xdmf file reader
        bool use_restart;					/// Use restart file for initialization
	std::string restart_data_file;				/// Restart data file
        bool time_averaging_active;				/// Activate time averaging
        bool reset_time_averaging;				/// Reset time averaging

        /// Timers information
        bool print_timers;					/// Print timers information
	std::string timers_information_file;			/// Timers information file

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
        DistributedArray P_0_field;				/// 3-D previous field of P

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
        DistributedArray work_vis_rhoe_flux;			/// 3-D viscous fluxes of work rhoe

        /// Source terms
        DistributedArray f_rhou_field;				/// 3-D field of f_rhou
        DistributedArray f_rhov_field;				/// 3-D field of f_rhov
        DistributedArray f_rhow_field;				/// 3-D field of f_rhow
        DistributedArray f_rhoE_field;				/// 3-D field of f_rhoE

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
        DistributedArray avg_mu_field;				/// 3-D field of time-averaged mu
        DistributedArray avg_kappa_field;			/// 3-D field of time-averaged kappa
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
        DistributedArray rmsf_mu_field;				/// 3-D field of root-mean-square-fluctuation mu
        DistributedArray rmsf_kappa_field;			/// 3-D field of root-mean-square-fluctuation kappa
        DistributedArray rmsf_c_v_field;			/// 3-D field of root-mean-square-fluctuation c_v
        DistributedArray rmsf_c_p_field;			/// 3-D field of root-mean-square-fluctuation c_p
        DistributedArray R_reynolds_uu_field;			/// 3-D field of Reynolds-averaged turbulent stress tensor overline{u'u'}
        DistributedArray R_reynolds_uv_field;			/// 3-D field of Reynolds-averaged turbulent stress tensor overline{u'v'}
        DistributedArray R_reynolds_uw_field;			/// 3-D field of Reynolds-averaged turbulent stress tensor overline{u'w'}
        DistributedArray R_reynolds_vv_field;			/// 3-D field of Reynolds-averaged turbulent stress tensor overline{v'v'}
        DistributedArray R_reynolds_vw_field;			/// 3-D field of Reynolds-averaged turbulent stress tensor overline{v'w'}
        DistributedArray R_reynolds_ww_field;			/// 3-D field of Reynolds-averaged turbulent stress tensor overline{w'w'}
        DistributedArray R_favre_uu_field;			/// 3-D field of Favre-averaged turbulent stress tensor overline{rhou''u''}/bar{rho}
        DistributedArray R_favre_uv_field;			/// 3-D field of Favre-averaged turbulent stress tensor overline{rhou''v''}/bar{rho}
        DistributedArray R_favre_uw_field;			/// 3-D field of Favre-averaged turbulent stress tensor overline{rhou''w''}/bar{rho}
        DistributedArray R_favre_vv_field;			/// 3-D field of Favre-averaged turbulent stress tensor overline{rhov''v''}/bar{rho}
        DistributedArray R_favre_vw_field;			/// 3-D field of Favre-averaged turbulent stress tensor overline{rhov''w''}/bar{rho}
        DistributedArray R_favre_ww_field;			/// 3-D field of Favre-averaged turbulent stress tensor overline{rhow''w''}/bar{rho}

	////////// THERMODYNAMIC MODEL, TRANSPORT COEFFICIENTS, RIEMANN SOLVER, EXPLICIT RUNGE-KUTTA METHOD //////////
	////////// COMPUTATIONAL DOMAIN, PARALLEL TOPOLOGY, WRITER/READER, PARALLEL TIMER //////////
        BaseThermodynamicModel *thermodynamics;			/// Thermodynamic model
        BaseTransportCoefficients *transport_coefficients;	/// Transport coefficients
        //BaseRiemannSolver *riemann_solver;			/// Riemann solver
        BaseExplicitRungeKuttaMethod *runge_kutta_method;	/// Runge-Kutta method
        ComputationalDomain *mesh;				/// Computational domain
        ParallelTopology *topo;					/// Parallel topology
        ManagerHDF5 *writer_reader;				/// HDF5 data writer/reader
        ParallelTimer *timers;					/// Parallel timer

	////////// VERSION STUFF //////////
        /// Version numbers consist of three numbers separated by dots; i.e., 1.2.3.
        /// The leftmost number (1) is the major version.
        /// The middle number (2) is the minor version.
        /// The rightmost number (3) is the revision, but it may also refer to a "point release" or "subminor version".
	
        /// Version number (updated 18/05/2022)
	std::string version_number = "2.0.0";			/// Version number	

    private:

};

//////////// BaseRiemannSolver CLASS //////////
//class BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        BaseRiemannSolver();					/// Default constructor
//        virtual ~BaseRiemannSolver();				/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//       
//        /// Calculate waves speed
//        virtual void calculateWavesSpeed(double &S_L, double &S_R, const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R);
//
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        virtual double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type) = 0;
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// DivergenceFluxApproximateRiemannSolver CLASS //////////
//class DivergenceFluxApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        DivergenceFluxApproximateRiemannSolver();						/// Default constructor
//        virtual ~DivergenceFluxApproximateRiemannSolver();					/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//       
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// MurmanRoeFluxApproximateRiemannSolver CLASS //////////
//class MurmanRoeFluxApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        MurmanRoeFluxApproximateRiemannSolver();					/// Default constructor
//        virtual ~MurmanRoeFluxApproximateRiemannSolver();				/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//       
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// KgpFluxApproximateRiemannSolver CLASS //////////
//class KgpFluxApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        KgpFluxApproximateRiemannSolver();						/// Default constructor
//        virtual ~KgpFluxApproximateRiemannSolver();					/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//       
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// ShimaFluxApproximateRiemannSolver CLASS //////////
//class ShimaFluxApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        ShimaFluxApproximateRiemannSolver();						/// Default constructor
//        virtual ~ShimaFluxApproximateRiemannSolver();					/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//       
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// HllApproximateRiemannSolver CLASS //////////
//class HllApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        HllApproximateRiemannSolver();							/// Default constructor
//        virtual ~HllApproximateRiemannSolver();						/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//        
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// HllcApproximateRiemannSolver CLASS //////////
//class HllcApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        HllcApproximateRiemannSolver();							/// Default constructor
//        virtual ~HllcApproximateRiemannSolver();					/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//        
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};
//
//////////// HllcPlusApproximateRiemannSolver CLASS //////////
//class HllcPlusApproximateRiemannSolver : public BaseRiemannSolver {
//   
//    public:
//
//        ////////// CONSTRUCTORS & DESTRUCTOR //////////
//        HllcPlusApproximateRiemannSolver();					/// Default constructor
//        virtual ~HllcPlusApproximateRiemannSolver();				/// Destructor
//
//	////////// GET FUNCTIONS //////////
//
//	////////// SET FUNCTIONS //////////
//
//	////////// METHODS //////////
//        
//        /// Calculate intercell flux ... var_type corresponds to: 0 for rho, 1-3 for rhouvw, 4 for rhoE
//        double calculateIntercellFlux(const double &rho_L, const double &rho_R, const double &u_L, const double &u_R, const double &v_L, const double &v_R, const double &w_L, const double &w_R, const double &E_L, const double &E_R, const double &P_L, const double &P_R, const double &a_L, const double &a_R, const int &var_type);
//
//    protected:
//
//        ////////// PARAMETERS //////////
//
//    private:
//
//};

////////// BaseExplicitRungeKuttaMethod CLASS //////////
class BaseExplicitRungeKuttaMethod {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        BaseExplicitRungeKuttaMethod();				/// Default constructor
        virtual ~BaseExplicitRungeKuttaMethod();		/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
       
        /// Set stage coefficients: y_rk+1 = rk_a*y_n + rk_b*y_rk + rk_c*delta_t*f( y_rk )
        virtual void setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage) = 0;

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// RungeKutta1Method CLASS //////////
class RungeKutta1Method : public BaseExplicitRungeKuttaMethod {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        RungeKutta1Method();				/// Default constructor
        virtual ~RungeKutta1Method();			/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Set stage coefficients: y_rk+1 = rk_a*y_n + rk_b*y_rk + rk_c*delta_t*f( y_rk )
        void setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// StrongStabilityPreservingRungeKutta2Method CLASS //////////
class StrongStabilityPreservingRungeKutta2Method : public BaseExplicitRungeKuttaMethod {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        StrongStabilityPreservingRungeKutta2Method();				/// Default constructor
        virtual ~StrongStabilityPreservingRungeKutta2Method();			/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Set stage coefficients: y_rk+1 = rk_a*y_n + rk_b*y_rk + rk_c*delta_t*f( y_rk )
        void setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage);

    protected:

        ////////// PARAMETERS //////////

    private:

};

////////// StrongStabilityPreservingRungeKutta3Method CLASS //////////
class StrongStabilityPreservingRungeKutta3Method : public BaseExplicitRungeKuttaMethod {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        StrongStabilityPreservingRungeKutta3Method();				/// Default constructor
        virtual ~StrongStabilityPreservingRungeKutta3Method();			/// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////
        
        /// Set stage coefficients: y_rk+1 = rk_a*y_n + rk_b*y_rk + rk_c*delta_t*f( y_rk )
        void setStageCoefficients(double &rk_a, double &rk_b, double &rk_c, const int &rk_time_stage);

    protected:

        ////////// PARAMETERS //////////

    private:

};

#endif /*_RHEA_FLOW_SOLVER_*/

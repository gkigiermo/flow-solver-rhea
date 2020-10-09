#ifndef _MY_RHEA_
#define _MY_RHEA_

////////// INCLUDES //////////
#include "src/FlowSolverRHEA.hpp"

////////// NAMESPACES //////////
using namespace std;

////////// myRHEA CLASS //////////
class myRHEA : public FlowSolverRHEA {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        myRHEA(const string configuration_file) : FlowSolverRHEA(configuration_file) {};	/// Parametrized constructor
        virtual ~myRHEA() {};									/// Destructor

	////////// SOLVER METHODS //////////
        
        /// Set initial conditions: u, v, w, P and T ... needs to be modified/overwritten according to the problem under consideration
        void setInitialConditions();

        /// Calculate rhou, rhov, rhow and rhoE source terms ... needs to be modified/overwritten according to the problem under consideration
        void calculateSourceTerms();

        /// Calculate inviscid fluxes ... needs to be modified/overwritten for this problem
        void calculateInviscidFluxes();

    protected:

    private:

};

#endif /*_MY_RHEA_*/

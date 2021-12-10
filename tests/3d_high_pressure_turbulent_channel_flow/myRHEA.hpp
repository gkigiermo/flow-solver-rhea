#ifndef _MY_RHEA_
#define _MY_RHEA_

////////// INCLUDES //////////
#include "../flowsolverrhea/src/FlowSolverRHEA.hpp"

////////// NAMESPACES //////////
using namespace std;

////////// myRHEA CLASS //////////
class myRHEA : public FlowSolverRHEA {
   
    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        myRHEA(const string configuration_file) : FlowSolverRHEA(configuration_file) {};	/// Parametrized constructor
        virtual ~myRHEA() {};									/// Destructor

	////////// SOLVER METHODS //////////
       
	/// Execute (aggregated method) RHEA
        void execute();
 
        /// Set initial conditions: u, v, w, P and T ... needs to be modified/overwritten according to the problem under consideration
        void setInitialConditions();

        /// Calculate rhou, rhov, rhow and rhoE source terms ... needs to be modified/overwritten according to the problem under consideration
        void calculateSourceTerms();

    protected:

    private:

};

#endif /*_MY_RHEA_*/

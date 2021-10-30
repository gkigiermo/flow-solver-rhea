#ifndef _ROOT_FINDING_MINIMIZATION_HPP_
#define _ROOT_FINDING_MINIMIZATION_HPP_

////////// INCLUDES //////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

////////// NAMESPACES //////////
using namespace std;

////////// CLASS DECLARATION //////////
class BaseRootFindingMinimization;			/// BaseRootFindingMinimization
class NewtonRaphson;					/// NewtonRaphson

////////// FUNCTION DECLARATION //////////

////////// BaseRootFindingMinimization CLASS //////////
class BaseRootFindingMinimization {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        BaseRootFindingMinimization();					/// Default constructor
        virtual ~BaseRootFindingMinimization();                         /// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////

	/// Numerical Recipes in C++ SQR algorithm
        inline double SQR(const double &a) { return a*a; };	

        /// Runs the minimization algorithm, calculating the minimum (fxmin) and its independent variables (xmin)
        virtual void solve(double &fxmin, vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance) = 0; 
        
        /// Evaluates the functions (residuals) in position x
        virtual void function_vector(vector<double> &x, vector<double> &fx) {}; 

    protected:

        ////////// PARAMETERS //////////

    private:

};



////////// NewtonRaphson CLASS //////////

class NewtonRaphson : public BaseRootFindingMinimization {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        //NewtonRaphson();						/// Default constructor 
        NewtonRaphson(vector<double> &fvec_);				/// Parametrized constructor
        virtual ~NewtonRaphson();               			/// Destructor    
 
	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////

        /// Runs the minimization algorithm, calculating the minimum (fxmin) and its independent variables (xmin)
        void solve(double &fxmin, vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance);

        /// Performs minimization linear search
        void lnsrch(vector<double> &xold, vector<double> &p, vector<double> &x, vector< vector<double> > &fjac, vector<int> &indx);

        /// Numerical recipes in C++ LU decomposition method
	void lubksb(vector< vector<double> > &a, vector<int> &indx, vector<double> &b);

        /// Numerical recipes in C++ LU decomposition method
        void ludcmp(vector< vector<double> > &a, vector<int> &indx, double &d);

        /// Numerical recipes in C++ Jacobian method
        void fdjac(vector<double> &x, vector< vector<double> > &df);

    protected:

        ////////// PARAMETERS //////////

	/// Vector of functions to be zeroed
	vector<double> &fvec;

};

#endif /*_ROOT_FINDING_MINIMIZATION_HPP_*/

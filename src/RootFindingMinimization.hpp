#ifndef _ROOT_FINDING_MINIMIZATION_HPP_
#define _ROOT_FINDING_MINIMIZATION_HPP_

////////// INCLUDES //////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>

////////// CLASS DECLARATION //////////
class BaseRootFindingMinimization;			/// BaseRootFindingMinimization
class Brent;						/// Brent
class NewtonRaphson;					/// NewtonRaphson
class Broyden;						/// Broyden

////////// FUNCTION DECLARATION //////////


////////// BaseRootFindingMinimization CLASS //////////
class BaseRootFindingMinimization {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        //BaseRootFindingMinimization();				/// Default constructor
        BaseRootFindingMinimization(std::vector<double> &fvec_);	/// Parametrized constructor
        virtual ~BaseRootFindingMinimization();                         /// Destructor

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////

        /// Numerical Recipes in C++ shift algorithm
        inline void shft3(double &a, double &b, double &c, const double &d) { a = b; b = c; c = d; };

	// Numerical Recipes in C++ sign algorithm
        inline double SIGN(const double &a, const double &b) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); };

	/// Numerical Recipes in C++ SQR algorithm
        inline double SQR(const double &a) { return a*a; };	

        /// Numerical recipes in C++ linear search
        void lnsrch(std::vector<double> &xold, double &fold, std::vector<double> &g, std::vector<double> &p, std::vector<double> &x,double &f, double &stpmax, bool &check,std::vector< std::vector<double> > &fjac);

        /// Numerical recipes in C++ LU decomposition method
	void lubksb(std::vector< std::vector<double> > &a, std::vector<int> &indx, std::vector<double> &b);

        /// Numerical recipes in C++ LU decomposition method
        void ludcmp(std::vector< std::vector<double> > &a, std::vector<int> &indx, double &d);

        /// Numerical recipes in C++ Jacobian method
        void fdjac(std::vector<double> &x, std::vector< std::vector<double> > &df);
        
        /// Numerical recipes in C++ QR decomposition method
        void qrdcmp(std::vector< std::vector<double> > &a, std::vector<double> &c, std::vector<double> &d, bool &sing);
        
        /// Numerical recipes in C++ rsolv method
        void rsolv(std::vector< std::vector<double> > &a, std::vector<double> &d, std::vector<double> &b);
        
        /// Numerical recipes in C++ qrupdt method
        void qrupdt(std::vector< std::vector<double> > &r, std::vector< std::vector<double> > &qt, std::vector<double> &u, std::vector<double> &v);
        
        /// Numerical recipes in C++ rotate method
        void rotate(std::vector< std::vector<double> > &r, std::vector< std::vector<double> > &qt, const int &i, const double &a, const double &b);

        /// Numerical recipes in C++ fmin method
        double fmin(std::vector<double> &x);

        /// Runs the minimization algorithm, calculating the minimum (fxmin) and its independent variables (xmin)
        virtual void solve(double &fxmin, std::vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance) = 0; 
        
        /// Evaluates the functions (residuals) in position x
        virtual void function_vector(std::vector<double> &x, std::vector<double> &fx) {}; 

    protected:

        ////////// PARAMETERS //////////

	/// Vector of functions to be zeroed
	std::vector<double> &fvec;

    private:

};


////////// Brent CLASS //////////

class Brent : public BaseRootFindingMinimization {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        //Brent();							/// Default constructor 
        Brent(std::vector<double> &fvec_);				/// Parametrized constructor
        virtual ~Brent();               				/// Destructor    

	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////
        void set_ax_bx_cx(const double &ax_, const double &bx_, const double &cx_);

	////////// METHODS //////////

        /// Runs the minimization algorithm, calculating the minimum (fxmin) and its independent variables (xmin)
        void solve(double &fxmin, std::vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance);

    protected:

        ////////// PARAMETERS //////////

        /// Bracketing triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx))
        double ax, bx, cx;

};


////////// NewtonRaphson CLASS //////////

class NewtonRaphson : public BaseRootFindingMinimization {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        //NewtonRaphson();						/// Default constructor 
        NewtonRaphson(std::vector<double> &fvec_);			/// Parametrized constructor
        virtual ~NewtonRaphson();               			/// Destructor    
 
	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////

        /// Runs the minimization algorithm, calculating the minimum (fxmin) and its independent variables (xmin)
        void solve(double &fxmin, std::vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance);

    protected:

        ////////// PARAMETERS //////////

};


////////// Broyden CLASS //////////

class Broyden : public BaseRootFindingMinimization {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        //Broyden();							/// Default constructor 
        Broyden(std::vector<double> &fvec_);				/// Parametrized constructor
        virtual ~Broyden();               				/// Destructor    
 
	////////// GET FUNCTIONS //////////

	////////// SET FUNCTIONS //////////

	////////// METHODS //////////

        /// Runs the minimization algorithm, calculating the minimum (fxmin) and its independent variables (xmin)
        void solve(double &fxmin, std::vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance);

    protected:

        ////////// PARAMETERS //////////

};

#endif /*_ROOT_FINDING_MINIMIZATION_HPP_*/

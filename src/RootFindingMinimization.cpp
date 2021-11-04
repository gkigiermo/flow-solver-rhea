#include "RootFindingMinimization.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
#define _TAYLORED_FDJAC_ 0
#define _ACTIVATE_COUT_ 0


////////// BaseRootFindingMinimization CLASS //////////

BaseRootFindingMinimization::BaseRootFindingMinimization() {};
        
BaseRootFindingMinimization::~BaseRootFindingMinimization() {};


////////// NewtonRaphson CLASS //////////

//NewtonRaphson::NewtonRaphson() : BaseRootFindingMinimization() {};

NewtonRaphson::NewtonRaphson(vector<double> &fvec_) : BaseRootFindingMinimization(), fvec(fvec_) {};

NewtonRaphson::~NewtonRaphson() {};                                                                           

void NewtonRaphson::solve(double &fxmin, vector<double> &xmin, const int &max_iter, int &iter, const double &tolerance) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    //const int MAXITS=200;
    const int MAXITS=max_iter;
    //const double EPS=numeric_limits<double>::epsilon();
    const double EPS=tolerance;
    //const double TOLF=1.0e-8, TOLX=EPS, STPMX=100.0, TOLMIN=1.0e-12;
    const double TOLF=EPS, TOLX=EPS;
    int i;
    //double d,f,fold,stpmax,sum,temp,test;
    double d,sum,temp,test;

    int n = xmin.size();
    vector<int> indx(n);
    vector<double> p(n),xold(n);
    vector < vector < double > > fjac(n,vector<double> (n));

    function_vector(xmin,fvec);
    test = 0.0;
    for(i = 0; i < n; i++) {
      if( fabs( fvec[i] ) > test ) fxmin = test = fabs( fvec[i] );
    }
    if( test < 0.01*TOLF ) {
      return;
    }
    sum = 0.0;
    for(i = 0; i < n; i++) sum += SQR( xmin[i] );
    //stpmax = STPMX*max( sqrt( sum ), double( n ) );
    for(iter = 0; iter < MAXITS; iter++) {
      fdjac(xmin,fjac);
      for(i = 0; i < n; i++) xold[i] = xmin[i];
      //fold = f;
      for(i = 0; i < n; i++) p[i] = -fvec[i];
      ludcmp(fjac,indx,d);
      lubksb(fjac,indx,p);
      lnsrch(xold,p,xmin,fjac,indx);
      test = 0.0;
      for(i = 0; i < n; i++) {
        if( fabs( fvec[i] ) > test ) fxmin = test = fabs( fvec[i] );
      }
      if( test < TOLF ) {
#if _ACTIVATE_COUT_
        cout << "Newton-Raphson's minimization TOLF error: " << test << " ( iteration: " << iter << " )" << endl;
#endif
	return;
      }
      test = 0.0;
      for(i = 0; i < n; i++) {
	temp = ( fabs( xmin[i] - xold[i] ) )/max( fabs(xmin[i]), 1.0 );
	if( temp > test ) fxmin = test = temp;
      }
      if( test < TOLX ) {
#if _ACTIVATE_COUT_
        cout << "Newton-Raphson's minimization TOLX error: " << test << " ( iteration: " << iter << " )" << endl;
#endif
	return;
      }
    }
#if _ACTIVATE_COUT_
    cout << "MAXITS exceeded in Newton-Raphson" << endl;
#endif

};

void NewtonRaphson::lnsrch(vector<double> &xold,vector<double> &p,vector<double> &x,vector< vector<double> > &fjac, vector<int> &indx) {

    /// Taylored linear search:

    const double sigma   = 0.01;
    const double Tau     = 0.1;
    //const double alamMin = 1e-6;
    //const double alamMin = 5.0e-2;
    const double alamMin = 1.0e-1;
    const int max_iter   = 10;

    int i, j, n = xold.size();

    double f0 = 0.0;
    for(i = 0; i < n; ++i) f0 += p[i]*p[i];
    f0 = sqrt(f0/n);
    f0 *= 0.5*f0;

    vector<double> myp(n);
    for(i = 0; i < n; ++i) myp[i] = p[i];

    double alam = 1.0;
    double f, fMax, tmp;
    for( j = 0; j < max_iter; j++ ) {
       /// Compute f( x + Lambda*dx )
       for(i = 0; i < n; ++i) x[i] = xold[i] + alam*p[i];
       for(i = 0; i < n; i++) myp[i] = -fvec[i];
       lubksb(fjac,indx,myp);
       f = 0.0;
       for(i = 0; i < n; ++i) f += myp[i]*myp[i];
       f = sqrt(f/n);
       f *= 0.5*f;

       fMax = ( 1.0 - 2.0*alam*sigma )*f0;
       if( f > fMax ) {
          tmp = alam*alam*f0/( ( 2.0*alam - 1.0 )*f0 + f );
          alam = ( Tau*alam > tmp ) ? Tau*alam : tmp;
          if( alam < alamMin ) { 
             alam = 10.0*alamMin;
             for(i = 0; i < n; ++i) x[i] = xold[i] + alam*p[i];
             break;
          }
       } else {
          // Lamdba is good
          //for(int i = 0; i < n; ++i) x[i] = xold[i] + alam*p[i];
          break;
       }
    }

};

void NewtonRaphson::lubksb(vector< vector<double> > &a, vector<int> &indx, vector<double> &b) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    int i,ii=0,ip,j;
    double sum;

    int n=a.size();
    for(i=0;i<n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if(ii != 0) {
	for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
      } else if (sum != 0.0) {
	ii=i+1;
      }
      b[i]=sum;
    }
    for(i=n-1;i>=0;i--) {
      sum=b[i];
      for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
    }

};

void NewtonRaphson::ludcmp(vector< vector<double> > &a, vector<int> &indx, double &d) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    const double TINY=1.0e-20;
    int i,j,k,imax=0;
    double big,dum,sum,temp;

    int n=a.size();
    vector<double> vv(n);
    d=1.0;
    for(i=0;i<n;i++) {
      big=0.0;
      for(j=0;j<n;j++) {
	if ((temp=fabs(a[i][j])) > big) big=temp;
      }
#if _ACTIVATE_COUT_
      if(big == 0.0) cout << "Singular matrix in routine ludcmp" <<endl;
#endif
      vv[i]=1.0/big;
    }
    for(j=0;j<n;j++) {
      for(i=0;i<j;i++) {
	sum=a[i][j];
	for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
      }
      big=0.0;
      for(i=j;i<n;i++) {
	sum=a[i][j];
	for(k=0;k<j;k++) sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
	if((dum=vv[i]*fabs(sum)) >= big) {
	  big=dum;
	  imax=i;
	}
      }
      if(j != imax) {
	for(k=0;k<n;k++) {
	  dum=a[imax][k];
	  a[imax][k]=a[j][k];
	  a[j][k]=dum;
	}
	d = -d;
	vv[imax]=vv[j];
      }
      indx[j]=imax;
      if(a[j][j] == 0.0) a[j][j]=TINY;
      if(j != n-1) {
	dum=1.0/(a[j][j]);
	for(i=j+1;i<n;i++) a[i][j] *= dum;
      }
    }

};

void NewtonRaphson::fdjac(vector<double> &x, vector< vector<double> > &df) {

    /// Newton-Raphson method using approximated derivatives:
    /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
    /// Numerical recipes in C++.
    /// Cambridge University Press, 2001.

    // Update fvec
    function_vector(x,fvec);

#if _TAYLORED_FDJAC_
    const double EPS   = 1.0e-6;
    const double DELTA = 1.0e-14;
    int i,j;
    //double h,h_p,h_m,temp;
    double h,temp;

    int n=x.size();
    vector<double> f(n);
    //vector<double> f_p(n);
    //vector<double> f_m(n);
    for(j=0;j<n;j++) {
      temp=x[j];
      h=EPS*fabs(temp) + DELTA;
      x[j]=temp+h;
      h=x[j]-temp;
      //h_p=x[j]-temp;
      function_vector(x,f);
      //function_vector(x,f_p);
      //x[j]=temp-h;
      //h_m=temp-x[j];
      //function_vector(x,f_m);
      x[j]=temp;
      for(i=0;i<n;i++) {
        df[i][j]=(f[i]-fvec[i])/h;
        //df[i][j]=(f_p[i]-f_m[i])/(h_p+h_m);
      }
    }
#else
    const double EPS   = 1.0e-8;
    const double DELTA = 1.0e-14;
    int i,j;
    double h,temp;

    int n = x.size();
    vector<double> f(n);

    for(j = 0; j < n; j++) {
      temp = x[j];
      //h = EPS*fabs( temp );
      //if(h == 0.0) h = EPS;
      h = EPS*fabs( temp ) + DELTA;
      x[j] = temp + h;
      //h = x[j] - temp;
      function_vector(x,f);
      x[j] = temp;
      for(i = 0; i < n; i++) {
        df[i][j] = ( f[i] - fvec[i] )/h;
      }
    }   
#endif

};
